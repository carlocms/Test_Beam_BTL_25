#include "interface/TOFHIRThresholdZero.h"
#include "interface/AnalysisUtils.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "interface/Na22SpectrumAnalyzer.h"
#include "interface/Na22SpectrumAnalyzerSingleBar.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <algorithm>
#include <iterator>
#include <dirent.h>

#include "TFile.h"
#include "TChain.h"
#include "TH1F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraphErrors.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TRandom3.h"


// Versione dello step 1 per uso Amplitude walk correction studies with MCP informatio as the external Time reference.
// Author: Carlo Giraldin - INFN Trieste - February 2026




int main(int argc, char** argv)
{
  setTDRStyle();
  float cpu[2]{0}, mem[2]={0}, vsz[2]={0}, rss[2]={0};  
  

  if( argc < 2 )
  {
    std::cout << ">>> moduleCharacterization_step1::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }
  
  
  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);
  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);
  

  //--- open files and make the tree chain
  std::string inputDir = opts.GetOpt<std::string>("Input.inputDir_wMCP");
  std::string fileBaseName = opts.GetOpt<std::string>("Input.fileBaseName");
  std::string runs = opts.GetOpt<std::string>("Input.runs");
  int maxEntries = opts.GetOpt<int>("Input.maxEntries");
  int usePedestals = opts.GetOpt<int>("Input.usePedestals");
  std::string source = opts.GetOpt<std::string>("Input.sourceName");
  int useTrackInfo = opts.GetOpt<int>("Input.useTrackInfo");
  float my_step1 = opts.GetOpt<float>("Input.vov");
  bool Debug_on = opts.GetOpt<bool>("Input.debug_mod");
  
  int DUTasic = opts.GetOpt<int>("Channels.DUTasic");
  int REFasic = opts.GetOpt<int>("Channels.REFasic");

  std::string discCalibrationFile = opts.GetOpt<std::string>("Input.discCalibration");
  TOFHIRThresholdZero thrZero(discCalibrationFile,1);

  TChain* tree = new TChain("data","data");
  TChain* tree_MCP = new TChain("MCP","MCP");
  


  std::stringstream ss(runs); 
  std::string token;
  while( std::getline(ss,token,',') )
  {
    std::stringstream ss2(token);
    std::string token2;
    int runMin = -1;
    int runMax = -1;
    while( std::getline(ss2,token2,'-') )
    {
      if( runMin != -1 && runMax == -1 ) runMax = atoi(token2.c_str());
      if( runMin == -1 ) runMin = atoi(token2.c_str());
    }
    if( runMax == -1 ) runMax = runMin;
    
    for(int run = runMin; run <= runMax; ++run) {

      // -- analyze only spills at a chosen OV to speed up analysis
      // - list of files in run folder
      DIR *dir_ptr;
      struct dirent *diread;
      std::vector<std::string> filenames;
      std::string directory_path = Form("%s/%04d/",inputDir.c_str(),run);
      std::cout << directory_path.c_str()<<std::endl;

      if ((dir_ptr = opendir(directory_path.c_str())) != nullptr) {
        while ((diread = readdir(dir_ptr)) != nullptr) {
	  //std::cout << diread->d_name << std::endl;
	  std::string fname(diread->d_name);
	  filenames.push_back(fname);
	}
	closedir(dir_ptr);    
      }	

      for (auto fname: filenames) {
	std::cout << fname.c_str() << std::endl;;

	if (fname == ".") continue;
	if (fname == "..") continue;

	// -- check if Vov selected
	bool addFile = true;
	TFile *f = TFile::Open((directory_path+fname).c_str());
	TTree *tmpTree = f->Get<TTree>("data");
	float step1;
	tmpTree->SetBranchAddress("step1",&step1);
	tmpTree->GetEntry(0);
	if (my_step1 > 0 && step1!=my_step1) addFile = false;
	delete tmpTree;
	f->Close();
	
	if (addFile){
	  std::cout << ">>> step1 = " << step1 << " --> Adding file: " << fname.c_str()<< std::endl;
	  tree->Add((directory_path+fname).c_str());
	  tree_MCP->Add((directory_path+fname).c_str());
	}
      }

      
    }
  }
  
     
  //--- define channels (read mapping from the configuration file)
  std::vector<unsigned int> channelMapping = opts.GetOpt<std::vector<unsigned int> >("Channels.channelMapping");
  
  int chL[16];
  int chR[16];
  for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){

      chL[iBar] = channelMapping[iBar*2+0]+32*DUTasic;
      chR[iBar] = channelMapping[iBar*2+1]+32*DUTasic;
      
    std::cout << "Bar: " << iBar << "   chL: "<< chL[iBar] << "    chR: " <<chR[iBar] <<std::endl;
  }
  

  //--- define branches from TTree data:
  float step1, step2;
  int channelIdx[256];
  std::vector<unsigned short> *qfine = 0;
  std::vector<float> *tot = 0;
  std::vector<float> *energy = 0;
  std::vector<long long> *time = 0;
  std::vector<float>* qT1 = 0;
  std::vector<unsigned short>* t1fine = 0;
  std::vector<unsigned short> *t1coarse = 0;
  
  int nhits;
  float x, y;
  
  tree -> SetBranchStatus("*",0);
  tree -> SetBranchStatus("step1",  1); tree -> SetBranchAddress("step1",  &step1);
  tree -> SetBranchStatus("step2",  1); tree -> SetBranchAddress("step2",  &step2);
  
  tree -> SetBranchStatus("channelIdx",  1); tree -> SetBranchAddress("channelIdx",  channelIdx);
  
  tree -> SetBranchStatus("qfine",  1); tree -> SetBranchAddress("qfine",   &qfine);  
  tree -> SetBranchStatus("tot",    1); tree -> SetBranchAddress("tot",       &tot);
  tree -> SetBranchStatus("energy", 1); tree -> SetBranchAddress("energy", &energy);
  tree -> SetBranchStatus("time",   1); tree -> SetBranchAddress("time",     &time);
  
  tree -> SetBranchStatus("qT1",       1); tree -> SetBranchAddress("qT1",      &qT1);
  tree -> SetBranchStatus("t1fine",    1); tree -> SetBranchAddress("t1fine",   &t1fine);
  tree -> SetBranchStatus("t1coarse",    1); tree -> SetBranchAddress("t1coarse",   &t1coarse);  

  if ( !opts.GetOpt<std::string>("Input.sourceName").compare("TB") &&  useTrackInfo ){
    tree -> SetBranchStatus("nhits_WC", 1); tree -> SetBranchAddress("nhits_WC",  &nhits);
    tree -> SetBranchStatus("x_WC", 1);     tree -> SetBranchAddress("x_WC",          &x);
    tree -> SetBranchStatus("y_WC", 1);     tree -> SetBranchAddress("y_WC",          &y);
  }
  

  //--- define branches from TTree MCP:
  long long mcp_index;
  double mcp_peak_time;
  double mcp_peak_amp;
  double mcp_phi_peak;
  double mcp_phi_peak_from_edge;
  double mcp_trigger_time;
  double mcp_trigger_offset_ps;
  double mcp_phi_trigger;
  double mcp_phi_trigger_from_edge;
  double mcp_t0_abs_ps;
	
  tree_MCP->SetBranchAddress("index",&mcp_index);
  tree_MCP->SetBranchAddress("peak_time",&mcp_peak_time);
  tree_MCP->SetBranchAddress("peak_amp",&mcp_peak_amp);
  tree_MCP->SetBranchAddress("phi_peak",&mcp_phi_peak);
  tree_MCP->SetBranchAddress("phi_peak_from_edge",&mcp_phi_peak_from_edge);
  tree_MCP->SetBranchAddress("trigger_time",&mcp_trigger_time);
  tree_MCP->SetBranchAddress("trigger_offset_ps",&mcp_trigger_offset_ps);
  tree_MCP->SetBranchAddress("phi_trigger",&mcp_phi_trigger);
  tree_MCP->SetBranchAddress("phi_trigger_from_edge",&mcp_phi_trigger_from_edge);
  tree_MCP->SetBranchAddress("t0_abs_ps",&mcp_t0_abs_ps);



  //--- get plot settings
  std::vector<float> Vov = opts.GetOpt<std::vector<float> >("Plots.Vov");
  std::vector<int> energyBins = opts.GetOpt<std::vector<int> >("Plots.energyBins");
  std::vector<int> energyMins = opts.GetOpt<std::vector<int> >("Plots.energyMins");
  std::vector<int> energyMaxs = opts.GetOpt<std::vector<int> >("Plots.energyMaxs");
  
  std::map<float,int> map_energyBins;
  std::map<float,int> map_energyMins;
  std::map<float,int> map_energyMaxs;
  for(unsigned int ii = 0; ii < Vov.size(); ++ii)
    {
      map_energyBins[Vov[ii]] = energyBins[ii];
      map_energyMins[Vov[ii]] = energyMins[ii];
      map_energyMaxs[Vov[ii]] = energyMaxs[ii];
    }
  
  //float vetoEnergyThreshold = opts.GetOpt<float>("Plots.vetoEnergyThreshold");   


  // -- read minimum energy for each bar from file
  int vetoOtherBars = opts.GetOpt<int>("Cuts.vetoOtherBars");
  std::string minEnergiesFileName = opts.GetOpt<std::string>("Cuts.minEnergiesFileName");
  std::map < std::pair<int, float>, float> minE; 
  if (minEnergiesFileName != "") {
    std::ifstream minEnergiesFile;
    minEnergiesFile.open(minEnergiesFileName);
    std::string line;
    int bar;
    float ov;
    float value;
    while ( minEnergiesFile.good() ){
      getline(minEnergiesFile, line);
      std::istringstream ss(line);
      ss >> bar >> ov >> value; 
      minE[std::make_pair(bar,ov)] = value; 
      std::cout<< bar <<  "   " << ov << "  " << minE[std::make_pair(bar,ov)] <<std::endl;
    }
  }
  else{
    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){
      for(unsigned int ii = 0; ii < Vov.size(); ++ii){
	minE[std::make_pair(iBar, Vov[ii])] = map_energyMins[Vov[ii]];
      }
    }
  }
  

  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameStep1_withMCP");
  TFile* outFile = TFile::Open(Form("%s",outFileName.c_str()),"RECREATE");
  outFile -> cd();
  
  std::map<int,TTree*> outTrees;

  std::map<int,TH1F*> h1_qfineL;
  std::map<int,TH1F*> h1_qfineR;  
  std::map<int,TH1F*> h1_totL;
  std::map<int,TH1F*> h1_totR;
  std::map<int,TH1F*> h1_energyL;
  std::map<int,TH1F*> h1_energyR;
  std::map<int,TH1F*> h1_energyLR;
  std::map<int,TH1F*> h1_energyLR_ext;
  std::map<int,TH1F*> h1_MeanTimeREF;
  std::map<int,TH1F*> h1_maxBar_DUT;
  std::map<int,TH1F*> h1_maxBar_REF;

  std::map<int,TCanvas*> c;
  std::map<int,std::vector<float>*> rangesLR;
  std::map<int,bool> acceptEvent;





  // -- Coincidence pre loop
  if( !opts.GetOpt<std::string>("Coincidence.status").compare("yes") &&
      ( !opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("Na22") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("TB") ||
	!opts.GetOpt<std::string>("Input.sourceName").compare("keepAll") ) )
    {
      float energyL_ext;
      float energyR_ext;
      int chL_ext = opts.GetOpt<float>("Coincidence.chL");//NB: gli passo direttamente da cfg il ch barra 8 +64
      int chR_ext = opts.GetOpt<float>("Coincidence.chR");//NB: gli passo direttamente da cfg il ch barra 8 +64
      
      int nEntries = tree->GetEntries();
      if( maxEntries > 0 ) nEntries = maxEntries;
      for(int entry = 0; entry < nEntries; ++entry){
	tree -> GetEntry(entry);
	acceptEvent[entry] = false;
	if( entry%200000 == 0 ){
	  std::cout << "\n>>> external bar loop: reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << std::endl;
	  TrackProcess(cpu, mem, vsz, rss);
	}
      
	float Vov = step1;
        float vth1 = float(int(step2/10000)-1);
        float vth2 = float(int((step2-10000*(vth1+1))/100.)-1);
        float vthE = float(int((step2-10000*(vth1+1)) - 100*(vth2+1))-1);
        float vth = 0.;
	
	if( entry%20000 == 0 ){
	  std::cout << step2 << " ith1: " << vth1 << " ith2: " << vth2 << " E: " << vthE << std::endl;
	}
	if(!opts.GetOpt<std::string>("Input.vth").compare("vth1"))  { vth = vth1;}
	if(!opts.GetOpt<std::string>("Input.vth").compare("vth2"))  { vth = vth2;}
	if(!opts.GetOpt<std::string>("Input.vth").compare("vthE"))  { vth = vthE;}

	// select only one OV
	if (my_step1 > 0  && my_step1 != step1) continue;

	if (channelIdx[chL_ext] <0 || channelIdx[chR_ext] <0) continue; //Condizione che impone che la barra 8 REF deve essere accesa
	

	// --- calculate energy sum/Nbars for module - useful to remove showering events/cross talk
	if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB")){

	  auto it = std::find( channelMapping.begin(), channelMapping.end(), chL_ext);       	
	  int thisBar = 0;
	  if (it != channelMapping.end()) thisBar = int((it - channelMapping.begin()))/2;

	  float energySumArray = 0.;
	  int   nActiveBarsArray = 0;
	  
	  // check adjacent bars
	  /*for( int iBar = thisBar - 2; iBar < thisBar+3; ++iBar) {                                                                                                                            
	    if (iBar<0 || iBar>15 ) continue;                                                                                                                                                   
	    if (iBar == thisBar) continue;
	    //int chL_iext = channelMapping[iBar*2+0];// array0 for coincidence is hard coded... - to be fixed
	    //int chR_iext = channelMapping[iBar*2+1];// array0 for coincidence is hard coded... - to be fixed
	    int chL_iext = channelMapping[iBar*2+0]+64;// array1 for coincidence is hard coded... - to be fixed
	    int chR_iext = channelMapping[iBar*2+1]+64;// array1 for coincidence is hard coded... - to be fixed
	    float energyL_iext = (*energy)[channelIdx[chL_iext]];              
	    float energyR_iext = (*energy)[channelIdx[chR_iext]]; 
	    float totL_iext    = 0.001*(*tot)[channelIdx[chL_iext]];              
	    float totR_iext    = 0.001*(*tot)[channelIdx[chR_iext]]; 
	    if ( totL_iext > 0 && totL_iext < 100 && totR_iext > 0 && totR_iext < 100   ){
	      float energyMean=(energyL_iext+energyR_iext)/2;
	      if (energyMean>200 && energyMean < 1024){
		energySumArray+=energyMean;
		nActiveBarsArray+=1;
	      }
	    }
	  }
	  if (nActiveBarsArray > 0 ) continue; // veto signals in the adjacent bars
	  */


	  for(int iBar = 0; iBar < int(channelMapping.size())/2; ++iBar) {
	    int chL_iext = channelMapping[iBar*2+0] + REFasic*32;
	    int chR_iext = channelMapping[iBar*2+1] + REFasic*32;

	   // int chL_iext = channelMapping[iBar*2+0];// module under test is array1, coincidence channel in array0 
	   // int chR_iext = channelMapping[iBar*2+1];// module under test is array1, coincidence channel in array0 
	   // if (opts.GetOpt<int>("Channels.array")==0) {
	   //   int chL_iext = channelMapping[iBar*2+0]+64;// module under test is array0, coincidence channel in array1
	   //   int chR_iext = channelMapping[iBar*2+1]+64;// module under test is array0, coincidence channel in array1
	   // }
	    float energyL_iext = (*energy)[channelIdx[chL_iext]];              
	    float energyR_iext = (*energy)[channelIdx[chR_iext]]; 
	    float totL_iext    = 0.001*(*tot)[channelIdx[chL_iext]];              
	    float totR_iext    = 0.001*(*tot)[channelIdx[chR_iext]]; 
	    if ( totL_iext > -10 && totL_iext < 100 && totR_iext > -10 && totR_iext < 100   ){
	      float energyMean=(energyL_iext+energyR_iext)/2;
	      if (energyMean>0){
		energySumArray+=energyMean;
		nActiveBarsArray+=1;
	      }
	    }
	  }
	  if (nActiveBarsArray > 3 ) continue; //Max number REF module active bars -> DEFAULT 3
	}
	
	energyL_ext = (*energy)[channelIdx[chL_ext]];
	energyR_ext = (*energy)[channelIdx[chR_ext]];
	
	int index( (10000*int(Vov*100.)) + (100*vth) + 99 );
	
	//--- create histograms, if needed
	if( h1_energyLR_ext[index] == NULL ){
	  c[index] = new TCanvas(Form("c1_Vov%.2f_th%02.0f",Vov,vth), Form("c1_Vov%.2f_th%02.0f",Vov,vth));
	  c[index] -> cd();
	  c[index] ->SetLogy();
	  h1_energyLR_ext[index] = new TH1F(Form("h1_energy_external_barL-R_Vov%.2f_th%02.0f",Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	}
	
	acceptEvent[entry] = true;
	
	h1_energyLR_ext[index] -> Fill(0.5*(energyL_ext + energyR_ext));    
	
      }


      
      std::cout << std::endl;
      
      if (!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") || !opts.GetOpt<std::string>("Input.sourceName").compare("Na22")){		
	for( auto index : h1_energyLR_ext){
	  rangesLR[index.first] = new std::vector<float>;
	  rangesLR[index.first]->push_back(30);
	  rangesLR[index.first]->push_back(950);
	  
	  if(!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar")) Na22SpectrumAnalyzerSingleBar(index.second,rangesLR[index.first]);
	  if(!opts.GetOpt<std::string>("Input.sourceName").compare("Na22")) Na22SpectrumAnalyzer(index.second,rangesLR[index.first]);
	}
      }
      
      if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB")){
	for( auto index : h1_energyLR_ext){
	  rangesLR[index.first] = new std::vector<float>;


	  float Vov = float ((int(index.first /10000))/100.);
	  float vth1 = float(int((index.first-Vov*10000*100)/100.));
	  float vth2 = float(int((step2-10000*(vth1+1))/100.)-1);
          float vthE = float(int((step2-10000*(vth1+1)) - 100*(vth2+1))-1);
          float vth = 0;
       
	  //std::cout << step2 << " ith1:1 " << vth1 << " ith2: " << vth2 << " E: " << vthE << std::endl;
       
	  if(!opts.GetOpt<std::string>("Input.vth").compare("vth1"))  { vth = vth1;}
          if(!opts.GetOpt<std::string>("Input.vth").compare("vth2"))  { vth = vth2;}
	  if(!opts.GetOpt<std::string>("Input.vth").compare("vthE"))  { vth = vthE;}
	  
	  index.second->GetXaxis()->SetRangeUser(200,900);
	  //if( opts.GetOpt<int>("Channels.array") == 0){
	  //  //	    index.second->GetXaxis()->SetRangeUser(50,900);
	  //  index.second->GetXaxis()->SetRangeUser(200,900);
	  //}
	  //if( opts.GetOpt<int>("Channels.array") == 1){
	  //  index.second->GetXaxis()->SetRangeUser(200,900);
	  //}

	  float max = index.second->GetBinCenter(index.second->GetMaximumBin());
	  index.second->GetXaxis()->SetRangeUser(0,1024);
	  
	  TF1* f_pre = new TF1(Form("fit_energy_coincBar_Vov%.2f_vth1_%02.0f",Vov,vth), "[0]*TMath::Landau(x,[1],[2])", 0, 1000.); 
	  f_pre -> SetRange(max*0.70, max*1.5);
	  f_pre -> SetLineColor(kBlack);
          f_pre -> SetLineWidth(2);
          f_pre -> SetParameters(index.second->Integral(index.second->GetMaximumBin(), index.second->GetNbinsX())/10, max, 0.1*max);
	  f_pre -> SetParLimits(1, 0, 9999);
	  f_pre -> SetParLimits(2, 0, 9999);
	  index.second->Fit(f_pre, "QRS+");      
	  
	  if (f_pre->GetParameter(1)>20)
	    rangesLR[index.first] -> push_back( 0.70*f_pre->GetParameter(1));
	  else
	    rangesLR[index.first] -> push_back( 20 );
	  

	  rangesLR[index.first] -> push_back( 950 );
	  
	  std::cout << "Vov = " << Vov << "  vth1 = " << vth1 << "   vth2 = " << vth2 
		    << "    Coincidence bar - energy range:  " << rangesLR[index.first]->at(0) << " - " << rangesLR[index.first]->at(1)<< std::endl;
	}
      }
    }
  
  
  






  //------------------------
  //--- 1st loop over events
  
  ModuleEventClass anEvent;

  // DUT:  
  unsigned short qfineL[16];
  unsigned short qfineR[16];    
  float totL[16];
  float totR[16];
  long long timeL[16];
  long long timeR[16];
  unsigned short t1fineL[16]; 
  unsigned short t1fineR[16]; 
  unsigned short t1coarseL[16];
  unsigned short t1coarseR[16];
  float qT1L[16]; 
  float qT1R[16]; 
  float energyL[16];
  float energyR[16];

  // REF(bar 7):
  long long timeL_ext;
  long long timeR_ext;
  float energyL_ext;
  float energyR_ext;
  float totL_ext;
  float totR_ext;
  unsigned short t1fineL_ext;
  unsigned short t1fineR_ext;
  unsigned short t1coarseL_ext;
  unsigned short t1coarseR_ext;

  // MCP data:
  double mcp_peak_time_evt;
  double mcp_peak_amp_evt;
  double mcp_phi_peak_evt;
  double mcp_phi_peak_from_edge_evt;
  double mcp_trigger_time_evt; 
  double mcp_trigger_offset_ps_evt;
  double mcp_phi_trigger_evt;
  double mcp_phi_trigger_from_edge_evt;
  double mcp_t0_abs_ps_evt;






  // --- Event loop with both data and MCP tree:

   int nEntries = tree->GetEntries();
   if( maxEntries > 0 ) nEntries = maxEntries;

   int mcpEntry = 0;
   int nEntries_MCP = tree_MCP->GetEntries();

   std::cout << "|| nEntries: " << nEntries << "||  nEntries_MCP: " << nEntries_MCP << std::endl;

   int pass_count = 0; 

   for(int entry = 0; entry < nEntries; ++entry) {

      bool hasMCP = false;   
      Long64_t globalEntry = entry;
      tree->GetEntry(globalEntry);
      int treeNumber = tree->GetTreeNumber();
      Long64_t localEntry = globalEntry - tree->GetTreeOffset()[treeNumber];
      TFile* currentFile = tree->GetCurrentFile();
      std::string fileName = currentFile ? currentFile->GetName() : "UNKNOWN";

      tree_MCP->GetEntry(mcpEntry);

      if(mcp_index == localEntry) 
      {
	      hasMCP = true;
	      mcpEntry++;
		      
	      if (Debug_on) {	      
		      std::cout
			      << "[DEBUG] file: " << fileName << "\n"
			      << "        treeNumber: " << treeNumber
			      << " | globalEntry: " << globalEntry
			      << " | localEntry: " << localEntry
			      << " | mcpEntry: " << mcpEntry
			      << " | mcp_index: " << mcp_index
			      << std::endl;
	      }
      }

      //-------------------------------------------------------------------------------
      
      if(hasMCP==false) continue; //coincidence with MCP required
     
     //-------------------------------------------------------------------------------


      if (useTrackInfo && nhits > 0 &&  (x < -100 || y < -100 ) ) continue;

      //std::cout << " ok, MCP selection" << std::endl;
      pass_count ++;

      float Vov = step1;
      float vth1 = float(int(step2/10000)-1);
      float vth2 = int((step2-10000*(vth1+1))/100.)-1;
      float vthE = float(int((step2-10000*(vth1+1)) - 100*(vth2+1))-1);
      float vth = 0;
      std::string vthMode = opts.GetOpt<std::string>("Input.vth");
      if(!opts.GetOpt<std::string>("Input.vth").compare("vth1"))  { vth = vth1;}
      if(!opts.GetOpt<std::string>("Input.vth").compare("vth2"))  { vth = vth2;}
      if(!opts.GetOpt<std::string>("Input.vth").compare("vthE"))  { vth = vthE;}
      // float vthe = float(int((step2-10000*vth1-step2-100*vth2)/1)-1);
      if( entry%20000 == 0 ){
      std::cout << step2 << " ith1: " << vth1 << " ith2: " << vth2 << " E: " << vthE << std::endl;
      }
    
     // select only one OV
     if (my_step1 > 0  && my_step1 != step1) continue;

    
     // --- check coincidence with another channel 
     if(!opts.GetOpt<std::string>("Coincidence.status").compare("yes"))
       {
	if(!acceptEvent[entry] ) continue;
	
	int chL_ext = opts.GetOpt<int>("Coincidence.chL");
	int chR_ext = opts.GetOpt<int>("Coincidence.chR");
	energyL_ext = (*energy)[channelIdx[chL_ext]];
	energyR_ext = (*energy)[channelIdx[chR_ext]];
	
	int label = (10000*int(Vov*100.)) + (100*vth) + 99;
	int eBin = opts.GetOpt<int>("Coincidence.peak511eBin");
	float avEn = 0.5 * ( energyL_ext + energyR_ext);
	if ( (!opts.GetOpt<std::string>("Input.sourceName").compare("Na22SingleBar") || !opts.GetOpt<std::string>("Input.sourceName").compare("Na22")) &&  avEn > rangesLR[label]-> at(eBin)) {
	  continue;
	}
	
	if ( (!opts.GetOpt<std::string>("Input.sourceName").compare("TB")) && ( avEn < rangesLR[label]-> at(0) || avEn > rangesLR[label]-> at(1) ) ) {
	  continue;
	}
      
  

        // REF(bar 7):
        timeL_ext = (*time)[channelIdx[chL_ext]];
        timeR_ext = (*time)[channelIdx[chR_ext]];
        t1fineL_ext = (*t1fine)[channelIdx[chL_ext]];
        t1fineR_ext = (*t1fine)[channelIdx[chR_ext]];
	t1coarseL_ext = (*t1coarse)[channelIdx[chL_ext]];
	t1coarseR_ext = (*t1coarse)[channelIdx[chR_ext]];
        totL_ext = 0.001*(*tot)[channelIdx[chL_ext]]; 
        totR_ext = 0.001*(*tot)[channelIdx[chR_ext]];
   

        // === Find max-energy bar in REF (after coincidence cuts) ===
        int maxBar_REF = -1;
        float maxE_REF = -999.;

        int index_REF( (10000*int(Vov*100.)) + (100*vth) );
        for(int iBar = 0; iBar < 16; iBar++){
           int chL_ref = channelMapping[iBar*2+0] + REFasic*32;
           int chR_ref = channelMapping[iBar*2+1] + REFasic*32;

           float Emean_ref = 0.5 * ( (*energy)[channelIdx[chL_ref]] + (*energy)[channelIdx[chR_ref]] );
    
	   if(Emean_ref > maxE_REF){
              maxE_REF = Emean_ref;
              maxBar_REF = iBar;
	   }
	}
    
	if (maxBar_REF < -1) continue;
        //if(maxBar_REF != 7) continue;
        
        if( h1_maxBar_REF[index_REF] == NULL ) {
	    	h1_maxBar_REF[index_REF] = new TH1F(Form("h1_maxBar_REF_Vov%.2f_th%02.0f",Vov,vth),"",16,-0.5,15.5);
	}

     h1_maxBar_REF[index_REF]->Fill(maxBar_REF);

   }//if coincidence == yes

 




    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar){
      
      if (channelIdx[chL[iBar]] >=0 && channelIdx[chR[iBar]] >=0){
	  
        qfineL[iBar]=(*qfine)[channelIdx[chL[iBar]]];
	qfineR[iBar]=(*qfine)[channelIdx[chR[iBar]]];
	totL[iBar]=0.001*(*tot)[channelIdx[chL[iBar]]];
	totR[iBar]=0.001*(*tot)[channelIdx[chR[iBar]]];
	energyL[iBar]=(*energy)[channelIdx[chL[iBar]]];
	energyR[iBar]=(*energy)[channelIdx[chR[iBar]]];
	timeL[iBar]=(*time)[channelIdx[chL[iBar]]];
	timeR[iBar]=(*time)[channelIdx[chR[iBar]]];
	t1fineL[iBar]=(*t1fine)[channelIdx[chL[iBar]]];
	t1fineR[iBar]=(*t1fine)[channelIdx[chR[iBar]]];
        t1coarseL[iBar]=(*t1coarse)[channelIdx[chL[iBar]]];
        t1coarseR[iBar]=(*t1coarse)[channelIdx[chR[iBar]]];

	}
      else
	{
	  qfineL[iBar]=-10;
	  qfineR[iBar]=-10;
	  totL[iBar]=-9999;
	  totR[iBar]=-9999;
	  energyL[iBar]=-10;
	  energyR[iBar]=-10;
	  timeL[iBar]=-10;
	  timeR[iBar]=-10;
	  t1fineL[iBar]=-10;
	  t1fineR[iBar]=-10;
	} 
    }// end loop over bars
    
    
    
    
    
    
    
    int maxEn=0;
    int maxBar=0;

    float energySumArray = 0;
    int   nActiveBarsArray = 0;
    int nBarsVeto[16];
	  
    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar)
      {
	nBarsVeto[iBar] = 0;
	
	if (totL[iBar]>-10 && totR[iBar]>-10 && totL[iBar]<100 && totR[iBar]<100)
	  {
	  
	    float energyMean=(energyL[iBar]+energyR[iBar])/2;
	    if (energyL[iBar]>0 && energyR[iBar]>0 && energyMean > 0){
	      //if (energyMean > 0){
	      energySumArray+=energyMean;
	      nActiveBarsArray+=1;
	    }
	    

	    // check energy in adjacent bars
	    for (int jBar = int(iBar) - 2; jBar < int(iBar) + 3; ++jBar){
	      if (jBar == int(iBar)) continue;
	      if (jBar < 0 || jBar > 15 ) continue;
	      if (totL[jBar]<-10 || totL[jBar]>100) continue;
	      if (totR[jBar]<-10 || totR[jBar]>100) continue;
	      float en = (energyL[jBar]+energyR[jBar])/2;
	      //if ( en > minE[std::make_pair(jBar, Vov)] && en<1024 ){
	      if ( en > minE[std::make_pair(jBar, Vov)] && minE[std::make_pair(jBar, Vov)]>1 && en<1024 ){
		nBarsVeto[iBar]+=1;
	      }
	    }
	    
	    // find max bar
	    if(energyMean>maxEn){
	      maxEn = energyMean;
	      maxBar = iBar;

	    }
	  }
      }// end loop over bars

    
    for(unsigned int iBar = 0; iBar < channelMapping.size()/2; ++iBar)                                                                                                                                  
      {
        if (totL[iBar]>-10 && totR[iBar]>-10 && totL[iBar]<100 && totR[iBar]<100){  
    
	int index( (10000*int(Vov*100.)) + (100*vth) + iBar );
	
	
	//--- create histograms, if needed
	if( h1_totL[index] == NULL )
	{
	  h1_qfineL[index] = new TH1F(Form("h1_qfine_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",512,-0.5,511.5);
	  h1_qfineR[index] = new TH1F(Form("h1_qfine_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",512,-0.5,511.5);
	  
	  h1_totL[index] = new TH1F(Form("h1_tot_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",400,-5.,35.);
	  h1_totR[index] = new TH1F(Form("h1_tot_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",400,-5.,35.);
	  
	  h1_energyL[index] = new TH1F(Form("h1_energy_bar%02dL_Vov%.2f_th%02.0f",iBar,Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	  h1_energyR[index] = new TH1F(Form("h1_energy_bar%02dR_Vov%.2f_th%02.0f",iBar,Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	  
	  outTrees[index] = new TTree(Form("data_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),Form("data_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth));
	  outTrees[index] -> Branch("event",&anEvent);
	  
	  h1_energyLR[index] = new TH1F(Form("h1_energy_bar%02dL-R_Vov%.2f_th%02.0f",iBar,Vov,vth),"",map_energyBins[Vov],map_energyMins[Vov],map_energyMaxs[Vov]);
	}	
	 

 
	if( h1_MeanTimeREF[iBar] == NULL )
	{	
		h1_MeanTimeREF[iBar] = new TH1F(Form("h1_MeanTimeREF_bar%02dL-R",iBar),"",200,-10.,10);
	}
      
	




	
	//--- fill histograms for each bar for Co60 & TB analysis
	if( !opts.GetOpt<std::string>("Input.sourceName").compare("TB"))
	  {
	  if( totL[iBar] <= -10. || totR[iBar] <= -10. ) continue;
	  if( totL[iBar] >= 50. ||  totR[iBar] >= 50.) continue;
	  if( ( thrZero.GetThresholdZero(chL[iBar],vthMode) + vth) > 63. ) continue;
          if( ( thrZero.GetThresholdZero(chR[iBar],vthMode) + vth) > 63. ) continue;

	  //if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB") && (energySumArray > 900 || nActiveBarsArray > 5)) continue; // to remove showering events
	  //if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB") && (vetoOtherBars && nBarsVeto[iBar] > 0)) continue; // to remove showering events/ cross talk
	  //if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB") && (vetoOtherBars && nActiveBarsArray > 5)) continue; // to remove showering events
	  //if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB") && (vetoOtherBars && nActiveBarsArray > 3)) continue; // to remove showering events
	  int maxActiveBars = 3; //DEFAULT 3
	  if (!opts.GetOpt<std::string>("Input.sourceName").compare("TB") && (vetoOtherBars && nActiveBarsArray > maxActiveBars)) continue; // to remove showering events
	  

	  h1_qfineL[index] -> Fill( qfineL[iBar] );
	  h1_totL[index] -> Fill( totL[iBar]  );
	  h1_energyL[index] -> Fill( energyL[iBar] );
	  h1_qfineR[index] -> Fill( qfineR[iBar] );
	  h1_totR[index] -> Fill( totR[iBar] );
	  h1_energyR[index] -> Fill( energyR[iBar] );
	  h1_energyLR[index] -> Fill(0.5*(energyL[iBar]+energyR[iBar]));
	  h1_MeanTimeREF[iBar] -> Fill((0.5*(timeL_ext+timeR_ext))/(10e14));

	  
	  //  REF(bar 7):
	  anEvent.timeL_ext = timeL_ext;
          anEvent.timeR_ext = timeR_ext;
          anEvent.t1fineL_ext = t1fineL_ext;
          anEvent.t1fineR_ext = t1fineR_ext;
          anEvent.totL_ext = totL_ext;
          anEvent.totR_ext = totR_ext;
	  anEvent.energyL_ext = energyL_ext;
	  anEvent.energyR_ext = energyR_ext;
	  anEvent.t1coarseL_ext = t1coarseL_ext;
	  anEvent.t1coarseR_ext = t1coarseR_ext;


          // DUT data:
	  anEvent.barID = iBar;
	  anEvent.Vov = Vov;
	  anEvent.vth1 = vth;
	  anEvent.energyL = energyL[iBar];
	  anEvent.energyR = energyR[iBar];
	  anEvent.totL = totL[iBar];
	  anEvent.totR = totR[iBar];
	  anEvent.timeL = timeL[iBar];
	  anEvent.timeR = timeR[iBar];
	  anEvent.t1fineL = t1fineL[iBar];
	  anEvent.t1fineR = t1fineR[iBar];
	  if(useTrackInfo){
	    anEvent.nhits = nhits;
	    anEvent.x = x;
	    anEvent.y = y;
	  }
	  else{
	    anEvent.nhits = -1;
	    anEvent.x = -999.;
	    anEvent.y = -999.;
	  }



	  //MCP data:
	  mcp_peak_time_evt= mcp_peak_time;
          mcp_peak_amp_evt= mcp_peak_amp;
          mcp_phi_peak_evt= mcp_phi_peak;
	  mcp_phi_peak_from_edge_evt = mcp_phi_peak_from_edge;
	  mcp_trigger_time_evt = mcp_trigger_time;
	  mcp_trigger_offset_ps_evt = mcp_trigger_offset_ps;
	  mcp_phi_trigger_evt = mcp_phi_trigger;
	  mcp_phi_trigger_from_edge_evt = mcp_phi_trigger_from_edge;
          mcp_t0_abs_ps_evt = mcp_t0_abs_ps;


	  anEvent.mcp_peak_time = mcp_peak_time_evt;
          anEvent.mcp_peak_amp = mcp_peak_amp_evt;
          anEvent.mcp_phi_peak = mcp_phi_peak_evt;
          anEvent.mcp_phi_peak_from_edge = mcp_phi_peak_from_edge_evt;
          anEvent.mcp_trigger_time = mcp_trigger_time_evt;
	  anEvent.mcp_trigger_offset_ps = mcp_trigger_offset_ps_evt;
	  anEvent.mcp_phi_trigger = mcp_phi_trigger_evt;
	  anEvent.mcp_phi_trigger_from_edge = mcp_phi_trigger_from_edge_evt;
	  anEvent.mcp_t0_abs_ps = mcp_t0_abs_ps_evt;


	  outTrees[index] -> Fill();
	  }	  
	}
	}// -- end loop over bars



    // --- for TB max bar plots
    if( !opts.GetOpt<std::string>("Input.sourceName").compare("TB")) {

	    int index_1( (10000*int(Vov*100.)) + (100*vth) );
	    if (maxBar == 0) continue;
	    if( h1_maxBar_DUT[index_1] == NULL )
	    {
		    h1_maxBar_DUT[index_1] = new TH1F(Form("h1_maxBar_DUT_Vov%.2f_th%02.0f",Vov,vth),"",16,-0.5,15.5);
	    }

	    h1_maxBar_DUT[index_1] -> Fill( maxBar);
    }

	   

	  } // --- end loop over events
  
// std::cout << "pass_count:" << pass_count << std::endl;

  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
