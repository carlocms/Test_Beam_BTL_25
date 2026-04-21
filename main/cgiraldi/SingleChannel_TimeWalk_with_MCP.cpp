#include "interface/AnalysisUtils.h"
#include "interface/Na22SpectrumAnalyzer.h"
//#include "interface/Na22SpectrumAnalyzerSingleBar.h"
#include "interface/Na22SpectrumAnalyzerSingleBar_TOFHIR2.h"
//#include "interface/Na22SpectrumAnalyzerModule_TOFHIR2.h"
#include "interface/Co60SpectrumAnalyzer_2Peaks.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <time.h>
#include <stdio.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <algorithm>
#include <iterator>

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
#include "TLegend.h"
#include "TSpectrum.h"
#include <TPaveStats.h>


Double_t langaufun(Double_t *x, Double_t *par)
{
  
  //Fit parameters:
  //par[0]=Width (scale) parameter of Landau density
  //par[1]=Most Probable (MP, location) parameter of Landau density
  //par[2]=Total area (integral -inf to inf, normalization constant)
  //par[3]=Width (sigma) of convoluted Gaussian function
  //
  //In the Landau distribution (represented by the CERNLIB approximation),
  //the maximum is located at x=-0.22278298 with the location parameter=0.
  //This shift is corrected within this function, so that the actual
  //maximum is identical to the MP parameter.
  
  // Numeric constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location
  
  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas
  
  // Variables
  Double_t xx;
  Double_t mpc;
  Double_t fland;
  Double_t sum = 0.0;
  Double_t xlow,xupp;
  Double_t step;
  Double_t i;
  
  
  // MP shift correction
  mpc = par[1] - mpshift * par[0];
  
  // Range of convolution integral
  xlow = x[0] - sc * par[3];
  xupp = x[0] + sc * par[3];
  
  step = (xupp-xlow) / np;
  
  // Convolution integral of Landau and Gaussian by sum
  for(i=1.0; i<=np/2; i++) {
    xx = xlow + (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
    
    xx = xupp - (i-.5) * step;
    fland = TMath::Landau(xx,mpc,par[0]) / par[0];
    sum += fland * TMath::Gaus(x[0],xx,par[3]);
  }
  
  return (par[2] * step * sum * invsq2pi / par[3]);
}





//---- find energy bins
void GetEnergyBins(TH1F *h, std::vector<float> *r, std::map<int, float> & b){

  for(unsigned int i = 1; i < r->size(); i++){
    TH1F *binHisto = new TH1F ( "binHisto", "binHisto", h -> FindBin(r->at(i)) - h->FindBin(r-> at(i-1)), r-> at(i-1), r->at(i));
    int j = 1;
    for (int bin = h->FindBin(r->at(i-1)) ; bin < h -> FindBin(r->at(i))+1 ; bin++){
      binHisto -> SetBinContent( j, h->GetBinContent(bin));
      j++;
    }
    b[i] = binHisto -> GetMean();
    binHisto -> Delete();
  }
}


//---- Find fit range profile
std::pair<float,float> GetDynamicFitRange(TProfile* prof, float marginFraction = 0.05)
{
    if (!prof) return {0., 1.};

    int firstBin = -1;
    int lastBin  = -1;

    // find first and last non-empty bin
    for (int i = 1; i <= prof->GetNbinsX(); ++i)
    {
        if (prof->GetBinEntries(i) > 0)
        {
            if (firstBin < 0) firstBin = i;
            lastBin = i;
        }
    }

    // default range if no bins populated
    if (firstBin < 0 || lastBin < firstBin)
        return {0., 1.};

    // compute x-range
    float xmin = prof->GetBinLowEdge(firstBin);
    float xmax = prof->GetBinLowEdge(lastBin + 1);

    // apply margin
    float margin = marginFraction * (xmax - xmin);
    xmin += margin;
    xmax -= margin;

    return {xmin, xmax};
}



float GetProfileBinCorrection(TProfile* prof, float x)
{
    if (!prof) return 0.0;

    int bin = prof->FindBin(x);

    if (bin < 1 || bin > prof->GetNbinsX())
        return 0.0;

    if (prof->GetBinEntries(bin) <= 0)
        return 0.0;

    return prof->GetBinContent(bin);
}








// ============  ********* MAIN  ************ =============== //
int main(int argc, char** argv)
{
  setTDRStyle();
  float cpu[2]{0}, mem[2]={0}, vsz[2]={0}, rss[2]={0};

  gErrorIgnoreLevel = kError;

  typedef std::numeric_limits<double> dbl;
        std::cout.precision(dbl::max_digits10);
  if( argc < 2 )
  {
    std::cout << ">>> SingleChannel_TimeWalk::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }

  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);

  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);


  //--- get parameters
  // --- From Loop 1:
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir_withMCP");
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop1_tot_histo/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop1_energy_TW/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop1_MeanTimeREF/",plotDir.c_str()));
 
  
  // --- From Loop 2:
  system(Form("mkdir -p %s/Loop2_t1fineL/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_t1fineR/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_mean_phase_bar/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_t1fineL_vs_t1fineR/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_MCP_phase/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_phase_MeanBar_vs_MCP/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_DeltaPhase_meanBar_MCP/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_phase_barL/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_phase_barR/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_phase_bar_L_vs_R/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_DeltaPhase_Bar_LR/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_REF_global/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_DeltaPhase_REF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_NumberEventsBar/",plotDir.c_str()));
 
  // --- From Loop3:
  system(Form("mkdir -p %s/Loop3_Profile_DeltaPhase_vs_MeanEnergyBar/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_Scatter_DeltaPhase_vs_MeanEnergyBar/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_Scatter_DeltaPhaseR_vs_EnergyBarR/",plotDir.c_str())); 
  system(Form("mkdir -p %s/Loop3_Scatter_DeltaPhaseL_vs_EnergyBarL/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_DeltaPhase_meanBar_DUT_REF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_DeltaPhase_MCP_MeanBarDUT/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_DeltaPhase_MCP_MeanBarREF/",plotDir.c_str()));

  //TOTAL:
  system(Form("mkdir -p %s/TOTAL_DeltaPhase_MCP_MeanBar_total/",plotDir.c_str()));
  system(Form("mkdir -p %s/TOTAL_Scatter_DeltaPhase_vs_MeanEnergyBar_total/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_DeltaPhase_MCP_MeanBar_corr/",plotDir.c_str()));

  // --- From Loop 4
  system(Form("mkdir -p %s/Loop4_TOTAL_DeltaPhase_MCP_MeanBar_total_corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_Scatter_DeltaPhase_vs_MeanEnergyBarDUT_corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_Profile_DeltaPhase_vs_MeanEnergyBar_corr/",plotDir.c_str()));



  std::vector<std::string> LRLabels;
  LRLabels.push_back("L");
  LRLabels.push_back("R");
  LRLabels.push_back("L-R");

  bool Debug_on = opts.GetOpt<bool>("Input.debug_mod"); 
  std::vector<float> Vov = opts.GetOpt<std::vector<float> >("Plots.Vov");
  std::vector<int> energyMins = opts.GetOpt<std::vector<int> >("Plots.energyMins");
  std::vector<int> energyMaxs = opts.GetOpt<std::vector<int> >("Plots.energyMaxs");

  std::map<float,int> map_energyMins;
  std::map<float,int> map_energyMaxs;
  for(unsigned int ii = 0; ii < Vov.size(); ++ii)
    {
      map_energyMins[Vov[ii]] = energyMins[ii];
      map_energyMaxs[Vov[ii]] = energyMaxs[ii];
    }

  int useTrackInfo = opts.GetOpt<int>("Input.useTrackInfo");


  // -- read minimum energy for each bar from file
  std::string minEnergiesFileName = opts.GetOpt<std::string>("Cuts.minEnergiesFileName");
  std::map < std::pair<int, float>, float> minE;
  std::cout << minEnergiesFileName << std::endl;
  if( minEnergiesFileName != "" )
    {
      std::ifstream minEnergiesFile;
      minEnergiesFile.open(minEnergiesFileName);
      std::string line;
      int bar;
      float ov;
      float value;
      while( minEnergiesFile.good() )
        {
          getline(minEnergiesFile, line);
          std::istringstream ss(line);
          ss >> bar >> ov >> value;
          minE[std::make_pair(bar,ov)] = value;
          //std::cout << "minEnergies:   bar " << bar <<  "   Vov " << ov << "   E " << minE[std::make_pair(bar,ov)] <<std::endl;
        }
    }
  else
    {
      for(unsigned int iBar = 0; iBar < 16; ++iBar)
        for(unsigned int ii = 0; ii < Vov.size(); ++ii)
          minE[std::make_pair(iBar, Vov[ii])] = map_energyMins[Vov[ii]];
    }


  //--- open files
  std::string step1FileName= opts.GetOpt<std::string>("Input.step1FileName_wMCP");
  TFile* inFile = TFile::Open(step1FileName.c_str(),"READ");

  std::map<std::string,TTree*> trees;

  std::map<std::string,int> VovLabels;
  std::map<std::string,int> thLabels;
  std::vector<std::string> stepLabels;
  std::map<std::string,float> map_Vovs;
  std::map<std::string,float> map_ths;

  TList* list = inFile -> GetListOfKeys();
  TIter next(list);
  TObject* object = 0;
  while( (object = next()) )
  {
    std::string name(object->GetName());
    std::vector<std::string> tokens = GetTokens(name,'_');
    std::size_t found;

    found = name.find("data_");
    //tree
    if( found!=std::string::npos )
    {
      std::string label(Form("%s_%s_%s",tokens[1].c_str(),tokens[2].c_str(),tokens[3].c_str()));
      trees[label] = (TTree*)( inFile->Get(name.c_str()) );
    }
    found = name.find("h1_energy_b");
    if( found!=std::string::npos )
    {
     //Vov e th
      std::string stepLabel = tokens[3]+"_"+tokens[4];
      VovLabels[tokens[3]] += 1;
      thLabels[tokens[4]] += 1;
      stepLabels.push_back(stepLabel);
      std::string string_Vov = tokens[3];
      string_Vov.erase(0,3);
      map_Vovs[stepLabel] = atof(string_Vov.c_str());
      std::string string_th = tokens[4];
      string_th.erase(0,2);
      map_ths[stepLabel] = atof(string_th.c_str());
    }
  }
  std::sort(stepLabels.begin(),stepLabels.end());
  stepLabels.erase(std::unique(stepLabels.begin(),stepLabels.end()),stepLabels.end());


  //--- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameSingleChannel_withMCP");
  TFile* outFile = TFile::Open(outFileName.c_str(),"RECREATE");
  outFile->cd();


  //std::map<double,TH1F*> h1_mean_time_ext;
  
  //Loop 2:
  std::map<double,TH1F*> h1_t1fineL;
  std::map<double,TH1F*> h1_t1fineR;
  std::map<double,TH2F*> h2_t1fineL_vs_t1fineR;
  std::map<double,TH1F*> h1_mean_phase_bar;
  std::map<double,TH1F*> h1_MPC_phase;
  std::map<double,TH2F*> h2_phase_MeanBar_vs_MCP;
  std::map<double,TH1F*> h1_DeltaPhase_meanBar_MCP;
  std::map<double,TH1F*> h1_DeltaPhase_Bar_LR;
  std::map<double,TH1F*> h1_phase_barL;
  std::map<double,TH1F*> h1_phase_barR;
  std::map<double,TH2F*> h2_phase_bar_L_vs_R;
  std::map<double,TH1F*> h1_DeltaPhase_REF;
  TH1F* h1_NumberEventsBar;
  TH1F* h1_DeltaPhaseLR_REF_global;

  //Loop 3:
  std::map<double,TProfile*> p1_DeltaPhase_vs_MeanEnergyBar;
  std::map<double,TH2F*> h2_DeltaPhase_vs_MeanEnergyBar;
  std::map<double,TProfile*> p1_DeltaPhaseL_vs_EnergyBarL;
  std::map<double,TProfile*> p1_DeltaPhaseR_vs_EnergyBarR;
  std::map<double,TH2F*> h2_DeltaPhaseR_vs_EnergyBarR;
  std::map<double,TH2F*> h2_DeltaPhaseL_vs_EnergyBarL;
  std::map<double,TH1F*> h1_DeltaPhase_meanBar_DUT_REF;
  std::map<double,TH1F*> h1_DeltaPhase_MCP_MeanBarDUT;
  std::map<double,TH1F*> h1_DeltaPhase_MCP_MeanBarREF;
  
  //total:
  std::map<double,TH1F*> h1_DeltaPhase_MCP_MeanBar_total;
  std::map<double,TH2F*> h2_DeltaPhase_vs_MeanEnergyBar_total;
  std::map<double,TProfile*> p1_DeltaPhase_vs_MeanEnergyBar_total;


  //Loop 4:
  std::map<double,TH1F*> h1_DeltaPhase_MCP_MeanBar_total_corr;
  std::map<double,TH2F*> h2_DeltaPhase_vs_MeanEnergyBar_total_corr;
  std::map<double,TProfile*> p1_DeltaPhase_vs_MeanEnergyBar_total_corr;
  std::map<double,TH1F*> h1_DeltaPhase_MCP_MeanBarDUT_corr;
  std::map<double,TProfile*> p1_DeltaPhase_vs_MeanEnergyBarDUT_corr;
  std::map<double,TH2F*> h2_DeltaPhase_vs_MeanEnergyBarDUT_corr;


  std::map<std::string, std::map<int, std::vector<float>*> > ranges; //ranges[LRlabel][index]
  std::map<std::string, std::map<int, std::map<std::string,std::pair<float,float> > > > peaks;  //peaks[LRlabel][index][energyPeak]
  std::map<std::string, std::map<int, std::map<int,float> > > energyBin; // energyBin[LRlabel][index]

  std::map<int,TF1*>  f_langaus; // f_langaus[index]
  std::map<int,TF1*>  f_gaus; // f_gaus[index]
  std::map<int,TF1*>  f_landau; // f_gaus[index]
  std::map<float,int>  Vov_LandauMin; //Vov_LandauMin[Vov]
  Vov_LandauMin[1.0] = 0;
  Vov_LandauMin[1.5] = 50;
  Vov_LandauMin[2.0] = 100;
  Vov_LandauMin[2.5] = 180;
  Vov_LandauMin[3.0] = 250;
  Vov_LandauMin[3.5] = 300;
  Vov_LandauMin[6.0] = 300;


  //--- get plot settings
  TCanvas* c;
  TCanvas* c2;
  float* vals = new float[6];
  TLatex* latex;
  TLatex* latex_tot;
  TH1F* histo;
  TProfile* prof;
  TH2F* h2;



  //------------------
  //--- draw 1st plots
  std::string source = opts.GetOpt<std::string>("Input.sourceName");
  std::string Na22 = "Na22";
  std::string Na22SingleBar = "Na22SingleBar";
  std::string Co60 = "Co60";
  std::string Co60SumPeak = "Co60SumPeak";
  std::string Laser = "Laser";
  std::string TB = "TB";
  std::string keepAll = "keepAll";
  std::vector<int> barList = opts.GetOpt<std::vector<int> >("Plots.barList");// list of bars to be analyzed read from cfg


  for(auto stepLabel : stepLabels)
    {
      //TrackProcess(cpu, mem, vsz, rss);

      float Vov = map_Vovs[stepLabel];
      float vth1 = map_ths[stepLabel];
      std::string VovLabel(Form("Vov%.2f",Vov));
      std::string thLabel(Form("th%02.0f",vth1));

      // -----------
      // --- external bar
      // -- draw energy
      c = new TCanvas(Form("c_energy_external_Vov%.2f_th%02.0f",Vov,vth1), Form("c_energy_external_Vov%.2f_th%02.0f",Vov,vth1));
      gPad -> SetLogy();
      histo = (TH1F*)( inFile->Get(Form("h1_energy_external_barL-R_Vov%.2f_th%02.0f", Vov, vth1)));
      if ( histo )
	{
	  histo -> SetTitle(";energy [a.u.];entries");
	  histo -> SetLineColor(kRed);
	  histo -> SetLineWidth(2);
	  histo -> Draw();

	  TF1 *ftemp = histo->GetFunction( ((histo->GetListOfFunctions()->FirstLink())->GetObject())->GetName() );
	  if (ftemp != NULL){
	    ftemp->SetNpx(1000);
	    ftemp->SetLineWidth(2);
	    ftemp->SetLineColor(1);
	    ftemp->Draw("same");
	  }
	  c -> Print(Form("%s/Loop1_energy_TW/c_energy_external__Vov%.2f_th%02.0f.png",plotDir.c_str(), Vov, vth1));
	  c -> Print(Form("%s/Loop1_energy_TW/c_energy_external__Vov%.2f_th%02.0f.pdf",plotDir.c_str(), Vov, vth1));
	  delete c;
	  delete ftemp;
	}

      //--------------------------------------------------------
      // --- loop over bars
      for(int iBar = 0; iBar < 16; ++iBar) {

	bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
	if (!barFound) continue;

	int index( (10000*int(Vov*100.)) + (100*vth1) + iBar );

	// -- loop over L, R, LR
	for(auto LRLabel : LRLabels ) {

	  //label histo
	  std::string label(Form("bar%02d%s_%s",iBar,LRLabel.c_str(),stepLabel.c_str()));

	  latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d%s}{V_{OV} = %.2f V, th. = %d DAC}",iBar,LRLabel.c_str(),Vov,int(vth1)));
	  if (LRLabel == "L-R") {
	    latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	  }
	  latex -> SetNDC();
	  latex -> SetTextFont(42);
	  latex -> SetTextSize(0.04);
	  latex -> SetTextColor(kRed);

	  // -- draw ToT histograms only for L, R
	  if (LRLabel == "R" || LRLabel == "L")
	    {
	      // -- draw ToT histograms
	      c = new TCanvas(Form("c_tot_%s",label.c_str()),Form("c_tot_%s",label.c_str()));
	      gPad -> SetLogy();

	      histo = (TH1F*)( inFile->Get(Form("h1_tot_%s",label.c_str())) );
	      if (histo)
		{
		  histo -> SetTitle(";TbT [ns];entries");
		  histo -> SetLineColor(kRed);
		  float max = histo->GetBinCenter(histo->GetMaximumBin());
		  histo -> Draw();
		  // TLine* line_totAcc1 = new TLine(cut_totAcc[chID][Vov],histo->GetMinimum(),cut_totAcc[chID][Vov],histo->GetMaximum());
		  // line_totAcc1 -> SetLineColor(kBlack);
		  // line_totAcc1 -> Draw("same");
		  latex -> Draw("same");

		  latex_tot = new TLatex(0.40,0.75,Form("max pos = %.2f ns", max));
                  latex_tot -> SetNDC();
                  latex_tot -> SetTextFont(42);
                  latex_tot -> SetTextSize(0.04);
                  latex_tot -> SetTextColor(kRed);
		  latex_tot -> Draw("same"); 
		  
		  histo -> Write();
		  c -> Print(Form("%s/Loop1_tot_histo/c_tot__%s.png",plotDir.c_str(),label.c_str()));
		  c -> Print(Form("%s/Loop1_tot_histo/c_tot__%s.pdf",plotDir.c_str(),label.c_str()));
		  delete c;
		}
	    }//end if LRLabel L || R

	  
	  // -- draw energy histograms
	  c = new TCanvas(Form("c_energy_%s",label.c_str()),Form("c_energy_%s",label.c_str()));
	  gPad -> SetLogy();

	  histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label.c_str())) );
	  if( !histo ) continue;
	  histo -> SetTitle(";energy [a.u.];entries");
	  histo -> SetLineColor(kRed);
	  histo -> SetLineWidth(2);
	  histo -> Draw();


	  // --- look for peaks and define energy ranges
	  if( source.compare(Na22) && source.compare(Na22SingleBar) && source.compare(Co60) && source.compare(Co60SumPeak) && source.compare(Laser) && source.compare(TB) && source.compare(keepAll) )
	    {
	      std::cout << " Source not found !!! " << std::endl;
	      return(0);
	    }

	  ranges[LRLabel][index] = new std::vector<float>;

	  // -- if MIP peak, we don't use the spectrum analyzers - just langaus fit to the energy peak.
	  if(!source.compare(TB)){
	    /*
	    //TF1* f_langaus = new TF1("f_langaus", langaufun, 300.,1000.,4);
	    f_langaus[index] = new TF1(Form("fit_energy_bar%02d_Vov%.2f_vth1_%02.0f",iBar,Vov,vth1),langaufun,Vov_LandauMin[Vov],1000.,4);
	    f_langaus[index] -> SetParameters(f_landau->GetParameter(2),f_landau->GetParameter(1),histo->Integral(histo->FindBin(300.),histo->FindBin(1000.))*histo->GetBinWidth(1),10.);
	    f_langaus[index] -> SetLineColor(kBlack);
	    f_langaus[index] -> SetLineWidth(2);
	    histo -> Fit(f_langaus[index],"QNRS+");
	    f_langaus[index] -> Draw("same");

	    ranges[LRLabel][index] -> push_back( 0.80*f_langaus[index]->GetParameter(1));
	    ranges[LRLabel][index] -> push_back( histo -> GetBinCenter(500) );
	    */

	    // if( opts.GetOpt<int>("Channels.array") == 1){
	    //   histo->GetXaxis()->SetRangeUser(minE[std::make_pair(iBar, Vov)], 950);
	    // }
	    // if( opts.GetOpt<int>("Channels.array") == 0){
	    //   histo->GetXaxis()->SetRangeUser(minE[std::make_pair(iBar, Vov)], 950);
	    // }
	    float max = histo->GetBinCenter(histo->GetMaximumBin());
	    //histo->GetXaxis()->SetRangeUser(0,1024);
	    histo->GetXaxis()->SetRangeUser(minE[std::make_pair(iBar, Vov)], 950); // minE to avoid fitting noise

	    f_gaus[index] = new TF1(Form("fit_energy_bar%02d%s_Vov%.2f_vth1_%02.0f",iBar,LRLabel.c_str(),Vov,vth1), "gaus", max-50, max+50);
	    f_gaus[index]->SetParameters(histo->GetMaximumBin(), max, 70);
	    histo->Fit(f_gaus[index], "QRS");
	    f_gaus[index]->SetRange(f_gaus[index]->GetParameter(1)-f_gaus[index]->GetParameter(2), f_gaus[index]->GetParameter(1)+f_gaus[index]->GetParameter(2));
	    histo->Fit(f_gaus[index], "QRS");
	    f_gaus[index] -> SetLineColor(kBlue);
	    f_gaus[index] -> SetLineWidth(2);
	    f_gaus[index] -> SetLineStyle(2);
	    //f_gaus[index] -> Draw("same");

	    //ranges[LRLabel][index] -> push_back( 0.80*f_gaus[index]->GetParameter(1));
	    //ranges[LRLabel][index] -> push_back( histo -> GetBinCenter(700) ); // to avoid sturation

	    f_landau[index] = new TF1(Form("f_landau_bar%02d%s_Vov%.2f_vth1_%02.0f", iBar,LRLabel.c_str(),Vov,vth1),"[0]*TMath::Landau(x,[1],[2])", 0,1000.);
	    float xmin = max * 0.65;
	    float xmax = std::min(max*2.5, 940.);
	    //float xmin = 0;
	    //float xmax = 1000;
	    f_landau[index] -> SetRange(xmin,xmax);
	    f_landau[index] -> SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, max, 0.1*max);
	    f_landau[index] -> SetParLimits(1,0,9999);
	    f_landau[index] -> SetParLimits(2,0,9999);
	    histo -> Fit(f_landau[index],"QRS");
	    if ( f_landau[index]->GetParameter(1) > 0 ){
	      xmin = f_landau[index]->GetParameter(1) - 2 * std::abs(f_landau[index]->GetParameter(2));
	      if (xmin < minE[std::make_pair(iBar, Vov)]) xmin = minE[std::make_pair(iBar, Vov)] ;
	      xmax = std::min(f_landau[index]->GetParameter(1) * 2.5, 940.);
	      f_landau[index] -> SetRange(xmin, xmax);
	      f_landau[index] -> SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, f_landau[index]->GetParameter(1), 0.1*f_landau[index]->GetParameter(1));
	    }
	    histo -> Fit(f_landau[index],"QRS");

	    f_landau[index] -> SetLineColor(kBlack);
	    f_landau[index] -> SetLineWidth(2);
	    f_landau[index] -> Draw("same");

	    if ( f_landau[index]->GetNDF() >0 && f_landau[index]->GetParameter(1) > minE[std::make_pair(iBar, Vov)] &&
		 (f_landau[index]->GetParameter(1) - 2.0 * std::abs(f_landau[index]->GetParameter(2))) >=  minE[std::make_pair(iBar, Vov)] &&
		 (f_landau[index]->GetParameter(1) - 2.0 * std::abs(f_landau[index]->GetParameter(2))) < 950) {
	      //ranges[LRLabel][index] -> push_back( 0.75*f_landau[index]->GetParameter(1));
	      //ranges[LRLabel][index] -> push_back( 0.60*f_landau[index]->GetParameter(1));
	      ranges[LRLabel][index] -> push_back( f_landau[index]->GetParameter(1) - 2.0 * std::abs(f_landau[index]->GetParameter(2)));
	      //ranges[LRLabel][index] -> push_back( 0 );
	    }
	    else
	      ranges[LRLabel][index] -> push_back( minE[std::make_pair(iBar, Vov)] ); //

	    if ( LRLabel=="L-R" && int(vth1)==10) std::cout << iBar << "  " << Vov  << "  " << ranges[LRLabel][index] ->at(0) <<std::endl;


	    ranges[LRLabel][index] -> push_back( std::min(f_landau[index]->GetParameter(1)*2.0, 940.)); // tight selection around the MIP peak
	    //ranges[LRLabel][index] -> push_back( 940 ); // use the entire mip spectrum -> DEFAULT CHOICE




	    for(auto range: (*ranges[LRLabel][index])){
	      TLine* line = new TLine(range,0.,range, histo->GetMaximum());
	      line -> SetLineWidth(2);
	      line -> SetLineStyle(7);
	      line -> Draw("same");
	    }

	    GetEnergyBins(histo, ranges[LRLabel][index], energyBin[LRLabel][index]);
	  }// end MIP (TB)



	  // -- draw and print energy plots

	  histo->GetXaxis()->SetRangeUser(0,1024);
	  latex -> Draw("same");
	  outFile -> cd();
	  histo->Write();
	  c -> Update();
	  c -> Print(Form("%s/Loop1_energy_TW/c_energy__%s.png",plotDir.c_str(),label.c_str()));
	  c -> Print(Form("%s/Loop1_energy_TW/c_energy__%s.pdf",plotDir.c_str(),label.c_str()));
	  delete c;
	  delete latex;


/*

	  if ( LRLabel=="L-R"){
	  // -- draw MeanTimeREF plots:
	  c = new TCanvas(Form("c_MeanTimeREF_%s",label.c_str()),Form("c_MeanTimeREF_%s",label.c_str()));
          gPad -> SetLogy();

          histo = (TH1F*)( inFile->Get(Form("h1_MeanTimeREF_%s",label.c_str())) );
          if( !histo ) continue;
          histo -> SetTitle(";Time [ps];entries");
          histo -> SetLineColor(kRed);
          histo -> SetLineWidth(2);
          histo -> Draw();
	  latex -> Draw("same");
          histo -> Write();
                  
	  c -> Print(Form("%s/MeanTimeREF/c_MeanTimeREF_%s.png",plotDir.c_str(),label.c_str()));        
	  c -> Print(Form("%s/MeanTimeREF/c_MeanTimeREF_%s.pdf",plotDir.c_str(),label.c_str()));

	  delete c;
          delete latex;
	  }
*/



	}// end loop over L, R, L-R labels


      }// -- end loop over bars

    } // -- end loop over stepLabels

  // ---  end 1st plots






  //------------------------
  //--- 2nd loop over events
  std::map<int,std::map<int,bool> > accept;

  for(auto mapIt : trees)
    {
      ModuleEventClass* anEvent = new ModuleEventClass();

      mapIt.second -> SetBranchAddress("event",&anEvent);

      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
	{
	  if( entry%100000 == 0 ) {
	    std::cout << ">>> 2nd loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
	    //TrackProcess(cpu, mem, vsz, rss);
	  }

	  mapIt.second -> GetEntry(entry);

	  bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
	  if (!barFound) continue;


	  int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );


	  int indexREF( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1));

	  accept[index1][entry] = false;

	  if(!ranges["L-R"][index1] ) continue;

	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

	  if( energyBinAverage < 1 ) continue;

	  accept[index1][entry] = true;

	  double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
	  


	  double MCP_phase = anEvent->mcp_phi_peak;
          double BarPhase_L = std::fmod(anEvent->timeL - (6250 * anEvent->t1coarseL), 6250.0);
	  double BarPhase_R = std::fmod(anEvent->timeR - (6250 * anEvent->t1coarseR), 6250.0);
	  long long phaseBarREF_L = std::fmod(anEvent->timeL_ext - (6250 * anEvent->t1coarseL_ext), 6250.0);
          long long phaseBarREF_R = std::fmod(anEvent->timeR_ext - (6250 * anEvent->t1coarseR_ext), 6250.0);
          double mean_phase_bar =  0.5*(BarPhase_L + BarPhase_R);
	  double Delta_phase_meanBar_MCP = mean_phase_bar - MCP_phase;
	  double Delta_BarPhaseL_BarPhaseR = BarPhase_L - BarPhase_R; 

	  if(Debug_on) {
	         std::cout
		 << "\n[DEBUG] =========================================\n"
		 << "  MCP_phase        : " << anEvent->mcp_phi_peak << "\n" 
		 << "  BarPhase_L       : " << BarPhase_L << "\n"
		 << "  BarPhase_R       : " << BarPhase_R << "\n"
		 << "=========================================\n"
                 << std::endl;
	  }

	  if( h1_mean_phase_bar[index2] == NULL )
	    {
	      std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));

	      //t1fine DUT bars histograms:
	      h1_t1fineL[index2] = new TH1F(Form("h1_t1fineL_%s",labelLR_energyBin.c_str()),"",50,0.,1000.);
              h1_t1fineR[index2] = new TH1F(Form("h1_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0.,1000.);
	      h2_t1fineL_vs_t1fineR[index2] = new TH2F(Form("h2_t1fineL_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",25,0.,1000.,25,0.,1000.);
              
	      //phase [ps] DUT bars histograms::
	      h1_phase_barL[index2] = new TH1F(Form("h1_phase_barL_%s",labelLR_energyBin.c_str()),"",250,0.,10000.);
              h1_phase_barR[index2] = new TH1F(Form("h1_phase_barR_%s",labelLR_energyBin.c_str()),"",250,0.,10000.);
              h2_phase_bar_L_vs_R[index2] = new TH2F(Form("h2_phase_bar_L_vs_R_%s",labelLR_energyBin.c_str()),"",100,0.,7000.,100,0.,7000.);


	      h1_mean_phase_bar[index2] = new TH1F(Form("h1_mean_phase_bar_%s",labelLR_energyBin.c_str()),"",50,0.,1000.);
              h1_MPC_phase[index2] = new TH1F(Form("h1_MCP_phase_%s",labelLR_energyBin.c_str()),"",250,-10000.,10000.);
              h2_phase_MeanBar_vs_MCP[index2] = new TH2F(Form("h2_phase_MeanBar_vs_MCP_%s",labelLR_energyBin.c_str()),"",100,0.,7000.,100,0.,7000.);
	  
	      h1_DeltaPhase_meanBar_MCP[index2] = new TH1F(Form("h1_DeltaPhase_meanBar_MCP_%s",labelLR_energyBin.c_str()),"",250,-6500.,6500.);     
              h1_DeltaPhase_Bar_LR[index2] = new TH1F(Form("h1_DeltaPhase_Bar_LR_%s",labelLR_energyBin.c_str()),"",250,-6500.,6500.);
	    }

	  // -- raw Delta T histograms
	  h1_t1fineL[index2]->Fill(anEvent->t1fineL);      
	  h1_t1fineR[index2]->Fill(anEvent->t1fineR);
	  h2_t1fineL_vs_t1fineR[index2]->Fill(anEvent->t1fineL,anEvent->t1fineR);
	  
	  h1_phase_barL[index2] ->Fill(BarPhase_L);
          h1_phase_barR[index2] ->Fill(BarPhase_R);
	  h2_phase_bar_L_vs_R[index2] ->Fill(BarPhase_L,BarPhase_R);
	
	  h1_mean_phase_bar[index2]->Fill(mean_phase_bar);
	  h1_MPC_phase[index2]->Fill(MCP_phase);
	  h2_phase_MeanBar_vs_MCP[index2]->Fill(mean_phase_bar,MCP_phase);
	  
	  h1_DeltaPhase_meanBar_MCP[index2]->Fill(Delta_phase_meanBar_MCP);
          h1_DeltaPhase_Bar_LR[index2]->Fill(Delta_BarPhaseL_BarPhaseR);
	


	  // -- Delta Phase REF bar histograms
          if( h1_DeltaPhase_REF[indexREF] == NULL )
            {
              std::string label_ext(Form("Vov%.2f_th%02d",anEvent->Vov,anEvent->vth1));
              std::string label_ext_global(Form("Vov"));
              h1_DeltaPhase_REF[indexREF] = new TH1F(Form("h1_DeltaPhase_REF_%s",label_ext.c_str()),"",2000,-12000,12000.);
             // h1_DeltaTimeLR_REF[indexREF] = new TH1F(Form("h1_DeltaTimeLR_REF_%s",label_ext.c_str()),"",2000,-12000,12000.);
            }
          h1_DeltaPhase_REF[indexREF] -> Fill(phaseBarREF_L  - phaseBarREF_R);
          //h1_DeltaTimeLR_REF[indexREF] ->Fill(DeltaT_REF);





	  // GLOBAL histogram (no indexREF dependence)
          if( h1_DeltaPhaseLR_REF_global == NULL )
          {
                  //h1_DeltaTimeLR_REF_global = new TH1F("h1_DeltaTimeLR_REF_global",";#DeltaT_{REF} = t_{L} - t_{R} [ps];Entries",2000,-12000,12000.);
                  h1_DeltaPhaseLR_REF_global = new TH1F("h1_DeltaPhaseLR_REF_global",";#Delta #phi_{REF} = #phi_{L} - #phi_{R} [ps];Entries",2000,-12000,12000.);
          }
          //h1_DeltaTimeLR_REF_global->Fill(DeltaT_REF);
          if(phaseBarREF_L!=0 && phaseBarREF_R!=0) {
          h1_DeltaPhaseLR_REF_global->Fill(phaseBarREF_L - phaseBarREF_R);

          }


	  if( h1_NumberEventsBar  == NULL )
          {

	  h1_NumberEventsBar = new TH1F("h1_NumberEventsBar",";DUT bar number;Entries",16,-0.5,15.5);
	  }
	  h1_NumberEventsBar->Fill(anEvent->barID);




	} // end loop over entries
    }



  //------------------
  //--- draw 2nd plots

   TF1* fitFunc_DeltaPhaseLR_REF_global;



   std::string label(Form("global"));
                
   //Delta Phase ChL-ChR REF bar -- Global
   c = new TCanvas(Form("c_DeltaPhaseLR_REF_global_%s",label.c_str()),Form("c_DeltaPhaseLR_REF_global_%s",label.c_str()));
   histo = h1_DeltaPhaseLR_REF_global;
   //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
   histo -> GetXaxis() -> SetRangeUser(-1000,1000);
   histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
   histo -> SetTitle(Form(";#Delta #phi = (#phi_{R} - #phi_{L});entries"));
   histo -> SetLineColor(kRed);
   histo -> SetLineWidth(2);
   histo -> Draw();
   histo -> Write();

   fitFunc_DeltaPhaseLR_REF_global = new TF1(Form("fitFunc_DeltaPhaseLR_REF_global_%s",label.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
   histo -> Fit(fitFunc_DeltaPhaseLR_REF_global,"QNRS");
   histo -> Fit(fitFunc_DeltaPhaseLR_REF_global,"QSR+","",fitFunc_DeltaPhaseLR_REF_global->GetParameter(1)-2.*fitFunc_DeltaPhaseLR_REF_global->GetParameter(2),fitFunc_DeltaPhaseLR_REF_global->GetParameter(1)+2.*fitFunc_DeltaPhaseLR_REF_global->GetParameter(2));
   histo -> Fit(fitFunc_DeltaPhaseLR_REF_global,"QSR+","",fitFunc_DeltaPhaseLR_REF_global->GetParameter(1)-2.*fitFunc_DeltaPhaseLR_REF_global->GetParameter(2),fitFunc_DeltaPhaseLR_REF_global->GetParameter(1)+2.*fitFunc_DeltaPhaseLR_REF_global->GetParameter(2));

   fitFunc_DeltaPhaseLR_REF_global -> SetLineColor(kBlue);
   fitFunc_DeltaPhaseLR_REF_global -> SetLineWidth(2);             
   fitFunc_DeltaPhaseLR_REF_global -> Draw("same");

   latex = new TLatex(0.20,0.85,Form("#splitline{REF Module - Bar 8}""{Mean = %.3f, stdDev = %.3f}", histo->GetMean(), histo->GetRMS()));
   latex -> SetNDC();
   latex -> SetTextFont(42);
   latex -> SetTextSize(0.04);
   latex -> SetTextColor(kBlack);
   latex -> Draw("same");

   c -> Print(Form("%s/Loop2_REF_global/c_DeltaPhaseLR_REF_global_%s.png",plotDir.c_str(),label.c_str()));
   c -> Print(Form("%s/Loop2_REF_global/c_DeltaPhaseLR_REF_global_%s.pdf",plotDir.c_str(),label.c_str()));
   delete latex;
   delete c;



   //h1_NumberEventsBar -- Global
   c = new TCanvas(Form("c_h1_NumberEventsBar_%s",label.c_str()),Form("c_h1_NumberEventsBar_%s",label.c_str()));
   histo = h1_NumberEventsBar;
   //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
   //histo -> GetXaxis() -> SetRangeUser(-1000,1000);
   histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
   histo -> SetTitle(Form(";DUT bar number);entries"));
   histo -> SetLineColor(kRed);
   histo -> SetLineWidth(2);
   histo -> Draw();
   histo -> Write();

   latex = new TLatex(0.20,0.85,Form("#splitline{DUT module}""{Counts per bar}"));
   latex -> SetNDC();
   latex -> SetTextFont(42);
   latex -> SetTextSize(0.04);
   latex -> SetTextColor(kBlack);
   latex -> Draw("same");


   c -> Print(Form("%s/Loop2_NumberEventsBar/c_NumberEventsBar_%s.png",plotDir.c_str(),label.c_str()));
   c -> Print(Form("%s/Loop2_NumberEventsBar/c_NumberEventsBar_%s.pdf",plotDir.c_str(),label.c_str()));
   delete latex;
   delete c;





//--------------------------------------------------------



for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth1 = map_ths[stepLabel];

      for(int iBar = 0; iBar < 16; ++iBar)
        {
          bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
          if (!barFound) continue;

          std::string labelLR(Form("bar%02d_%s",iBar,stepLabel.c_str()));

          int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );

          if( !ranges["L-R"][index1] ) continue;

          int nEnergyBins = ranges["L-R"][index1]->size()-1;

          for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
            {
              //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
              double index2( 10000000*iEnergyBin+index1 );


	      // -- DUT bar t1Fine histograms:
	      if (!h1_t1fineL[index2]) continue;

	      std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));

	      c = new TCanvas(Form("c_t1fineL_%s",labelLR_energyBin.c_str()),Form("c_t1fineL_%s",labelLR_energyBin.c_str()));
	      histo = h1_t1fineL[index2];
	      //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
	      //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
	      histo -> SetTitle(Form(";t1fine_{L};entries"));
	      histo -> SetLineColor(kRed);
	      histo -> SetLineWidth(2);
	      histo -> Draw();
	      histo -> Write();

	      latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kRed);
	      latex -> Draw("same");

	      c -> Print(Form("%s/Loop2_t1fineL/c_t1fineL_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	      c -> Print(Form("%s/Loop2_t1fineL/c_t1fineL_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	      delete latex;
	      delete c;



              if(!h1_t1fineR[index2]) continue;
	      c = new TCanvas(Form("c_t1fineR_%s",labelLR_energyBin.c_str()),Form("c_t1fineR_%s",labelLR_energyBin.c_str()));
              histo = h1_t1fineR[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";t1fine_{R};entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop2_t1fineR/c_t1fineR_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop2_t1fineR/c_t1fineR_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;




              if(!h1_mean_phase_bar[index2]) continue;
	      c = new TCanvas(Form("c_mean_phase_bar_%s",labelLR_energyBin.c_str()),Form("c_mean_phase_bar_%s",labelLR_energyBin.c_str()));
              histo = h1_mean_phase_bar[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";< t1Fine >_{L-R};entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d mean}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop2_mean_phase_bar/c_mean_phase_bar_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop2_mean_phase_bar/c_mean_phase_bar_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;



	      //TH2F phase bar L vs phase bar R:
	      if(!h2_t1fineL_vs_t1fineR[index2]) continue;

               c = new TCanvas(Form("c_t1fineL_vs_t1fineR_%s",labelLR_energyBin.c_str()),Form("c_t1fineL_vs_t1fineR_%s",labelLR_energyBin.c_str()));
               c -> SetGridy();

               h2 = h2_t1fineL_vs_t1fineR[index2];
               //h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) -600., h2 -> GetMean(2)+600);
               h2 -> SetTitle(Form(";t1fine_{L};t1fine_{R}"));
               h2 -> Draw("colz");

               latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
               latex -> SetNDC();
               latex -> SetTextFont(42);
               latex -> SetTextSize(0.04);
               latex -> SetTextColor(kRed);
               latex -> Draw("same");

               c -> Print(Form("%s/Loop2_t1fineL_vs_t1fineR/c_t1fineL_vs_t1fineR_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
               c -> Print(Form("%s/Loop2_t1fineL_vs_t1fineR/c_t1fineL_vs_t1fineR_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
               delete c;
               delete latex;



	      // phase MCP histogram:
              if(!h1_MPC_phase[index2]) continue;
	      c = new TCanvas(Form("c_MPC_phase_%s",labelLR_energyBin.c_str()),Form("c_MPC_phase_%s",labelLR_energyBin.c_str()));
              histo = h1_MPC_phase[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form("; phase MCP;entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop2_MCP_phase/c_MCP_phase_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop2_MCP_phase/c_MCP_phase_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;



	      // TH2F mean phase bar vs MCP phase
	      if(!h2_phase_MeanBar_vs_MCP[index2]) continue;

               c = new TCanvas(Form("c_phase_MeanBar_vs_MCP_%s",labelLR_energyBin.c_str()),Form("c_phase_MeanBar_vs_MCP_%s",labelLR_energyBin.c_str()));
               c -> SetGridy();

               h2 = h2_phase_MeanBar_vs_MCP[index2];
               //h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) -600., h2 -> GetMean(2)+600);
               h2 -> SetTitle(Form(";<#phi_{bar}>; #phi_{MCP}"));
               h2 -> Draw("colz");

               latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
               latex -> SetNDC();
               latex -> SetTextFont(42);
               latex -> SetTextSize(0.04);
               latex -> SetTextColor(kRed);
               latex -> Draw("same");

               c -> Print(Form("%s/Loop2_phase_MeanBar_vs_MCP/c_phase_MeanBar_vs_MCP_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
               c -> Print(Form("%s/Loop2_phase_MeanBar_vs_MCP/c_phase_MeanBar_vs_MCP_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
               delete c;
               delete latex;



	       
	      // DELta phase = (mean phase bar) - (pahse MCP):
              if(!h1_DeltaPhase_meanBar_MCP[index2]) continue;
              c = new TCanvas(Form("c_Delta_phase_meanBar_MCP_%s",labelLR_energyBin.c_str()),Form("c_Delta_phase_meanBar_MCP_%s",labelLR_energyBin.c_str()));
              histo = h1_DeltaPhase_meanBar_MCP[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#Delta #phi = (<#phi_{bar}> - #phi_{MCP});entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

	      latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d - V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}",iBar,Vov,int(vth1), histo->GetMean(), histo->GetRMS()));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop2_DeltaPhase_meanBar_MCP/c_DeltaPhase_meanBar_MCP_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop2_DeltaPhase_meanBar_MCP/c_DeltaPhase_meanBar_MCP_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;	       


	      // DELta phase = (phase bar L) - (phase bar R):
              if(!h1_DeltaPhase_Bar_LR[index2]) continue;
              c = new TCanvas(Form("c_DeltaPhase_Bar_LR_%s",labelLR_energyBin.c_str()),Form("c_DeltaPhase_Bar_LR_%s",labelLR_energyBin.c_str()));
              histo = h1_DeltaPhase_Bar_LR[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#Delta #phi = (#phi_{L} - #phi_{R});entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d - V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}",iBar,Vov,int(vth1), histo->GetMean(), histo->GetRMS()));
	      latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop2_DeltaPhase_Bar_LR/c_DeltaPhase_Bar_LR_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop2_DeltaPhase_Bar_LR/c_DeltaPhase_Bar_LR_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;


              // -- h1_phase_bar_L_v2 histograms:
              if (!h1_phase_barL[index2]) continue;

              c = new TCanvas(Form("c_phase_barL_%s",labelLR_energyBin.c_str()),Form("c_phase_barL_%s",labelLR_energyBin.c_str()));
              histo = h1_phase_barL[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#phi_{L} [ps];entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop2_phase_barL/c_phase_barL_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop2_phase_barL/c_phase_barL_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;


	      // -- h1_phase_bar_R_v2 histograms:
              if (!h1_phase_barR[index2]) continue;

              c = new TCanvas(Form("c_phase_barR_%s",labelLR_energyBin.c_str()),Form("c_phase_barR_%s",labelLR_energyBin.c_str()));
              histo = h1_phase_barR[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#phi_{R} [ps];entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop2_phase_barR/c_phase_barR_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop2_phase_barR/c_phase_barR_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;



	      //TH2F phase bar L vs phase bar R:
              if(!h2_phase_bar_L_vs_R[index2]) continue;

               c = new TCanvas(Form("c_phase_bar_L_vs_R_%s",labelLR_energyBin.c_str()),Form("c_phase_bar_L_vs_R_%s",labelLR_energyBin.c_str()));
               c -> SetGridy();

               h2 = h2_phase_bar_L_vs_R[index2];
               //h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) -600., h2 -> GetMean(2)+600);
               h2 -> SetTitle(Form(";#phi_{L} [ps];#phi_{R} [ps]"));
               h2 -> Draw("colz");

               latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
               latex -> SetNDC();
               latex -> SetTextFont(42);
               latex -> SetTextSize(0.04);
               latex -> SetTextColor(kRed);
               latex -> Draw("same");

               c -> Print(Form("%s/Loop2_phase_bar_L_vs_R/c_phase_bar_L_vs_R_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
               c -> Print(Form("%s/Loop2_phase_bar_L_vs_R/c_phase_bar_L_vs_R_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
               delete c;
               delete latex;



            } // --- end loop over energy bins

        } // --- end loop ober bars



      int indexREF( (10000*(Vov*100.)) + (100*vth1));
      std::string label_ext(Form("barExt_%s",stepLabel.c_str()));

      if(!h1_DeltaPhase_REF[indexREF]) continue;
      c = new TCanvas(Form("c_DeltaPhase_REF_%s",label_ext.c_str()),Form("c_DeltaPhase_REF_%s",label_ext.c_str()));
      histo = h1_DeltaPhase_REF[indexREF];
      histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
      histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
      histo -> SetTitle(Form(";#phi_{L} - #phi_{R} [ps];entries"));
      histo -> SetLineColor(kRed);
      histo -> SetLineWidth(2);
      //c -> SetLogy();
      histo -> Draw();
      histo -> Write();
      latex = new TLatex(0.20,0.85,Form("#splitline{bar 7 REF - V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}", Vov, int(vth1), histo->GetMean(), histo->GetRMS()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kBlack);
      latex -> Draw("same");

      c -> Print(Form("%s/Loop2_DeltaPhase_REF/DeltaPhase_REF_%s.png",plotDir.c_str(),label_ext.c_str()));
      c -> Print(Form("%s/Loop2_DeltaPhase_REF/DeltaPhase_REF_%s.pdf",plotDir.c_str(),label_ext.c_str()));
      delete c;
      delete latex;






    } // --- end loop over stepLabels



  //------------------------
  //--- 3nd loop over events
  for(auto mapIt : trees)
    {
      ModuleEventClass* anEvent = new ModuleEventClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);

      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
        {
          if( entry%100000 == 0 ){
            std::cout << ">>> TIME-WALK loop 3: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
            //TrackProcess(cpu, mem, vsz, rss);
          }
          mapIt.second -> GetEntry(entry);

          bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;
	  if(anEvent->barID > 4) continue; //bar selection

	  int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
          int indexREF( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1));
          int indexTotal( (100000*int(anEvent->Vov*100.))  + anEvent->barID );


          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );


	  double mean_t1fine_bar =  0.5*(anEvent->t1fineL + anEvent->t1fineL);
          double MCP_phase = anEvent->mcp_phi_peak;
          double BarPhase_L = std::fmod(anEvent->timeL - (6250 * anEvent->t1coarseL), 6250.0);
          double BarPhase_R = std::fmod(anEvent->timeR - (6250 * anEvent->t1coarseR), 6250.0);
          double phaseBarREF_L = std::fmod(anEvent->timeL_ext - (6250 * anEvent->t1coarseL_ext), 6250.0);
          double phaseBarREF_R = std::fmod(anEvent->timeR_ext - (6250 * anEvent->t1coarseR_ext), 6250.0);
	  double mean_phase_barDUT = 0.5*(BarPhase_L + BarPhase_R);
          double mean_phase_barREF = 0.5*(phaseBarREF_L + phaseBarREF_R);
          double Delta_phase_meanBar_MCP = MCP_phase - mean_phase_barDUT;
	  double Delta_phase_BarL_MCP = BarPhase_L - MCP_phase;
	  double Delta_phase_BarR_MCP = BarPhase_R - MCP_phase;
	  double mean_energy_barDUT =  0.5*(anEvent->energyL + anEvent->energyR);
	  double mean_energy_barREF =  0.5*(anEvent->energyL_ext + anEvent->energyR_ext);
 


          if( p1_DeltaPhase_vs_MeanEnergyBar[index2] == NULL )
            {
              std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));

              p1_DeltaPhase_vs_MeanEnergyBar[index2] = new TProfile(Form("p1_DeltaPhase_vs_MeanEnergyBar_%s",labelLR_energyBin.c_str()),"",40,ranges["L-R"][index1]->at(0),800.);
              //p1_deltaT_L_raw_vs_energyL[index2]->SetErrorOption("s"); //std DEv
              h2_DeltaPhase_vs_MeanEnergyBar[index2] = new TH2F(Form("h2_DeltaPhase_vs_MeanEnergyBar_%s",labelLR_energyBin.c_str()),"",40,ranges["L-R"][index1]->at(0),800., 350, -3500, 3500.);

	      p1_DeltaPhaseL_vs_EnergyBarL[index2] = new TProfile(Form("p1_DeltaPhaseL_vs_EnergyBarL_%s",labelLR_energyBin.c_str()),"",50,ranges["L-R"][index1]->at(0),800.);
	      h2_DeltaPhaseL_vs_EnergyBarL[index2] = new TH2F(Form("h2_DeltaPhaseL_vs_EnergyBarL_%s",labelLR_energyBin.c_str()),"",50,ranges["L-R"][index1]->at(0),800., 2000, -12000., 12000.);
	      p1_DeltaPhaseR_vs_EnergyBarR[index2] = new TProfile(Form("p1_DeltaPhaseR_vs_EnergyBarR_%s",labelLR_energyBin.c_str()),"",50,ranges["L-R"][index1]->at(0),800.);
              h2_DeltaPhaseR_vs_EnergyBarR[index2] = new TH2F(Form("h2_DeltaPhaseR_vs_EnergyBarR_%s",labelLR_energyBin.c_str()),"",50,ranges["L-R"][index1]->at(0),800., 2000, -12000., 12000.);
              h1_DeltaPhase_meanBar_DUT_REF[index2] = new TH1F(Form("h1_DeltaPhase_meanBar_DUT_REF_%s",labelLR_energyBin.c_str()),"",250,-6500.,6500.);

	      h1_DeltaPhase_MCP_MeanBarDUT[index2] = new TH1F(Form("h1_DeltaPhase_MCP_MeanBarDUT_%s",labelLR_energyBin.c_str()),"",350,-3500.,3500.);

	    }

	            
	  if( h1_DeltaPhase_MCP_MeanBarREF[indexREF] == NULL )
            {

	      std::string label_ext(Form("Vov%.2f_th%02d",anEvent->Vov,anEvent->vth1));
	      h1_DeltaPhase_MCP_MeanBarREF[indexREF] = new TH1F(Form("h1_DeltaPhase_MCP_MeanBarREF_%s",label_ext.c_str()),"",250,-6500.,6500.);

            }

	  

	  
//-------------------------------------------------------------
// TOTAL:

	  std::string label_total(Form("Vov%.2f_Bar%02d",anEvent->Vov,anEvent->barID));
	  if(h1_DeltaPhase_MCP_MeanBar_total[indexTotal] == NULL)
	  { 
		  h1_DeltaPhase_MCP_MeanBar_total[indexTotal] = new TH1F(Form("h1_DeltaPhase_MCP_MeanBar_total_%s",label_total.c_str()),"",350,-3500.,3500.);
		  p1_DeltaPhase_vs_MeanEnergyBar_total[indexTotal] = new TProfile(Form("p1_DeltaPhase_vs_MeanEnergyBar_total_%s",label_total.c_str()),"",20,ranges["L-R"][index1]->at(0),800.);
		  h2_DeltaPhase_vs_MeanEnergyBar_total[indexTotal] = new TH2F(Form("h2_DeltaPhase_vs_MeanEnergyBar_total_%s",label_total.c_str()),"",50,ranges["L-R"][index1]->at(0),800., 250, -7000, 7000.);
		  }

	  if(Delta_phase_meanBar_MCP > -2200. && Delta_phase_meanBar_MCP < -1400.){
	  p1_DeltaPhase_vs_MeanEnergyBar_total[indexTotal] -> Fill( mean_energy_barDUT, Delta_phase_meanBar_MCP );
	  h2_DeltaPhase_vs_MeanEnergyBar_total[indexTotal] -> Fill( mean_energy_barDUT, Delta_phase_meanBar_MCP );
	  h1_DeltaPhase_MCP_MeanBar_total[indexTotal] -> Fill(Delta_phase_meanBar_MCP);
	  }
//-------------------------------------------------------------


	  if(Delta_phase_meanBar_MCP > -2500. && Delta_phase_meanBar_MCP < -1300.){ 
          p1_DeltaPhase_vs_MeanEnergyBar[index2] -> Fill( mean_energy_barDUT, Delta_phase_meanBar_MCP );
	  h2_DeltaPhase_vs_MeanEnergyBar[index2] -> Fill( mean_energy_barDUT, Delta_phase_meanBar_MCP );
	  
	  p1_DeltaPhaseL_vs_EnergyBarL[index2] -> Fill( anEvent->energyL, Delta_phase_BarL_MCP );
	  h2_DeltaPhaseL_vs_EnergyBarL[index2] -> Fill( anEvent->energyL, Delta_phase_BarL_MCP );
	  p1_DeltaPhaseR_vs_EnergyBarR[index2] -> Fill( anEvent->energyR, Delta_phase_BarR_MCP );
	  h2_DeltaPhaseR_vs_EnergyBarR[index2] -> Fill( anEvent->energyR, Delta_phase_BarR_MCP );
 
	  h1_DeltaPhase_meanBar_DUT_REF[index2] -> Fill(mean_phase_barDUT - mean_phase_barREF);
	  h1_DeltaPhase_MCP_MeanBarDUT[index2] -> Fill(MCP_phase - mean_phase_barDUT);
	  h1_DeltaPhase_MCP_MeanBarREF[indexREF] -> Fill(MCP_phase - mean_phase_barREF);

	  }
	}

      std::cout << std::endl;
    }



//-----------------------
  //--- draw 3rd loop plots
std::map<double,TF1*> fitFunc_DeltaPhase_MCP_MeanBar_total;
std::map<double,TF1*> fitFunc_DeltaPhase_vs_MeanEnergyBar_total;
std::map<double,TF1*> fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT;
std::map<double,TF1*> fitFunc_DeltaPhase_MCP_MeanBar_DUT;
  for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth1 = map_ths[stepLabel];

      for(int iBar = 0; iBar < 16; ++iBar){

        bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
        if (!barFound) continue;

        std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
	std::string label_total(Form("Vov%.2f_Bar%02d",Vov,iBar));

	int indexTotal( (100000*int(Vov*100.))  + iBar );
        int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
        if( !ranges["L-R"][index1] ) continue;

        int nEnergyBins = ranges["L-R"][index1]->size()-1;
        
	for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
          {
            //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
            double  index2( 10000000*iEnergyBin+index1 );


            if(!h1_DeltaPhase_meanBar_DUT_REF[index2]) continue;

            std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));

            // -- Profile Delta phase vs mean Energy bar
            c = new TCanvas(Form("c_DeltaPhase_vs_MeanEnergyBar_%s",labelLR_energyBin.c_str()),Form("c_DeltaPhase_vs_MeanEnergyBar_%s",labelLR_energyBin.c_str()));

            prof = p1_DeltaPhase_vs_MeanEnergyBar[index2];
            prof -> SetTitle(Form(";<Energy_{Bar}>;#Delta #phi [ps]"));
            prof -> GetYaxis() -> SetRangeUser(-2200.,-1400.);
            prof -> Draw("");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

            //float fitXMin_L = ranges["L"][index1]->at(0);
            //float fitXMax_L = ranges["L"][index1]->at(1);

            fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT[index2] = new TF1(Form("fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT_%s",labelLR_energyBin.c_str()),"pol3",250.,600.);
            prof -> Fit(fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT[index2],"QRS+");
            fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT[index2] -> SetLineColor(kRed);
            fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT[index2] -> SetLineWidth(2);
            fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT[index2] -> Draw("same");

            c -> Print(Form("%s/Loop3_Profile_DeltaPhase_vs_MeanEnergyBar/c_DeltaPhase_vs_MeanEnergyBar_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop3_Profile_DeltaPhase_vs_MeanEnergyBar/c_DeltaPhase_vs_MeanEnergyBar_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;




            // -- profile + TH2F delta phase vs mean energy bar
            c = new TCanvas(Form("c_scatter_DeltaPhase_vs_MeanEnergyBar_%s",labelLR_energyBin.c_str()),Form("c_scatter_DeltaPhase_vs_MeanEnergyBar_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_DeltaPhase_vs_MeanEnergyBar[index2];
            h2->SetTitle(";<Energy_{Bar}> [a.u.];#Delta #phi");
            h2->GetYaxis()->SetRangeUser(- 2200., -1300.);
            h2->Draw("colz");

            prof = p1_DeltaPhase_vs_MeanEnergyBar[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop3_Scatter_DeltaPhase_vs_MeanEnergyBar/c_scatter_DeltaPhase_vs_MeanEnergyBar_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop3_Scatter_DeltaPhase_vs_MeanEnergyBar/c_scatter_DeltaPhase_vs_MeanEnergyBar_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    //h2_DeltaPhaseL_vs_EnergyBarL

	    c = new TCanvas(Form("c_scatter_DeltaPhaseL_vs_EnergyBarL_%s",labelLR_energyBin.c_str()),Form("c_scatter_DeltaPhaseL_vs_EnergyBarL_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_DeltaPhaseL_vs_EnergyBarL[index2];
            h2->SetTitle(";<Energy_{Bar} L> [a.u.];#Delta #phi");
            //h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_DeltaPhaseL_vs_EnergyBarL[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop3_Scatter_DeltaPhaseL_vs_EnergyBarL/c_scatter_DeltaPhaseL_vs_EnergyBarL_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop3_Scatter_DeltaPhaseL_vs_EnergyBarL/c_scatter_DeltaPhaseL_vs_EnergyBarL_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    //DeltaPhaseR_vs_EnergyBarR
	    c = new TCanvas(Form("c_scatter_DeltaPhaseR_vs_EnergyBarR_%s",labelLR_energyBin.c_str()),Form("c_scatter_DeltaPhaseR_vs_EnergyBarR_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_DeltaPhaseR_vs_EnergyBarR[index2];
            h2->SetTitle(";<Energy_{Bar} R> [a.u.];#Delta #phi");
            //h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_DeltaPhaseR_vs_EnergyBarR[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop3_Scatter_DeltaPhaseR_vs_EnergyBarR/c_scatter_DeltaPhaseR_vs_EnergyBarR_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop3_Scatter_DeltaPhaseR_vs_EnergyBarR/c_scatter_DeltaPhaseR_vs_EnergyBarR_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;





	    // DELta phase = (mean phase bar DUT) - (mean phase REF):
              if(!h1_DeltaPhase_meanBar_DUT_REF[index2]) continue;
              c = new TCanvas(Form("c_DeltaPhase_meanBar_DUT_REF_%s",labelLR_energyBin.c_str()),Form("c_DeltaPhase_meanBar_DUT_REF_%s",labelLR_energyBin.c_str()));
              histo = h1_DeltaPhase_meanBar_DUT_REF[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#Delta #phi = (<#phi_{DUT}> - #phi_{REF});entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d - V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}",iBar,Vov,int(vth1), histo->GetMean(), histo->GetRMS()));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_DeltaPhase_meanBar_DUT_REF/c_DeltaPhase_meanBar_DUT_REF_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop3_DeltaPhase_meanBar_DUT_REF/c_DeltaPhase_meanBar_DUT_REF_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;


	      // Delta phase = (phase MCP) - (mean phase DUT):
              if(!h1_DeltaPhase_MCP_MeanBarDUT[index2]) continue;
              c = new TCanvas(Form("c_DeltaPhase_MCP_MeanBarDUT_%s",labelLR_energyBin.c_str()),Form("c_DeltaPhase_MCP_MeanBarDUT_%s",labelLR_energyBin.c_str()));
              histo = h1_DeltaPhase_MCP_MeanBarDUT[index2];
              histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              
	      
	      fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2] = new TF1(Form("fitFunc_DeltaPhase_MCP_MeanBar_DUT_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
              histo -> Fit(fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2],"QNRS");
              histo -> Fit(fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2],"QSR+","",fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2]->GetParameter(1)-2.*fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2]->GetParameter(2),fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2]->GetParameter(1)+2.*fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2]->GetParameter(2));
              histo -> Fit(fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2],"QSR+","",fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2]->GetParameter(1)-2.*fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2]->GetParameter(2),fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2]->GetParameter(1)+2.*fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2]->GetParameter(2));

              fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2] -> SetLineColor(kBlack);
              fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2] -> SetLineWidth(2);
              fitFunc_DeltaPhase_MCP_MeanBar_DUT[index2] -> Draw("same");
	      
	      histo -> SetTitle(Form(";#Delta #phi = (#phi_{MCP} - <#phi_{DUT}>);entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d - V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}",iBar,Vov,int(vth1), histo->GetMean(), histo->GetRMS()));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_DeltaPhase_MCP_MeanBarDUT/c_DeltaPhase_MCP_MeanBarDUT_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop3_DeltaPhase_MCP_MeanBarDUT/c_DeltaPhase_MCP_MeanBarDUT_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;     



	      //TOTAL--------------------------------------------------------------------------

              // h1_DeltaPhase_MCP_MeanBar_total[indexTotal]
              if(!h1_DeltaPhase_MCP_MeanBar_total[indexTotal]) continue;
              c = new TCanvas(Form("c_DeltaPhase_MCP_MeanBar_total_%s",label_total.c_str()),Form("c_DeltaPhase_MCP_MeanBar_total_%s",label_total.c_str()));
              histo = h1_DeltaPhase_MCP_MeanBar_total[indexTotal];
              histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#Delta #phi = (#phi_{MCP} - <#phi_{DUT}>);entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

	      fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal] = new TF1(Form("fitFunc_DeltaPhase_MCP_MeanBar_total_%s",label_total.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
              histo -> Fit(fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal],"QNRS");
              histo -> Fit(fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal],"QSR+","",fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal]->GetParameter(1)-2.*fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal]->GetParameter(2),fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal]->GetParameter(1)+2.*fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal]->GetParameter(2));
              histo -> Fit(fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal],"QSR+","",fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal]->GetParameter(1)-2.*fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal]->GetParameter(2),fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal]->GetParameter(1)+2.*fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal]->GetParameter(2));

              fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal] -> SetLineColor(kBlack);
              fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal] -> SetLineWidth(2);
              fitFunc_DeltaPhase_MCP_MeanBar_total[indexTotal] -> Draw("same");

              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d - V_{OV} = %.2f V}""{Mean = %.3f, stdDev = %.3f}",iBar,Vov, histo->GetMean(), histo->GetRMS()));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/TOTAL_DeltaPhase_MCP_MeanBar_total/c_DeltaPhase_MCP_MeanBar_total_%s.png",plotDir.c_str(),label_total.c_str()));
              c -> Print(Form("%s/TOTAL_DeltaPhase_MCP_MeanBar_total/c_DeltaPhase_MCP_MeanBar_total_%s.pdf",plotDir.c_str(),label_total.c_str()));
              delete latex;
              delete c;





	    //--p1_DeltaPhase_vs_MeanEnergyBar_total
            c = new TCanvas(Form("c_p1_DeltaPhase_vs_MeanEnergyBar_total_%s",label_total.c_str()),Form("c_p1_DeltaPhase_vs_MeanEnergyBar_total_%s",label.c_str()));

            prof = p1_DeltaPhase_vs_MeanEnergyBar_total[indexTotal];
            prof -> SetTitle(Form(";<Energy_{Bar}> [a.u.]; #phi_{MCP} - <#phi_{Bar}>"));
            prof -> GetYaxis() -> SetRangeUser(-3000.,-1500.);
            prof -> Draw("");

            latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d }{V_{OV} = %.2f V}",iBar,Vov));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kBlack);
            latex -> Draw("same");

            //float fitXMin_L = ranges["L"][index1]->at(0);
            //float fitXMax_L = ranges["L"][index1]->at(1);

            fitFunc_DeltaPhase_vs_MeanEnergyBar_total[indexTotal] = new TF1(Form("fitFunc_DeltaPhase_vs_MeanEnergyBar_total_%s",label_total.c_str()),"pol3",280.,600.);
            prof -> Fit(fitFunc_DeltaPhase_vs_MeanEnergyBar_total[indexTotal],"QRS+");
            fitFunc_DeltaPhase_vs_MeanEnergyBar_total[indexTotal] -> SetLineColor(kRed);
            fitFunc_DeltaPhase_vs_MeanEnergyBar_total[indexTotal] -> SetLineWidth(2);
            fitFunc_DeltaPhase_vs_MeanEnergyBar_total[indexTotal] -> Draw("same");

            c -> Print(Form("%s/TOTAL_Scatter_DeltaPhase_vs_MeanEnergyBar_total/c_p1_DeltaPhase_vs_MeanEnergyBar_total_%s.png",plotDir.c_str(),label_total.c_str()));
            c -> Print(Form("%s/TOTAL_Scatter_DeltaPhase_vs_MeanEnergyBar_total/c_DeltaPhase_vs_MeanEnergyBar_total_%s.pdf",plotDir.c_str(),label_total.c_str()));
            delete latex;
            delete c;







	    //--h2_DeltaPhase_vs_MeanEnergyBar_total
            c = new TCanvas(Form("c_scatter_h2_DeltaPhase_vs_MeanEnergyBar_total_%s",label_total.c_str()),Form("c_scatter_h2_DeltaPhase_vs_MeanEnergyBar_total_%s",label_total.c_str()));
            c->SetGridy();

            h2 = h2_DeltaPhase_vs_MeanEnergyBar_total[indexTotal];
            h2->SetTitle(";<Energy_{Bar}> [a.u.]; #phi_{MCP} - <#phi_{Bar}>");
            //h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_DeltaPhase_vs_MeanEnergyBar_total[indexTotal];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V}",iBar,Vov));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/TOTAL_Scatter_DeltaPhase_vs_MeanEnergyBar_total/c_scatter_DeltaPhase_vs_MeanEnergyBar_total_%s.png",plotDir.c_str(),label_total.c_str()));
            c->Print(Form("%s/TOTAL_Scatter_DeltaPhase_vs_MeanEnergyBar_total/c_scatter_DeltaPhase_vs_MeanEnergyBar_total_%s.pdf",plotDir.c_str(),label_total.c_str()));
            delete latex;
            delete c;






	  }
      }


	      int indexREF( (10000*(Vov*100.)) + (100*vth1));
              std::string label_ext(Form("barExt_%s",stepLabel.c_str()));

	      // Delta phase = (phase MCP) - (mean phase REF):
              if(!h1_DeltaPhase_MCP_MeanBarREF[indexREF]) continue;
              c = new TCanvas(Form("c_DeltaPhase_MCP_MeanBarREF_%s",label_ext.c_str()),Form("c_DeltaPhase_MCP_MeanBarREF_%s",label_ext.c_str()));
              histo = h1_DeltaPhase_MCP_MeanBarREF[indexREF];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#Delta #phi = (#phi_{MCP} - <#phi_{REF}>);entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              latex = new TLatex(0.20,0.85,Form("#splitline{bar 07 REF - V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}",Vov,int(vth1), histo->GetMean(), histo->GetRMS()));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_DeltaPhase_MCP_MeanBarREF/c_DeltaPhase_MCP_MeanBarREF_%s.png",plotDir.c_str(),label_ext.c_str()));
              c -> Print(Form("%s/Loop3_DeltaPhase_MCP_MeanBarREF/c_DeltaPhase_MCP_MeanBarREF_%s.pdf",plotDir.c_str(),label_ext.c_str()));
              delete latex;
              delete c;

          
      
    }









  //------------------------
  //--- 4nd loop over events
  for(auto mapIt : trees)
    {
      ModuleEventClass* anEvent = new ModuleEventClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);

      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
        {
          if( entry%100000 == 0 ){
            std::cout << ">>> TIME-WALK loop 4: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
            //TrackProcess(cpu, mem, vsz, rss);
          }
          mapIt.second -> GetEntry(entry);

          bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;
	  if(anEvent->barID > 3) continue; //bar selection 

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
          int indexREF( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1));
          int indexTotal( (100000*int(anEvent->Vov*100.))  + anEvent->barID );
          

          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );


          double mean_t1fine_bar =  0.5*(anEvent->t1fineL + anEvent->t1fineL);
          double MCP_phase = anEvent->mcp_phi_peak;
          double BarPhase_L = std::fmod(anEvent->timeL - (6250 * anEvent->t1coarseL), 6250.0);
          double BarPhase_R = std::fmod(anEvent->timeR - (6250 * anEvent->t1coarseR), 6250.0);
          double phaseBarREF_L = std::fmod(anEvent->timeL_ext - (6250 * anEvent->t1coarseL_ext), 6250.0);
          double phaseBarREF_R = std::fmod(anEvent->timeR_ext - (6250 * anEvent->t1coarseR_ext), 6250.0);
          double mean_phase_barDUT = 0.5*(BarPhase_L + BarPhase_R);
          double mean_phase_barREF = 0.5*(phaseBarREF_L + phaseBarREF_R);
          double Delta_phase_meanBar_MCP = MCP_phase - mean_phase_barDUT;
          double Delta_phase_BarL_MCP = BarPhase_L - MCP_phase;
          double Delta_phase_BarR_MCP = BarPhase_R - MCP_phase;
          double mean_energy_barDUT =  0.5*(anEvent->energyL + anEvent->energyR);
          double mean_energy_barREF =  0.5*(anEvent->energyL_ext + anEvent->energyR_ext);





	  double TW_corr_meanBar = fitFunc_DeltaPhase_vs_MeanEnergyBar_total[indexTotal]->Eval(mean_energy_barDUT);
          double TW_corr_meanBar_DUT = fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT[index2]->Eval(mean_energy_barDUT);


	 
   	  if( h1_DeltaPhase_MCP_MeanBarDUT_corr[index2] == NULL )
            {
              std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));

              p1_DeltaPhase_vs_MeanEnergyBarDUT_corr[index2] = new TProfile(Form("p1_DeltaPhase_vs_MeanEnergyBarDUT_corr_%s",labelLR_energyBin.c_str()),"",40,ranges["L-R"][index1]->at(0),800.);
              h2_DeltaPhase_vs_MeanEnergyBarDUT_corr[index2] = new TH2F(Form("h2_DeltaPhase_vs_MeanEnergyBarDUT_corr_%s",labelLR_energyBin.c_str()),"",40,ranges["L-R"][index1]->at(0),800., 350, -3500, 3500.);

              h1_DeltaPhase_MCP_MeanBarDUT_corr[index2] = new TH1F(Form("h1_DeltaPhase_MCP_MeanBarDUT_corr_%s",labelLR_energyBin.c_str()),"",350,-3500.,3500.);

            }


	   
 	  if(Delta_phase_meanBar_MCP > -2500. && Delta_phase_meanBar_MCP < -1500.){
          p1_DeltaPhase_vs_MeanEnergyBarDUT_corr[index2] -> Fill( mean_energy_barDUT, Delta_phase_meanBar_MCP - TW_corr_meanBar_DUT);
          h2_DeltaPhase_vs_MeanEnergyBarDUT_corr[index2] -> Fill( mean_energy_barDUT, Delta_phase_meanBar_MCP - TW_corr_meanBar_DUT);
          h1_DeltaPhase_MCP_MeanBarDUT_corr[index2] -> Fill(Delta_phase_meanBar_MCP - TW_corr_meanBar_DUT);

          }


//-------------------------------------------------------------
// TOTAL:

          std::string label_total(Form("Vov%.2f_Bar%02d",anEvent->Vov,anEvent->barID));
          if(h1_DeltaPhase_MCP_MeanBar_total_corr[indexTotal] == NULL)
          {
                  h1_DeltaPhase_MCP_MeanBar_total_corr[indexTotal] = new TH1F(Form("h1_DeltaPhase_MCP_MeanBar_total_corr_%s",label_total.c_str()),"",350,-3500.,3500.);
                  p1_DeltaPhase_vs_MeanEnergyBar_total_corr[indexTotal] = new TProfile(Form("p1_DeltaPhase_vs_MeanEnergyBar_total_corr_%s",label_total.c_str()),"",20,ranges["L-R"][index1]->at(0),800.);
                  h2_DeltaPhase_vs_MeanEnergyBar_total_corr[indexTotal] = new TH2F(Form("h2_DeltaPhase_vs_MeanEnergyBar_total_corr_%s",label_total.c_str()),"",50,ranges["L-R"][index1]->at(0),800., 250, -7000, 7000.);
                  }

          if(Delta_phase_meanBar_MCP > -2200. && Delta_phase_meanBar_MCP < -1400.){
          p1_DeltaPhase_vs_MeanEnergyBar_total_corr[indexTotal] -> Fill( mean_energy_barDUT, Delta_phase_meanBar_MCP );
          h2_DeltaPhase_vs_MeanEnergyBar_total_corr[indexTotal] -> Fill( mean_energy_barDUT, Delta_phase_meanBar_MCP );
          h1_DeltaPhase_MCP_MeanBar_total_corr[indexTotal] -> Fill(Delta_phase_meanBar_MCP - TW_corr_meanBar);
          }
//-------------------------------------------------------------


        }

      std::cout << std::endl;
    }





  //-----------------------
  //--- draw 4rd loop plots

//std::map<double,TF1*> fitFunc_DeltaPhase_vs_MeanEnergyBar_total;
std::map<double,TF1*> fitFunc_DeltaPhase_MCP_MeanBar_total_corr;
std::map<double,TF1*> fitFunc_DeltaPhase_MCP_MeanBar_corr;


  for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth1 = map_ths[stepLabel];

      for(int iBar = 0; iBar < 16; ++iBar){

        bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
        if (!barFound) continue;

        std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
        std::string label_total(Form("Vov%.2f_Bar%02d",Vov,iBar));

        int indexTotal( (100000*int(Vov*100.))  + iBar );
        int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
        if( !ranges["L-R"][index1] ) continue;

        int nEnergyBins = ranges["L-R"][index1]->size()-1;

        for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
          {
            //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
            double  index2( 10000000*iEnergyBin+index1 );


	    std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));


	    if(iBar < 4) {
	                  
	    // h1_DeltaPhase_MCP_MeanBarDUT_corr[index2]
              if(!h1_DeltaPhase_MCP_MeanBarDUT_corr[index2]) continue;
              c = new TCanvas(Form("c_DeltaPhase_MCP_MeanBarDUT_corr_%s",labelLR_energyBin.c_str()),Form("c_DeltaPhase_MCP_MeanBarDUT_corr_%s",labelLR_energyBin.c_str()));
              histo = h1_DeltaPhase_MCP_MeanBarDUT_corr[index2];
              histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#Delta #phi = (#phi_{MCP} - <#phi_{DUT}>);entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

	      
              fitFunc_DeltaPhase_MCP_MeanBar_corr[index2] = new TF1(Form("fitFunc_DeltaPhase_MCP_MeanBar_corr_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
              histo -> Fit(fitFunc_DeltaPhase_MCP_MeanBar_corr[index2],"QNRS");
              histo -> Fit(fitFunc_DeltaPhase_MCP_MeanBar_corr[index2],"QSR+","",fitFunc_DeltaPhase_MCP_MeanBar_corr[index2]->GetParameter(1)-2.*fitFunc_DeltaPhase_MCP_MeanBar_corr[index2]->GetParameter(2),fitFunc_DeltaPhase_MCP_MeanBar_corr[index2]->GetParameter(1)+2.*fitFunc_DeltaPhase_MCP_MeanBar_corr[index2]->GetParameter(2));
              histo -> Fit(fitFunc_DeltaPhase_MCP_MeanBar_corr[index2],"QSR+","",fitFunc_DeltaPhase_MCP_MeanBar_corr[index2]->GetParameter(1)-2.*fitFunc_DeltaPhase_MCP_MeanBar_corr[index2]->GetParameter(2),fitFunc_DeltaPhase_MCP_MeanBar_corr[index2]->GetParameter(1)+2.*fitFunc_DeltaPhase_MCP_MeanBar_corr[index2]->GetParameter(2));

              fitFunc_DeltaPhase_MCP_MeanBar_corr[index2] -> SetLineColor(kBlack);
              fitFunc_DeltaPhase_MCP_MeanBar_corr[index2] -> SetLineWidth(2);
              fitFunc_DeltaPhase_MCP_MeanBar_corr[index2] -> Draw("same");

	       latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d - V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}",iBar,Vov,int(vth1), histo->GetMean(), histo->GetRMS()));
              //latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d - V_{OV} = %.2f V}""{Mean = %.3f, stdDev = %.3f}",iBar,Vov, histo->GetMean(), histo->GetRMS()));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop4_DeltaPhase_MCP_MeanBar_corr/c_DeltaPhase_MCP_MeanBarDUT_corr_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop4_DeltaPhase_MCP_MeanBar_corr/c_DeltaPhase_MCP_MeanBarDUT_corr_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;





	      // -- Profile Delta phase vs mean Energy bar
            c = new TCanvas(Form("c_DeltaPhase_vs_MeanEnergyBarDUT_corr_%s",labelLR_energyBin.c_str()),Form("c_DeltaPhase_vs_MeanEnergyBarDUT_corr_%s",labelLR_energyBin.c_str()));

            prof = p1_DeltaPhase_vs_MeanEnergyBarDUT_corr[index2];
            prof -> SetTitle(Form(";<Energy_{Bar}>;#Delta #phi"));
            prof -> GetYaxis() -> SetRangeUser(-500.,500.);
            prof -> Draw("");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

            //float fitXMin_L = ranges["L"][index1]->at(0);
            //float fitXMax_L = ranges["L"][index1]->at(1);
/*
            fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT[index2] = new TF1(Form("fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT_%s",labelLR_energyBin.c_str()),"pol2",250.,600.);
            prof -> Fit(fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT[index2],"QRS+");
            fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT[index2] -> SetLineColor(kRed);
            fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT[index2] -> SetLineWidth(2);
            fitFunc_DeltaPhase_vs_MeanEnergyBar_DUT[index2] -> Draw("same");
*/
            c -> Print(Form("%s/Loop4_Profile_DeltaPhase_vs_MeanEnergyBar_corr/p1_DeltaPhase_vs_MeanEnergyBarDUT_corr_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop4_Profile_DeltaPhase_vs_MeanEnergyBar_corr/p1_DeltaPhase_vs_MeanEnergyBarDUT_corr_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;






	    // -- profile + TH2F delta phase vs mean energy bar
            c = new TCanvas(Form("c_scatter_DeltaPhase_vs_MeanEnergyBarDUT_corr_%s",labelLR_energyBin.c_str()),Form("c_scatter_DeltaPhase_vs_MeanEnergyBarDUT_corr_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_DeltaPhase_vs_MeanEnergyBarDUT_corr[index2];
            h2->SetTitle(";<Energy_{Bar}> [a.u.];#Delta #phi");
            h2->GetYaxis()->SetRangeUser(-500., 500.);
            h2->Draw("colz");

            prof = p1_DeltaPhase_vs_MeanEnergyBarDUT_corr[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_Scatter_DeltaPhase_vs_MeanEnergyBarDUT_corr/c_scatter_DeltaPhase_vs_MeanEnergyBarDUT_corr_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_Scatter_DeltaPhase_vs_MeanEnergyBarDUT_corr/c_scatter_DeltaPhase_vs_MeanEnergyBarDUT_corr_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;





	    }







              //TOTAL--------------------------------------------------------------------------

              // h1_DeltaPhase_MCP_MeanBar_total_corr[indexTotal]
              if(!h1_DeltaPhase_MCP_MeanBar_total_corr[indexTotal]) continue;
              c = new TCanvas(Form("c_DeltaPhase_MCP_MeanBar_total_corr_%s",label_total.c_str()),Form("c_DeltaPhase_MCP_MeanBar_total_corr_%s",label_total.c_str()));
              histo = h1_DeltaPhase_MCP_MeanBar_total_corr[indexTotal];
              histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#Delta #phi = (#phi_{MCP} - <#phi_{DUT}>);entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();


	      fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal] = new TF1(Form("fitFunc_DeltaPhase_MCP_MeanBar_total_corr_%s",label_total.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
              histo -> Fit(fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal],"QNRS");
              histo -> Fit(fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal],"QSR+","",fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal]->GetParameter(1)-2.*fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal]->GetParameter(2),fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal]->GetParameter(1)+2.*fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal]->GetParameter(2));
              histo -> Fit(fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal],"QSR+","",fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal]->GetParameter(1)-2.*fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal]->GetParameter(2),fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal]->GetParameter(1)+2.*fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal]->GetParameter(2));

              fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal] -> SetLineColor(kBlack);
              fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal] -> SetLineWidth(2);
              fitFunc_DeltaPhase_MCP_MeanBar_total_corr[indexTotal] -> Draw("same");


              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d - V_{OV} = %.2f V}""{Mean = %.3f, stdDev = %.3f}",iBar,Vov, histo->GetMean(), histo->GetRMS()));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop4_TOTAL_DeltaPhase_MCP_MeanBar_total_corr/c_DeltaPhase_MCP_MeanBar_total_corr_%s.png",plotDir.c_str(),label_total.c_str()));
              c -> Print(Form("%s/Loop4_TOTAL_DeltaPhase_MCP_MeanBar_total_corr/c_DeltaPhase_MCP_MeanBar_tota_corrl_%s.pdf",plotDir.c_str(),label_total.c_str()));
              delete latex;
              delete c;



	  }
      }
    }








      int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
