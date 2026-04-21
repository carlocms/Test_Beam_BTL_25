#include "interface/AnalysisUtils.h"
#include "interface/Na22SpectrumAnalyzer.h"
//#include "interface/Na22SpectrumAnalyzerSingleBar.h"
#include "interface/Na22SpectrumAnalyzerSingleBar_TOFHIR2.h"
//#include "interface/Na22SpectrumAnalyzerModule_TOFHIR2.h"
#include "interface/Co60SpectrumAnalyzer_2Peaks.h"
#include "interface/FitUtils.h"
#include "interface/SetTDRStyle.h"
#include "interface/AmplitudeWalk_Utils.h"
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



// ============  ********* MAIN  ************ =============== //
int main(int argc, char** argv)
{
  setTDRStyle();
//  float cpu[2]{0}, mem[2]={0}, vsz[2]={0}, rss[2]={0};

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

  bool optional_plots = false;

  //--- get parameters
  // --- From Loop 1:
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir_TW");
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop1/tot_histo/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop1/energy_TW/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop1/MeanTimeREF/",plotDir.c_str()));
 
  // --- From Loop 2:
  system(Form("mkdir -p %s/Loop2/ToT/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2/deltaT_Raw/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2/energyRatioDUT_LR/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2/deltaT_vs_ToT/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2/EnergyRatio_DUT_REF/Left",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2/EnergyRatio_DUT_REF/Right",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2/deltaT_EneBin/",plotDir.c_str()));

  // --- From Loop 3:
  system(Form("mkdir -p %s/Loop3/DeltaT_raw_vs_EnergyRatioREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3/DeltaT_raw_vs_EnergyRatioREF_scatter/",plotDir.c_str())); 
  system(Form("mkdir -p %s/Loop3/DeltaT_raw_vs_MeanEnergyREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3/EneDUT_vs_MeanEneREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3/DeltaT_Raw_wOffset/",plotDir.c_str()));  
  system(Form("mkdir -p %s/Loop3/MeanfitFunc_deltaT_Raw/",plotDir.c_str()));

  // --- From Loop 4:  
  system(Form("mkdir -p %s/Loop4/TW_EneRatio/DeltaT_EneRatioREF_corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4/TW_EneRatio/Optional/DeltaT_EneRatioREFCorr_vs_energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4/TW_EneRatio/Optional/DeltaT_EneRatioREFCorr_vs_energy_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4/TW_EneRatio/DeltaT_EneRatioREFCorr_vs_energyRatioREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4/TW_EneRatio/DeltaT_EneRatioREFCorr_vs_t1fine/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4/TW_EneRatio/DeltaT_EneRatioREFCorr_vs_t1fine_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4/TW_EneRatio/Optional/DeltaT_EneRatioREFCorr_vs_MeanEneREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4/t1fineL_vs_t1fineR_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4/TW_EneDUT/TW_func_comparison/Side_L/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4/TW_EneDUT/TW_func_comparison/Side_R/",plotDir.c_str()));

  //--Loop 5:
  system(Form("mkdir -p %s/Loop5/TW_EneREF/DeltaT_EneRatioREF_t1fine_Corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5/TW_EneREF/DeltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5/TW_EneREF/DeltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF/",plotDir.c_str()));

  // -- Loop 6:
  system(Form("mkdir -p %s/Loop6_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr/",plotDir.c_str()));
 
  std::vector<std::string> LRLabels;
  LRLabels.push_back("L");
  LRLabels.push_back("R");
  LRLabels.push_back("L-R");

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
  std::string step1FileName= opts.GetOpt<std::string>("Input.step1FileName_TW");
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
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameSingleChannel_TW");
  TFile* outFile = TFile::Open(outFileName.c_str(),"RECREATE");
  outFile->cd();


  //std::map<double,TH1F*> h1_mean_time_ext;
  
  //Loop 2:
  std::map<double,TH1F*> h1_deltaT_R_raw;
  std::map<double,TH1F*> h1_deltaT_L_raw;
  std::map<double,TH1F*> h1_energyRatioL_REF;
  std::map<double,TH1F*> h1_energyRatioR_REF;
  std::map<double,TH1F*> h1_ToT_L_TW;
  std::map<double,TH1F*> h1_ToT_R_TW;
  std::map<double,TH1F*> h1_energyRatioDUT;
  std::map<double,TH2F*> h2_deltaT_vs_ToT_L;
  std::map<double,TH2F*> h2_deltaT_vs_ToT_R;
  std::map< std::pair<float,int>, TGraphErrors* > g_MeanfitFunc_deltaT_Raw_L;
  std::map< std::pair<float,int>, TGraphErrors* > g_MeanfitFunc_deltaT_Raw_R;

  //Loop 3:
  std::map<double,TProfile*> p1_deltaT_L_raw_vs_energyRatioREF;
  std::map<double,TProfile*> p1_deltaT_R_raw_vs_energyRatioREF;
  std::map<double,TH2F*> h2_deltaT_L_raw_vs_energyRatioREF;
  std::map<double,TH2F*> h2_deltaT_R_raw_vs_energyRatioREF;
  std::map<double,TH2F*> h2_deltaT_L_raw_vs_MeanEnergyREF;
  std::map<double,TH2F*> h2_deltaT_R_raw_vs_MeanEnergyREF;
  std::map<double,TProfile*> p1_deltaT_L_raw_vs_MeanEnergyREF;
  std::map<double,TProfile*> p1_deltaT_R_raw_vs_MeanEnergyREF;
  std::map<double,TH2F*> h2_EneDUT_L_vs_MeanEneREF;
  std::map<double,TH2F*> h2_EneDUT_R_vs_MeanEneREF;
  std::map<double,TH1F*> h1_deltaT_L_raw_wOffset;
  std::map<double,TH1F*> h1_deltaT_R_raw_wOffset;

  //Loop 4:
  std::map<double,TH1F*> h1_deltaT_L_EneRatioREFcorr;
  std::map<double,TH1F*> h1_deltaT_R_EneRatioREFcorr;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_L_vs_energyRatioREF;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_R_vs_energyRatioREF;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_L_vs_energyRatioREF;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_R_vs_energyRatioREF;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_L_vs_t1fineL;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_R_vs_t1fineR;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_L_vs_t1fineL;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_R_vs_t1fineR;

  std::map<double,TH2F*> h2_t1fineL_vs_t1fineR;
  std::map<double,TProfile*> p1_t1fineL_vs_t1fineR;
  
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_L_vs_energy;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_R_vs_energy;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_L_vs_energy;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_R_vs_energy;
 
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_L_vs_MeanEneREF;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_R_vs_MeanEneREF;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_L_vs_MeanEneREF;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_R_vs_MeanEneREF;

  //Loop 5:
  std::map<double,TH1F*> h1_deltaT_EneRatioREF_t1fine_Corr_L;
  std::map<double,TH1F*> h1_deltaT_EneRatioREF_t1fine_Corr_R;
  std::map<double,TProfile*> p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF;
  std::map<double,TProfile*> p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF;
  std::map<double,TH2F*> h2_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF;
  std::map<double,TH2F*> h2_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF;

  //Loop 6:
  std::map<double,TH1F*> h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L;
  std::map<double,TH1F*> h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R;


  //Delta T histograms per energy bin:
  std::map<double, bool> energyBinsDefined;
  std::map<double, std::array<float,91>> ene_BIN;
  const int n_BinEne = 90;
  std::map<double,TH1F*> h1_deltaT_L_ext_ENE; 


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
	  c -> Print(Form("%s/Loop1/energy_TW/c_energy_external__Vov%.2f_th%02.0f.png",plotDir.c_str(), Vov, vth1));
	  c -> Print(Form("%s/Loop1/energy_TW/c_energy_external__Vov%.2f_th%02.0f.pdf",plotDir.c_str(), Vov, vth1));
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
		  c -> Print(Form("%s/Loop1/tot_histo/c_tot__%s.png",plotDir.c_str(),label.c_str()));
		  c -> Print(Form("%s/Loop1/tot_histo/c_tot__%s.pdf",plotDir.c_str(),label.c_str()));
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
	  c -> Print(Form("%s/Loop1/energy_TW/c_energy__%s.png",plotDir.c_str(),label.c_str()));
	  c -> Print(Form("%s/Loop1/energy_TW/c_energy__%s.pdf",plotDir.c_str(),label.c_str()));
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

  std::cout << std::endl;
  std::cout << "[DRAW Loop 1] ---> Done! " << std::endl;
  std::cout << std::endl;




  //////////////////////////////
  //--- 2nd loop over events  //
  /////////////////////////////
 
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
	  }

	  mapIt.second -> GetEntry(entry);

	  bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
	  if (!barFound) continue;

	  int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

	  accept[index1][entry] = false;

	  if(!ranges["L-R"][index1] ) continue;

	  //--- SELECTION ON THE MEAN ENERGY REF BAR USED FOR EXTERNAL TIME REFERENCE:
	  if((0.5*(anEvent->energyL_ext + anEvent->energyR_ext)) > 320 || (0.5*(anEvent->energyL_ext + anEvent->energyR_ext)) < 240 ) continue;

	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

	  if( energyBinAverage < 1 ) continue;

	  accept[index1][entry] = true;

	  double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

	  long long timeMean_ext = 0.5 * (anEvent->timeL_ext + anEvent->timeR_ext);
          long long deltaT_L_raw = anEvent->timeL - timeMean_ext;
          long long deltaT_R_raw = anEvent->timeR - timeMean_ext;

	  float energyMeanREF = 0.5 *(anEvent->energyL_ext + anEvent->energyR_ext);

	  std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));

	  if( h1_deltaT_L_raw[index2] == NULL )  
	  {
	      //std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
	      h1_deltaT_L_raw[index2] = new TH1F(Form("h1_deltaT_L_raw_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
              h1_deltaT_R_raw[index2] = new TH1F(Form("h1_deltaT_R_raw_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
              h1_energyRatioDUT[index2] = new TH1F(Form("h1_energyRatioDUT_%s",labelLR_energyBin.c_str()),"",1000,0.,5.);
              h1_energyRatioL_REF[index2] = new TH1F(Form("h1_energyRatioL_REF_%s",labelLR_energyBin.c_str()),"",500,0.,5.);
              h1_energyRatioR_REF[index2] = new TH1F(Form("h1_energyRatioR_REF_%s",labelLR_energyBin.c_str()),"",500,0.,5.);  
	  }

	  // -- raw Delta T histograms
	  h1_deltaT_L_raw[index2]->Fill(deltaT_L_raw);      
	  h1_deltaT_R_raw[index2]->Fill(deltaT_R_raw);

	  // -- energy Ratio histograms
          h1_energyRatioDUT[index2] -> Fill( anEvent->energyR / anEvent->energyL );
          h1_energyRatioL_REF[index2] -> Fill( anEvent->energyL/energyMeanREF );
          h1_energyRatioR_REF[index2] -> Fill( anEvent->energyR/energyMeanREF );


          // -- ToT histograms
	  if(optional_plots)
	  {
	      if(h1_ToT_L_TW[index2] == NULL) {	  
	      h1_ToT_L_TW[index2] = new TH1F(Form("h1_ToT_L_TW_%s",labelLR_energyBin.c_str()),"",2000,-60,60.);
              h1_ToT_R_TW[index2] = new TH1F(Form("h1_ToT_R_TW_%s",labelLR_energyBin.c_str()),"",2000,-60,60.);
              h2_deltaT_vs_ToT_L[index2] = new TH2F(Form("h2_deltaT_vs_ToT_L_%s",labelLR_energyBin.c_str()),"",2000,-10, 20, 100, -12000., 12000.);
              h2_deltaT_vs_ToT_R[index2] = new TH2F(Form("h2_deltaT_vs_ToT_R_%s",labelLR_energyBin.c_str()),"",2000,-10, 20, 100, -12000., 12000.);
	      }
          		  
	  h1_ToT_L_TW[index2] -> Fill(anEvent->totL);
	  h1_ToT_L_TW[index2] -> Fill(anEvent->totR);
	  h2_deltaT_vs_ToT_L[index2] -> Fill( anEvent->totL,deltaT_L_raw);
	  h2_deltaT_vs_ToT_R[index2] -> Fill( anEvent->totR,deltaT_R_raw);
	  }


          //------------------------------------------------------------------------
          //Energy Histo per single energy bin (Used for checks)
	  //------------------------------------------------------------------------
	  if(optional_plots)
	  {

	  if (anEvent->barID == 8 && anEvent->vth1 == 10)
          {
            // Create energy bin boundaries only once for this index2
            if (!energyBinsDefined[index2])
            {
                float ene_min = ranges["L"][index1]->at(0);
                float ene_max = 800.;
                float step = (ene_max - ene_min) / n_BinEne;

                for (int i = 0; i <= n_BinEne; i++)
                    ene_BIN[index2][i] = ene_min + i * step;

                energyBinsDefined[index2] = true;
            }

            // Assign event to a specific energy bin
            for (int i = 0; i < n_BinEne; i++)
            {
                if (anEvent->energyL >= ene_BIN[index2][i] &&
                    anEvent->energyL <  ene_BIN[index2][i+1])
                {
                    // Unique index for this ENE bin
                    double index_ene2 = index2 * 100 + i;

                    // Create histogram for this energy bin
                    if (h1_deltaT_L_ext_ENE[index_ene2] == NULL)
                    {
                        std::string label(Form("bar08_Vov%.2f_th10_ENE%02d",anEvent->Vov, i));

                        h1_deltaT_L_ext_ENE[index_ene2] = new TH1F(Form("h1_deltaT_L_ext_%s", label.c_str()),"",40000,-120000, 120000);

                    }

                    h1_deltaT_L_ext_ENE[index_ene2]->Fill(deltaT_L_raw);
           
	   	}
	    }        
	  }
	  }
       //-------------------------------------------------------------------------------------
       //-------------------------------------------------------------------------------------

	  

	} // end loop over entries
    }


  std::cout << std::endl;
  std::cout << " [Loop 2] ---> Done! " << std::endl;
  std::cout << std::endl;




  //////////////////////////////
  //--- Draw plots loop 2    //
  /////////////////////////////

  std::map<double,TF1*> fitFunc_deltaT_Raw_L;
  std::map<double,TF1*> fitFunc_deltaT_Raw_R;
  std::map<double,TF1*> fitFunc_energyRatioDUT;

  std::map<double,float> CTRMeans_L;
  std::map<double,float> CTRSigmas_L;

  std::map<double,float> CTRMeans_R;
  std::map<double,float> CTRSigmas_R;

  ComputeCTRFromHistoMap(h1_deltaT_L_raw, CTRMeans_L, CTRSigmas_L);
  ComputeCTRFromHistoMap(h1_deltaT_R_raw, CTRMeans_R, CTRSigmas_R);

  std::map<std::pair<double,std::string>, double > offsets_deltaT_raw;
  std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> > mean_fitFunc_deltaT_Raw;


//----------------------------
         
if(optional_plots)
{
  for (auto plot : h1_deltaT_L_ext_ENE) {
    
	  double code = plot.first;
	  TH1F* histo = plot.second;

	  long long keyInt = llround(code);
	  int iBin = keyInt % 100;

	  c = new TCanvas(Form("c_deltaT_L_ext_ENE_bin%02d", iBin),Form("c_deltaT_L_ext_ENE_bin%02d", iBin));
	  histo->GetXaxis()->SetRangeUser(histo->GetMean() - 5.*histo->GetRMS(),histo->GetMean() + 5.*histo->GetRMS());
	  histo->SetMaximum(1.25 * histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean() - 2.*histo->GetRMS(),histo->GetMean() + 2.*histo->GetRMS()))));

          histo->SetTitle(";(Time_{left} - Time_{REF}) [ps]; entries");
          histo->SetLineColor(kRed);
          histo->SetLineWidth(2);
	  //c -> SetLogy();
          histo->Draw();

          histo->Write();

	    latex = new TLatex(0.40,0.85,Form("#splitline{bar 8, Bin %i}{#splitline{Mean = %f}""{std Dev = %.3f}}", iBin, histo->GetMean(), histo->GetRMS()));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

	   /* 
	    std::cout << std::fixed << std::setprecision(3);
	    std::cout
		    << " | bin="  << std::setw(2)  << iBin
		    << " | entries=" << std::setw(6) << histo->GetEntries()
		    << " | mean="    << std::setw(8) << histo->GetMean()
		    << " | err="     << std::setw(8) << histo->GetRMS()
		    << std::endl;
		    */
          c->Print(Form("%s/Loop2_deltaT_EneBin/deltaT_L_Energy_bin%02d.png",plotDir.c_str(), iBin));
          c->Print(Form("%s/Loop2_deltaT_EneBin/deltaT_L_Energy_bin%02d.pdf",plotDir.c_str(), iBin));
	  delete c;
	  delete latex;
  }
}
//-----------------------------------------------

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



	      // -- energy Ratio DUT: [Ene_R/Ene_L]
	      if (!h1_energyRatioDUT[index2]) continue;

	      std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));

              //---------------------------------------------------------------------
              // --- Energy channel R / Energy channel L (DUT)
              //---------------------------------------------------------------------

	      fitFunc_energyRatioDUT[index2] = DrawAndFitHistogram(
				               h1_energyRatioDUT[index2],
				               Form("c_energyRatio_%s", labelLR_energyBin.c_str()),
				               "energy_{right} / energy_{left}",
				               plotDir,"Loop2/energyRatioDUT_LR",
				               iBar,Vov,int(vth1));

              //---------------------------------------------------------------------
              // --- Delta T Raw (DUT - REF)
              //---------------------------------------------------------------------

	      // --- Channel L 
	      fitFunc_deltaT_Raw_L[index2] = DrawAndFitHistogram(
			                     h1_deltaT_L_raw[index2],
		      	                     Form("c_deltaT_L_raw_%s", labelLR_energyBin.c_str()),
	      		                     "#Delta T_{L} Raw [ps]",
      			                     plotDir,"Loop2/deltaT_Raw",
			                     iBar,Vov,int(vth1),"L");



              // --- Channel R              
	      fitFunc_deltaT_Raw_R[index2] = DrawAndFitHistogram(
			                     h1_deltaT_R_raw[index2],
		      	                     Form("c_deltaT_R_raw_%s", labelLR_energyBin.c_str()),
	      		                     "#Delta T_{R} Raw [ps]",
      			                     plotDir,"Loop2/deltaT_Raw",
			                     iBar,Vov,int(vth1),"R");


	       // --- Fill offsets_deltaT_raw map:
	       offsets_deltaT_raw[std::pair(index2,"L")] = fitFunc_deltaT_Raw_L[index2]->GetParameter(1);
	       offsets_deltaT_raw[std::pair(index2,"R")] = fitFunc_deltaT_Raw_R[index2]->GetParameter(1);

	       

              //---------------------------------------------------------------------
              // --- Energy Ratio (DUT / Mean REF)
              //---------------------------------------------------------------------

	      // --- Channel L:
	      DrawHistogram(h1_energyRatioL_REF[index2],
			     Form("c_energyRatioL_REF_%s", labelLR_energyBin.c_str()),
			     "(E_{L}/E_{REF})",
			     plotDir,"Loop2/EnergyRatio_DUT_REF/Left",
			     Form("energyRatioL_REF__%s", labelLR_energyBin.c_str()),
			     BuildLatexSplit(iBar, Vov, vth1, h1_energyRatioL_REF[index2], "L"));

	      // --- Channel R:
              DrawHistogram(h1_energyRatioR_REF[index2],
                            Form("c_energyRatioR_REF_%s", labelLR_energyBin.c_str()),
                            "(E_{R}/E_{REF})",
                            plotDir,"Loop2/EnergyRatio_DUT_REF/Right",
                            Form("energyRatioR_REF__%s", labelLR_energyBin.c_str()),
                            BuildLatexSplit(iBar, Vov, vth1, h1_energyRatioR_REF[index2], "R"));


              //---------------------------------------------------------------------
              // --- Additional plots ToT
              //---------------------------------------------------------------------
	       
	       if(optional_plots){

	       // -- draw ToT channels side L
               DrawHistogram(h1_ToT_L_TW[index2],
                             Form("c_ToT_L_TW_%s", labelLR_energyBin.c_str()),
                             "ToT [ns]",
                             plotDir,"Loop2/ToT",
                             Form("ToT_L_TW_%s", labelLR_energyBin.c_str()),
                             BuildLatexSplit(iBar, Vov, vth1, h1_ToT_L_TW[index2], "L"),
			     -1.5,1.5,true);

	       // -- draw ToT channels side R
               DrawHistogram(h1_ToT_R_TW[index2],
                             Form("c_ToT_R_TW_%s", labelLR_energyBin.c_str()),
                             "ToT [ns]",
                             plotDir,"Loop2/ToT",
                             Form("ToT_R_TW_%s", labelLR_energyBin.c_str()),
                             BuildLatexSplit(iBar, Vov, vth1, h1_ToT_R_TW[index2], "R"),
			     -1.5,1.5,true);



            //-------------------------------------------------------------------------------------
            // --- Delta T Raw vs Tot --> TH2F
            //-------------------------------------------------------------------------------------

            // --- Channel L:
            DrawTH2WithProfile(
                            h2_deltaT_vs_ToT_L[index2],
                            nullptr,
                            Form("c_deltaT_vs_ToT_L_%s", labelLR_energyBin.c_str()),
                            "ToT_{L} [ps]","#Delta T [ps]",
                            plotDir,"Loop2/deltaT_vs_ToT",
                            iBar,Vov,vth1,"L");

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_vs_ToT_R[index2],
                            nullptr,
                            Form("c_deltaT_vs_ToT_R_%s", labelLR_energyBin.c_str()),
                            "ToT_{R} [ps]","#Delta T [ps]",
                            plotDir,"Loop2/deltaT_vs_ToT",
                            iBar,Vov,vth1,"R");

	       }// optional_plots


            } // --- end loop over energy bins

        } // --- end loop ober bars


    } // --- end loop over stepLabels


  std::cout << std::endl;
  std::cout << "[DRAW Loop 2] ---> Done! " << std::endl;
  std::cout << std::endl;






  
  //////////////////////////////
  //--- 3rd loop over events  //
  /////////////////////////////
  
  for(auto mapIt : trees)
    {
      ModuleEventClass* anEvent = new ModuleEventClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);

      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
        {

          if( entry%100000 == 0 ) {
            std::cout << ">>> TIME-WALK loop 3: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
          }

          mapIt.second -> GetEntry(entry);

          bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

          long long timeMean_ext = 0.5 * (anEvent->timeL_ext + anEvent->timeR_ext);
	  auto keyL = std::make_pair(index2,"L");
	  auto keyR = std::make_pair(index2,"R");
          long long deltaT_L_raw = (anEvent->timeL - timeMean_ext) - offsets_deltaT_raw[keyL];
          long long deltaT_R_raw = (anEvent->timeR - timeMean_ext) - offsets_deltaT_raw[keyR];

	  float energyMeanREF = 0.5 *(anEvent->energyL_ext + anEvent->energyR_ext);

	  //float deltaT_L_raw_low = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
	  //float deltaT_L_raw_hig = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
	  //float deltaT_R_raw_low = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          //float deltaT_R_raw_hig = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);

	  float deltaT_L_raw_low = - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_L_raw_hig =  5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);

          float deltaT_R_raw_low = - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_R_raw_hig =  5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);


	  // --- Histograms:
	  std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));

          if(h1_deltaT_L_raw_wOffset[index2] == NULL )
            {

              h1_deltaT_L_raw_wOffset[index2] = new TH1F(Form("h1_deltaT_L_raw_wOffset_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
              h1_deltaT_R_raw_wOffset[index2] = new TH1F(Form("h1_deltaT_R_raw_wOffset_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	
              if(optional_plots){	      
	      h2_deltaT_L_raw_vs_MeanEnergyREF[index2] = new TH2F(Form("h2_deltaT_L_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),"",50,200.,800., 2000, -12000., 12000.);
              h2_deltaT_R_raw_vs_MeanEnergyREF[index2] = new TH2F(Form("h2_deltaT_R_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),"",50,200.,800., 2000, -12000., 12000.);
	      p1_deltaT_L_raw_vs_MeanEnergyREF[index2] = new TProfile(Form("p1_deltaT_L_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
	      p1_deltaT_R_raw_vs_MeanEnergyREF[index2] = new TProfile(Form("p1_deltaT_R_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
	      }

	      h2_EneDUT_L_vs_MeanEneREF[index2] = new TH2F(Form("h2_EneDUT_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",80,ranges["L"][index1]->at(0),800.,80,200.,800.);
              h2_EneDUT_R_vs_MeanEneREF[index2] = new TH2F(Form("h2_EneDUT_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",80,ranges["R"][index1]->at(0),800.,80,200.,800.);

	    }


	  if(p1_deltaT_L_raw_vs_energyRatioREF[index2] == NULL ) 
            {
      	      p1_deltaT_L_raw_vs_energyRatioREF[index2] = new TProfile(Form("p1_deltaT_L_raw_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5.);
              h2_deltaT_L_raw_vs_energyRatioREF[index2] = new TH2F(Form("h2_deltaT_L_raw_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);	   
	      p1_deltaT_R_raw_vs_energyRatioREF[index2] = new TProfile(Form("p1_deltaT_R_raw_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5.);
              h2_deltaT_R_raw_vs_energyRatioREF[index2] = new TH2F(Form("h2_deltaT_R_raw_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);
	    }	      


	  // --- Fill Histograms:
	  if((deltaT_L_raw > deltaT_L_raw_low ) && (deltaT_L_raw < deltaT_L_raw_hig))
	  {
             h1_deltaT_L_raw_wOffset[index2]->Fill(deltaT_L_raw);
	     p1_deltaT_L_raw_vs_energyRatioREF[index2] -> Fill( (anEvent->energyL)/(energyMeanREF), deltaT_L_raw);
	     h2_deltaT_L_raw_vs_energyRatioREF[index2] -> Fill( (anEvent->energyL)/(energyMeanREF), deltaT_L_raw );
	     h2_EneDUT_L_vs_MeanEneREF[index2] -> Fill (anEvent->energyL,energyMeanREF);
	  
	     if(optional_plots){ 
	     p1_deltaT_L_raw_vs_MeanEnergyREF[index2] -> Fill( energyMeanREF, deltaT_L_raw);
             h2_deltaT_L_raw_vs_MeanEnergyREF[index2] -> Fill( energyMeanREF, deltaT_L_raw);
	     }
	  }
	  
	  if((deltaT_R_raw > deltaT_R_raw_low ) && (deltaT_R_raw < deltaT_R_raw_hig))
	  {
	     h1_deltaT_R_raw_wOffset[index2]->Fill(deltaT_R_raw);
	     p1_deltaT_R_raw_vs_energyRatioREF[index2] -> Fill( (anEvent->energyR)/(energyMeanREF),deltaT_R_raw);
	     h2_deltaT_R_raw_vs_energyRatioREF[index2] -> Fill( (anEvent->energyR)/(energyMeanREF),deltaT_R_raw );
	     h2_EneDUT_R_vs_MeanEneREF[index2] -> Fill(anEvent->energyR,energyMeanREF);
	     
	     if(optional_plots){ 
             p1_deltaT_R_raw_vs_MeanEnergyREF[index2] -> Fill( energyMeanREF, deltaT_R_raw);
             h2_deltaT_R_raw_vs_MeanEnergyREF[index2] -> Fill( energyMeanREF, deltaT_R_raw);
	     }
	  }

	 
	}
    }

  std::cout << std::endl;
  std::cout << " [Loop 3] ---> Done! " << std::endl;
  std::cout << std::endl;






  //////////////////////////////
  //--- Draw plots loop 3     //
  //////////////////////////////

  std::map<double,TF1*> fitFunc_deltaT_EneLCorr_L;
  std::map<double,TF1*> fitFunc_deltaT_EneRCorr_R;
  std::map<double,TF1*> fitFunc_deltaT_L_raw_ERatioREF;
  std::map<double,TF1*> fitFunc_deltaT_R_raw_ERatioREF;
  std::map<double,TF1*> fitFunc_deltaT_Raw_wOffset_L;
  std::map<double,TF1*> fitFunc_deltaT_Raw_wOffset_R;

  std::map<double,float> CTRMeans_DeltaT_Raw_Offset_L;
  std::map<double,float> CTRSigmas_DeltaT_Raw_Offset_L;

  std::map<double,float> CTRMeans_DeltaT_Raw_Offset_R;
  std::map<double,float> CTRSigmas_DeltaT_Raw_Offset_R;

  ComputeCTRFromHistoMap(h1_deltaT_L_raw_wOffset, CTRMeans_DeltaT_Raw_Offset_L, CTRSigmas_DeltaT_Raw_Offset_L);
  ComputeCTRFromHistoMap(h1_deltaT_R_raw_wOffset, CTRMeans_DeltaT_Raw_Offset_R, CTRSigmas_DeltaT_Raw_Offset_R);

  for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth1 = map_ths[stepLabel];

      for(int iBar = 0; iBar < 16; ++iBar){

        bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
        if (!barFound) continue;

        std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));

        int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
        if( !ranges["L-R"][index1] ) continue;

        int nEnergyBins = ranges["L-R"][index1]->size()-1;

        for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
          {
            //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
            double  index2( 10000000*iEnergyBin+index1 );
            
	    
	    if(!h1_deltaT_L_raw_wOffset[index2]) continue;

            std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));


              //---------------------------------------------------------------------
              // --- Delta T Raw (DUT - REF) with offset applied
              //---------------------------------------------------------------------

              // --- Channel L 
              fitFunc_deltaT_Raw_wOffset_L[index2] = DrawAndFitHistogram(
                                             h1_deltaT_L_raw_wOffset[index2],
                                             Form("c_deltaT_L_raw_wOffset_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{L} Raw w/ Offset [ps]",
                                             plotDir,"Loop3/DeltaT_Raw_wOffset",
                                             iBar,Vov,int(vth1),"L");



              // --- Channel R              
              fitFunc_deltaT_Raw_wOffset_R[index2] = DrawAndFitHistogram(
                                             h1_deltaT_R_raw_wOffset[index2],
                                             Form("c_deltaT_R_raw_wOffset_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{R} Raw w/ Offset [ps]",
                                             plotDir,"Loop3/DeltaT_Raw_wOffset",
                                             iBar,Vov,int(vth1),"R");


              if (fitFunc_deltaT_Raw_wOffset_L[index2] && fitFunc_deltaT_Raw_wOffset_R[index2]) {
              mean_fitFunc_deltaT_Raw[std::make_tuple(iBar,Vov,vth1,"L")] = std::make_pair(fitFunc_deltaT_Raw_wOffset_L[index2]->GetParameter(1), fitFunc_deltaT_Raw_wOffset_L[index2]->GetParameter(2));
              mean_fitFunc_deltaT_Raw[std::make_tuple(iBar,Vov,vth1,"R")] = std::make_pair(fitFunc_deltaT_Raw_wOffset_R[index2]->GetParameter(1), fitFunc_deltaT_Raw_wOffset_L[index2]->GetParameter(2));
              }



            //---------------------------------------------------------------------
            // --- Delta T Raw vs Energy Ratio [DUT]/[Mean REF]
            //---------------------------------------------------------------------

            // --- Channel L:
            fitFunc_deltaT_L_raw_ERatioREF[index2] = DrawAndFitProfile(
                                                p1_deltaT_L_raw_vs_energyRatioREF[index2],
                                                CTRMeans_L,CTRSigmas_L,ranges,true,
                                                index1,index2,
                                                "pol3",
                                                Form("c_deltaT_L_raw_vs_energyRatioREF_%s", labelLR_energyBin.c_str()),
                                                "E_{L}/E_{REF}","#Delta T_{L} raw [ps]",
                                                plotDir,"Loop3/DeltaT_raw_vs_EnergyRatioREF",
                                                iBar, Vov, int(vth1),"L",
						true, 0.,3.5);
	    

            // --- Channel R:
            fitFunc_deltaT_R_raw_ERatioREF[index2] = DrawAndFitProfile(
                                                p1_deltaT_R_raw_vs_energyRatioREF[index2],
                                                CTRMeans_R,CTRSigmas_R,ranges,true,
                                                index1,index2,
                                                "pol3",
                                                Form("c_deltaT_R_raw_vs_energyRatioREF_%s", labelLR_energyBin.c_str()),
                                                "E_{R}/E_{REF}","#Delta T_{R} raw [ps]",
                                                plotDir,"Loop3/DeltaT_raw_vs_EnergyRatioREF",
                                                iBar, Vov, int(vth1),"R",
						true, 0.,3.5);
	    

            //---------------------------------------------------------------------
            // --- Delta T Raw vs Energy Ratio [DUT]/[Mean REF] TH2F + TProfile
            //---------------------------------------------------------------------

	    // --- Channel L:
            DrawTH2WithProfile(
                            h2_deltaT_L_raw_vs_energyRatioREF[index2],
                            p1_deltaT_L_raw_vs_energyRatioREF[index2],
                            Form("c_deltaT_L_raw_vs_energyRatioREF_scatter_%s", labelLR_energyBin.c_str()),
                            "E_{L}/E_{REF}","#Delta T_{L} [ps]",
                            plotDir,"Loop3/DeltaT_raw_vs_EnergyRatioREF_scatter",
                            iBar,Vov,vth1,"L",
			    true, 0., 3.5);

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_R_raw_vs_energyRatioREF[index2],
                            p1_deltaT_R_raw_vs_energyRatioREF[index2],
                            Form("c_deltaT_R_raw_vs_energyRatioREF_scatter_%s", labelLR_energyBin.c_str()),
                            "E_{R}/E_{REF}","#Delta T_{R} [ps]",
                            plotDir,"Loop3/DeltaT_raw_vs_EnergyRatioREF_scatter",
                            iBar,Vov,vth1,"R",
			    true, 0., 3.5);



            //-------------------------------------------------------------------------------------
            // --- Energy DUT vs Mean Energy REF --> TH2F
            //-------------------------------------------------------------------------------------

            // --- Channel L:
            DrawTH2WithProfile(
                            h2_EneDUT_L_vs_MeanEneREF[index2],
                            nullptr,
                            Form("c_h2_EneDUT_L_vs_MeanEneREF_%s", labelLR_energyBin.c_str()),
                            "Energy^{DUT}_{L} [a.u.]","<E^{REF}> [a.u.]",
                            plotDir,"Loop3/EneDUT_vs_MeanEneREF",
                            iBar,Vov,vth1,"L");

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_EneDUT_R_vs_MeanEneREF[index2],
                            nullptr,
                            Form("c_h2_EneDUT_R_vs_MeanEneREF_%s", labelLR_energyBin.c_str()),
                            "Energy^{DUT}_{R} [a.u.]","<E^{REF}> [a.u.]",
                            plotDir,"Loop3/EneDUT_vs_MeanEneREF",
                            iBar,Vov,vth1,"R");



	    

	    if(optional_plots){
            //-------------------------------------------------------------------------------------
            // --- Delta T Raw vs Mean Energy REF (reference bar for external time) TH2F + TProfile
            //-------------------------------------------------------------------------------------

            // --- Channel L:
	    DrawTH2WithProfile(
			    h2_deltaT_L_raw_vs_MeanEnergyREF[index2],
			    p1_deltaT_L_raw_vs_MeanEnergyREF[index2],
			    Form("c_deltaT_raw_vs_MeanEnergyREF_L_%s", labelLR_energyBin.c_str()),
			    "<E_{REF}>","#Delta T_{L} [ps]",
			    plotDir,"Loop3/DeltaT_raw_vs_MeanEnergyREF",
			    iBar,Vov,vth1,"L");

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_R_raw_vs_MeanEnergyREF[index2],
                            p1_deltaT_R_raw_vs_MeanEnergyREF[index2],
                            Form("c_deltaT_raw_vs_MeanEnergyREF_R_%s", labelLR_energyBin.c_str()),
                            "<E_{REF}>","#Delta T_{R} [ps]",
                            plotDir,"Loop3/DeltaT_raw_vs_MeanEnergyREF",
                            iBar,Vov,vth1,"R");
	    }//optional_plots


          }
      }
    }

  std::cout << std::endl;
  std::cout << "[DRAW Loop 3] ---> Done! " << std::endl;
  std::cout << std::endl;






//-------------------------------
//  Plot of the RESIDUAL factors
//-------------------------------

BuildGraphs(mean_fitFunc_deltaT_Raw,g_MeanfitFunc_deltaT_Raw_L,g_MeanfitFunc_deltaT_Raw_R);
PlotGraphs(g_MeanfitFunc_deltaT_Raw_L,g_MeanfitFunc_deltaT_Raw_R,"Mean_DeltaT_Raw_TGraphError","#mu_{#DeltaT} [ps] ",plotDir + "Loop3/MeanfitFunc_deltaT_Raw",-300., 300.);




  //////////////////////////////
  //--- 4th loop over events  //
  /////////////////////////////

  for(auto mapIt : trees)
    {
      ModuleEventClass* anEvent = new ModuleEventClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);

      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
        {
          if( entry%100000 == 0 ){
            std::cout << ">>> TIME-WALK loop 4: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
          }
          mapIt.second -> GetEntry(entry);

          bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
          //double indexBarID( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) );
	 
	  long long timeMean_ext = 0.5 * (anEvent->timeL_ext + anEvent->timeR_ext);
	  auto keyL = std::make_pair(index2,"L");
          auto keyR = std::make_pair(index2,"R");
          long long deltaT_L_raw = (anEvent->timeL - timeMean_ext) - offsets_deltaT_raw[keyL];
          long long deltaT_R_raw = (anEvent->timeR - timeMean_ext) - offsets_deltaT_raw[keyR];

	  //long long deltaT_L_raw = anEvent->timeL - timeMean_ext;
          //long long deltaT_R_raw = anEvent->timeR - timeMean_ext;
          //float t1fineMean = 0.5 * ( anEvent->t1fineR + anEvent->t1fineL );
          float energyMeanREF = 0.5 *(anEvent->energyL_ext + anEvent->energyR_ext);
 
	  //float deltaT_L_raw_low = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          //float deltaT_L_raw_hig = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
	  //float deltaT_R_raw_low = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          //float deltaT_R_raw_hig = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);	  
	  float deltaT_L_raw_low = - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_L_raw_hig =   5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_R_raw_low = - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_R_raw_hig =   5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);

	  // -- EneRatioREF correction term --> From fitFunc_deltaT_L/R_raw_ERatioREF:
	  float EneRatio_corr_L = fitFunc_deltaT_L_raw_ERatioREF[index2]->Eval(anEvent->energyL/energyMeanREF);
          float EneRatio_corr_R = fitFunc_deltaT_R_raw_ERatioREF[index2]->Eval(anEvent->energyR/energyMeanREF);



	  // --- Histograms:
	  std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));

          if( h1_deltaT_L_EneRatioREFcorr[index2] == NULL )
	  {

	    h1_deltaT_L_EneRatioREFcorr[index2] = new TH1F(Form("h1_deltaT_L_EneRatioREFcorr_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
            h1_deltaT_R_EneRatioREFcorr[index2] = new TH1F(Form("h1_deltaT_R_EneRatioREFcorr_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);

	    p1_deltaT_EneRatioREFCorr_L_vs_energyRatioREF[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_L_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5.);
            p1_deltaT_EneRatioREFCorr_R_vs_energyRatioREF[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_R_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5.);

            h2_deltaT_EneRatioREFCorr_L_vs_energyRatioREF[index2] = new TH2F(Form("h2_deltaT_EneRatioREFCorr_L_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);
            h2_deltaT_EneRatioREFCorr_R_vs_energyRatioREF[index2] = new TH2F(Form("h2_deltaT_EneRatioREFCorr_R_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);

            p1_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_L_vs_t1fineL_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
            p1_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_R_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
	    
	    h2_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2] = new TH2F(Form("h2_deltaT_EneRatioREFCorr_L_vs_t1fineL_%s",labelLR_energyBin.c_str()),"",50,0.,1000., 2000, -12000., 12000.);
	    h2_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2] = new TH2F(Form("h2_deltaT_EneRatioREFCorr_R_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0.,1000., 2000, -12000., 12000.);

            p1_t1fineL_vs_t1fineR[index2] = new TProfile(Form("p1_t1fineL_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0.,1000.);
            h2_t1fineL_vs_t1fineR[index2]= new TH2F(Form("h2_t1fineL_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0.,1000.,50,0.,1000.);
 	    

	    if(optional_plots){
            p1_deltaT_EneRatioREFCorr_L_vs_energy[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_L_vs_energy_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800.);
            p1_deltaT_EneRatioREFCorr_R_vs_energy[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_R_vs_energy_%s",labelLR_energyBin.c_str()),"",50,ranges["R"][index1]->at(0),800.);

            h2_deltaT_EneRatioREFCorr_L_vs_energy[index2] = new TH2F(Form("h2_deltaT_TWCorr_vs_energy_L_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800., 2000, -12000., 12000.);
            h2_deltaT_EneRatioREFCorr_R_vs_energy[index2] = new TH2F(Form("h2_deltaT_TWCorr_vs_energy_R_%s",labelLR_energyBin.c_str()),"",50,ranges["R"][index1]->at(0),800., 2000, -12000., 12000.);
	  
            p1_deltaT_EneRatioREFCorr_L_vs_MeanEneREF[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
            p1_deltaT_EneRatioREFCorr_R_vs_MeanEneREF[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
            h2_deltaT_EneRatioREFCorr_L_vs_MeanEneREF[index2] = new TH2F(Form("h2_deltaT_EneRatioREFCorr_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.,2000, -12000., 12000.);
            h2_deltaT_EneRatioREFCorr_R_vs_MeanEneREF[index2] = new TH2F(Form("h2_deltaT_EneRatioREFCorr_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.,2000, -12000., 12000.);
	    }//optional_plots
	    
	  }

	  
	  if((deltaT_L_raw > deltaT_L_raw_low ) && (deltaT_L_raw < deltaT_L_raw_hig)){
          
	  p1_t1fineL_vs_t1fineR[index2] -> Fill ( anEvent->t1fineL, anEvent->t1fineR );
	  h2_t1fineL_vs_t1fineR[index2] -> Fill(anEvent->t1fineL, anEvent->t1fineR);

	  h1_deltaT_L_EneRatioREFcorr[index2] -> Fill ( deltaT_L_raw - EneRatio_corr_L );
	  p1_deltaT_EneRatioREFCorr_L_vs_energyRatioREF[index2] -> Fill ( anEvent->energyL/energyMeanREF, deltaT_L_raw - EneRatio_corr_L );
	  h2_deltaT_EneRatioREFCorr_L_vs_energyRatioREF[index2] -> Fill ( anEvent->energyL/energyMeanREF, deltaT_L_raw - EneRatio_corr_L );
	  p1_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2] -> Fill(anEvent->t1fineL, deltaT_L_raw - EneRatio_corr_L);
	  h2_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2] -> Fill(anEvent->t1fineL, deltaT_L_raw - EneRatio_corr_L);
	  
	  if (optional_plots){
	  p1_deltaT_EneRatioREFCorr_L_vs_energy[index2] -> Fill ( anEvent->energyL, deltaT_L_raw - EneRatio_corr_L );
          h2_deltaT_EneRatioREFCorr_L_vs_energy[index2] -> Fill ( anEvent->energyL, deltaT_L_raw - EneRatio_corr_L );

          p1_deltaT_EneRatioREFCorr_L_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_L_raw - EneRatio_corr_L);
          h2_deltaT_EneRatioREFCorr_L_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_L_raw - EneRatio_corr_L);	  
	  }
	  
	  }


	  if((deltaT_R_raw > deltaT_R_raw_low ) && (deltaT_R_raw < deltaT_R_raw_hig)){	  
	  
	  h1_deltaT_R_EneRatioREFcorr[index2] -> Fill ( deltaT_R_raw - EneRatio_corr_R );
	  p1_deltaT_EneRatioREFCorr_R_vs_energyRatioREF[index2] -> Fill ( anEvent->energyR/energyMeanREF, deltaT_R_raw - EneRatio_corr_R );
          h2_deltaT_EneRatioREFCorr_R_vs_energyRatioREF[index2] -> Fill ( anEvent->energyR/energyMeanREF, deltaT_R_raw - EneRatio_corr_R );
	  p1_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2] -> Fill(anEvent->t1fineR, deltaT_R_raw - EneRatio_corr_R);
	  h2_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2] -> Fill(anEvent->t1fineR, deltaT_R_raw - EneRatio_corr_R);
	  
	  if (optional_plots){
          p1_deltaT_EneRatioREFCorr_R_vs_energy[index2] -> Fill ( anEvent->energyR, deltaT_R_raw - EneRatio_corr_R );
          h2_deltaT_EneRatioREFCorr_R_vs_energy[index2] -> Fill ( anEvent->energyR, deltaT_R_raw - EneRatio_corr_R );	  
	  
          p1_deltaT_EneRatioREFCorr_R_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_R_raw - EneRatio_corr_R);
          h2_deltaT_EneRatioREFCorr_R_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_R_raw - EneRatio_corr_R);	  
	  }
	 
	  }
	  
	}
    }

  std::cout << std::endl;
  std::cout << "[Loop 4] ---> Done! " << std::endl;
  std::cout << std::endl;





DrawOverlayFitsFromIndices(fitFunc_deltaT_EneLCorr_L,{13001204,13001205,13001206,13001207,13001208,13001209,13001210},"c_overlay_L_selected",plotDir+"/Loop4/TW_EneDUT/TW_func_comparison/Side_L","Side L comparison");
DrawOverlayFitsFromIndices(fitFunc_deltaT_EneRCorr_R,{13001204,13001205,13001206,13001207,13001208,13001209,13001210},"c_overlay_R_selected",plotDir+"Loop4/TW_EneDUT/TW_func_comparison/Side_R","Side R comparison");






  //////////////////////////////
  //--- Draw plots loop 4     //
  //////////////////////////////

  std::map<double,TF1*> fitFunc_deltaT_EneRatioREFcorr_L;
  std::map<double,TF1*> fitFunc_deltaT_EneRatioREFcorr_R;
  std::map<double,TF1*> fitFunc_deltaT_EneRatioREF_TW_L;
  std::map<double,TF1*> fitFunc_deltaT_EneRatioREF_TW_R;
  std::map<double,TF1*> fitFunc_deltaT_EneRCorrL_t1fineL;
  std::map<double,TF1*> fitFunc_deltaT_EneRCorrR_t1fineR;

  std::map<double,float> CTRMeans_ERatioREF_L;
  std::map<double,float> CTRSigmas_ERatioREF_L;
  std::map<double,float> CTRMeans_ERatioREF_R;
  std::map<double,float> CTRSigmas_ERatioREF_R;
  ComputeCTRFromHistoMap(h1_deltaT_L_EneRatioREFcorr, CTRMeans_ERatioREF_L, CTRSigmas_ERatioREF_L);
  ComputeCTRFromHistoMap(h1_deltaT_R_EneRatioREFcorr, CTRMeans_ERatioREF_R, CTRSigmas_ERatioREF_R);


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


              if (!h1_deltaT_L_EneRatioREFcorr[index2]) continue;
              std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));


              //---------------------------------------------------------------------
              // --- TH1F Delta T with Energy Ratio correction
              //---------------------------------------------------------------------

              // --- Channel L: 
              fitFunc_deltaT_EneRatioREFcorr_L[index2] = DrawAndFitHistogram(
                                             h1_deltaT_L_EneRatioREFcorr[index2],
                                             Form("c_h1_deltaT_L_EneRatioREFcorr_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{L} w/ Energy Ratio corr [ps]",
                                             plotDir,"Loop4/TW_EneRatio/DeltaT_EneRatioREF_corr",
                                             iBar,Vov,int(vth1),"L");

              // --- Channel R: 
              fitFunc_deltaT_EneRatioREFcorr_R[index2] = DrawAndFitHistogram(
                                             h1_deltaT_R_EneRatioREFcorr[index2],
                                             Form("c_h1_deltaT_R_EneRatioREFcorr_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{R} w/ Energy Ratio corr [ps]",
                                             plotDir,"Loop4/TW_EneRatio/DeltaT_EneRatioREF_corr",
                                             iBar,Vov,int(vth1),"R");


            //-------------------------------------------------------------------------------------
            // --- TH2F + TProfile: Delta T with Eneregy Ratio correction vs Energy Ratio
            //-------------------------------------------------------------------------------------

            // --- Channel L:
            DrawTH2WithProfile(
                            h2_deltaT_EneRatioREFCorr_L_vs_energyRatioREF[index2],
                            p1_deltaT_EneRatioREFCorr_L_vs_energyRatioREF[index2],
                            Form("c_deltaT_EneRatioREFCorr_L_vs_energyRatioREF_%s", labelLR_energyBin.c_str()),
                            "E_{L}^{DUT}/<E^{REF}> [a.u.]","#Delta T_{L} w/ Ene Ratio corr [ps]",
                            plotDir,"Loop4/TW_EneRatio/DeltaT_EneRatioREFCorr_vs_energyRatioREF",
                            iBar,Vov,vth1,"L",
			    true, 0., 3.5);

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_EneRatioREFCorr_R_vs_energyRatioREF[index2],
                            p1_deltaT_EneRatioREFCorr_R_vs_energyRatioREF[index2],
                            Form("c_deltaT_EneRatioREFCorr_R_vs_energyRatioREF_%s", labelLR_energyBin.c_str()),
                            "E_{R}^{DUT}/<E^{REF}> [a.u.]","#Delta T_{R} w/ Ene Ratio corr [ps]",
                            plotDir,"Loop4/TW_EneRatio/DeltaT_EneRatioREFCorr_vs_energyRatioREF",
                            iBar,Vov,vth1,"R",
			    true, 0., 3.5);


            //---------------------------------------------------------------------
            // --- TProfile Delta T with Energy Ratio correction vs t1Fine
            //---------------------------------------------------------------------

            // --- Channel L:
            fitFunc_deltaT_EneRCorrL_t1fineL[index2] = DrawAndFitProfile(
                                                p1_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2],
                                                CTRMeans_ERatioREF_L,CTRSigmas_ERatioREF_L,ranges,false,
                                                index1,index2,
                                                "pol3",
                                                Form("c_deltaT_EneRatioREFCorr_L_vs_t1fineL_%s", labelLR_energyBin.c_str()),
                                                "t1fine_{L}^{DUT}","#Delta T_{L} w/ ERatio corr [ps]",
                                                plotDir,"Loop4/TW_EneRatio/DeltaT_EneRatioREFCorr_vs_t1fine",
                                                iBar, Vov, int(vth1),"L");


            // --- Channel R:
            fitFunc_deltaT_EneRCorrR_t1fineR[index2] = DrawAndFitProfile(
                                                p1_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2],
                                                CTRMeans_ERatioREF_R,CTRSigmas_ERatioREF_R,ranges,false,
                                                index1,index2,
                                                "pol3",
                                                Form("c_deltaT_EneRatioREFCorr_R_vs_t1fineR_%s", labelLR_energyBin.c_str()),
                                                "t1fine_{R}^{DUT}","#Delta T_{R} w/ ERatio corr [ps]",
                                                plotDir,"Loop4/TW_EneRatio/DeltaT_EneRatioREFCorr_vs_t1fine",
                                                iBar, Vov, int(vth1),"R");


            //-------------------------------------------------------------------------------------
            // --- TH2F + TProfile: Delta T with Energy Ratio correction vs t1Fine DUT
            //-------------------------------------------------------------------------------------

            // --- Channel L:
            DrawTH2WithProfile(
                            h2_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2],
                            p1_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2],
                            Form("c_deltaT_EneRatioREFCorr_L_vs_t1fineL_scatter_%s", labelLR_energyBin.c_str()),
                            "t1fine_{L}^{DUT}","#Delta T_{L} w/ ERatio corr [ps]",
                            plotDir,"Loop4/TW_EneRatio/DeltaT_EneRatioREFCorr_vs_t1fine_scatter",
                            iBar,Vov,vth1,"L");

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2],
                            p1_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2],
                            Form("c_deltaT_EneRatioREFCorr_R_vs_t1fineR_scatter_%s", labelLR_energyBin.c_str()),
                            "t1fine_{R}^{DUT}","#Delta T_{R} w/ ERatio corr [ps]",
                            plotDir,"Loop4/TW_EneRatio/DeltaT_EneRatioREFCorr_vs_t1fine_scatter",
                            iBar,Vov,vth1,"R");


            //-------------------------------------------------------------------------------------
            // --- TH2F + TProfile: f1Fine ch L vs t1Fine ch R (DUT module)
            //-------------------------------------------------------------------------------------

            DrawTH2WithProfile(
                            h2_t1fineL_vs_t1fineR[index2],
                            p1_t1fineL_vs_t1fineR[index2],
                            Form("c_h2_t1fineL_vs_t1fineR_scatter_%s", labelLR_energyBin.c_str()),
                            "t1fine_{L}^{DUT}","t1fine_{R}^{DUT}",
                            plotDir,"Loop4/t1fineL_vs_t1fineR_scatter",
                            iBar,Vov,vth1,"L_R");





            if (optional_plots){

            //---------------------------------------------------------------------
            // --- TProfile Delta T with Energy Ratio correction vs Energy DUT
            //---------------------------------------------------------------------

            // --- Channel L:
            fitFunc_deltaT_EneRatioREF_TW_L[index2] = DrawAndFitProfile(
                                                p1_deltaT_EneRatioREFCorr_L_vs_energy[index2],
                                                CTRMeans_ERatioREF_L,CTRSigmas_ERatioREF_L,ranges,false,
                                                index1,index2,
                                                "pol3",
                                                Form("c_p1_deltaT_EneRatioREFCorr_L_vs_energy_%s", labelLR_energyBin.c_str()),
                                                "E_{L}^{DUT}","#Delta T_{L} w/ ERatio corr [ps]",
                                                plotDir,"Loop4/TW_EneRatio/Optional/DeltaT_EneRatioREFCorr_vs_energy",
                                                iBar, Vov, int(vth1),"L");

            // --- Channel R:
            fitFunc_deltaT_EneRatioREF_TW_R[index2] = DrawAndFitProfile(
                                                p1_deltaT_EneRatioREFCorr_R_vs_energy[index2],
                                                CTRMeans_ERatioREF_R,CTRSigmas_ERatioREF_R,ranges,false,
                                                index1,index2,
                                                "pol3",
                                                Form("c_p1_deltaT_EneRatioREFCorr_R_vs_energy_%s", labelLR_energyBin.c_str()),
                                                "E_{R}^{DUT}","#Delta T_{R} w/ ERatio corr [ps]",
                                                plotDir,"Loop4/TW_EneRatio/Optional/DeltaT_EneRatioREFCorr_vs_energy",
                                                iBar, Vov, int(vth1),"R");



            //-------------------------------------------------------------------------------------
            // --- TH2F + TProfile: Delta T with Energy Ratio correction vs vs Energy DUT
            //-------------------------------------------------------------------------------------

            // --- Channel L:
            DrawTH2WithProfile(
                            h2_deltaT_EneRatioREFCorr_L_vs_energy[index2],
                            p1_deltaT_EneRatioREFCorr_L_vs_energy[index2],
                            Form("c_deltaT_EneRatioREFCorr_L_vs_energy_scatter_%s", labelLR_energyBin.c_str()),
                            "E_{L}^{DUT}","#Delta T_{L} w/ ERatio corr [ps]",
                            plotDir,"Loop4/TW_EneRatio/Optional/DeltaT_EneRatioREFCorr_vs_energy_scatter",
                            iBar,Vov,vth1,"L");

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_EneRatioREFCorr_R_vs_energy[index2],
                            p1_deltaT_EneRatioREFCorr_R_vs_energy[index2],
                            Form("c_deltaT_EneRatioREFCorr_R_vs_energy_scatter_%s", labelLR_energyBin.c_str()),
                            "E_{R}^{DUT}","#Delta T_{R} w/ ERatio corr [ps]",
                            plotDir,"Loop4/TW_EneRatio/Optional/DeltaT_EneRatioREFCorr_vs_energy_scatter",
                            iBar,Vov,vth1,"R");



            //-------------------------------------------------------------------------------------
            // --- TH2F + TProfile: Delta T with Energy Ratio correction vs Mean Energy REF
            //-------------------------------------------------------------------------------------

            // --- Channel L:
            DrawTH2WithProfile(
                            h2_deltaT_EneRatioREFCorr_L_vs_MeanEneREF[index2],
                            p1_deltaT_EneRatioREFCorr_L_vs_MeanEneREF[index2],
                            Form("c_deltaT_EneRatioREFCorr_L_vs_MeanEneREF_%s", labelLR_energyBin.c_str()),
                            "<E^{REF}>","#Delta T_{L} w/ ERatio corr [ps]",
                            plotDir,"Loop4/TW_EneRatio/Optional/DeltaT_EneRatioREFCorr_vs_MeanEneREF",
                            iBar,Vov,vth1,"L");

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_EneRatioREFCorr_R_vs_MeanEneREF[index2],
                            p1_deltaT_EneRatioREFCorr_R_vs_MeanEneREF[index2],
                            Form("c_deltaT_EneRatioREFCorr_R_vs_MeanEneREF_%s", labelLR_energyBin.c_str()),
                            "<E^{REF}>","#Delta T_{R} w/ ERatio corr [ps]",
                            plotDir,"Loop4/TW_EneRatio/Optional/DeltaT_EneRatioREFCorr_vs_MeanEneREF",
                            iBar,Vov,vth1,"R");

	    }// optional_plots

            

	    }
	}

    }


  std::cout << std::endl;
  std::cout << "[DRAW Loop 4] ---> Done! " << std::endl;
  std::cout << std::endl;








  //////////////////////////////
  //--- 5th loop over events  //
  /////////////////////////////

  for(auto mapIt : trees)
    {
      ModuleEventClass* anEvent = new ModuleEventClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);

      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
        {
          if( entry%100000 == 0 ){
            std::cout << ">>> TIME-WALK loop 5: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
            //TrackProcess(cpu, mem, vsz, rss);
          }
          mapIt.second -> GetEntry(entry);

          bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

          long long timeMean_ext = 0.5 * (anEvent->timeL_ext + anEvent->timeR_ext);
	  auto keyL = std::make_pair(index2,"L");
          auto keyR = std::make_pair(index2,"R");
          long long deltaT_L_raw = (anEvent->timeL - timeMean_ext) - offsets_deltaT_raw[keyL];
          long long deltaT_R_raw = (anEvent->timeR - timeMean_ext) - offsets_deltaT_raw[keyR];

	  float energyMeanREF = 0.5 *(anEvent->energyL_ext + anEvent->energyR_ext);

          float deltaT_L_raw_low = - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_L_raw_hig =  5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_R_raw_low = - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_R_raw_hig =  5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);

	  // -- EneRatioREF correction (from delta T raw):
          float EneRatio_corr_L = fitFunc_deltaT_L_raw_ERatioREF[index2]->Eval(anEvent->energyL/energyMeanREF);
          float EneRatio_corr_R = fitFunc_deltaT_R_raw_ERatioREF[index2]->Eval(anEvent->energyR/energyMeanREF);

	  // -- EneRatio REF + t1fine corr
          //float EneR_corr_t1fine_L = fitFunc_deltaT_EneRCorrL_t1fineL[index2]->Eval(anEvent->t1fineL);
          //float EneR_corr_t1fine_R = fitFunc_deltaT_EneRCorrR_t1fineR[index2]->Eval(anEvent->t1fineR);

	  float EneR_corr_t1fine_L = GetProfileBinCorrection(p1_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2],anEvent->t1fineL); //--> corr bin x bin
	  float EneR_corr_t1fine_R = GetProfileBinCorrection(p1_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2],anEvent->t1fineR);


	  // -- Mean t1fine REF
	  float Mean_t1fineREF = 0.5*(anEvent->t1fineL_ext + anEvent->t1fineR_ext);


	  // --- Histograms:
          std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));


          if( h1_deltaT_EneRatioREF_t1fine_Corr_L[index2] == NULL )
          {	  
	  
	  h1_deltaT_EneRatioREF_t1fine_Corr_L[index2] = new TH1F(Form("h1_deltaT_EneRatioREF_t1fine_Corr_L_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
          h1_deltaT_EneRatioREF_t1fine_Corr_R[index2] = new TH1F(Form("h1_deltaT_EneRatioREF_t1fine_Corr_R_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	  
          p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2] = new TProfile(Form("p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
          h2_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2] = new TH2F(Form("h2_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0.,1000., 2000, -12000., 12000.);
	  p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2] = new TProfile(Form("p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
          h2_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2] = new TH2F(Form("h2_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0.,1000., 2000, -12000., 12000.);
	  
	  }
	

	  // --- Fill Histograms:
	  if((deltaT_L_raw > deltaT_L_raw_low ) && (deltaT_L_raw < deltaT_L_raw_hig)){
	  
          h1_deltaT_EneRatioREF_t1fine_Corr_L[index2] -> Fill(deltaT_L_raw - EneRatio_corr_L - EneR_corr_t1fine_L);
	  p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_L_raw - EneRatio_corr_L - EneR_corr_t1fine_L);
	  h2_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_L_raw - EneRatio_corr_L - EneR_corr_t1fine_L);
	  }


	  if((deltaT_R_raw > deltaT_R_raw_low ) && (deltaT_R_raw < deltaT_R_raw_hig)){
          
          h1_deltaT_EneRatioREF_t1fine_Corr_R[index2] -> Fill(deltaT_R_raw - EneRatio_corr_R - EneR_corr_t1fine_R);
          p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_R_raw - EneRatio_corr_R - EneR_corr_t1fine_R);
          h2_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_R_raw - EneRatio_corr_R - EneR_corr_t1fine_R);
	  }

        }
    }

  std::cout << std::endl;
  std::cout << "[Loop 5] ---> Done! " << std::endl;
  std::cout << std::endl;





  //////////////////////////////
  //--- Draw plots 5th loop   //
  /////////////////////////////


    std::map<double,TF1*> fitFunc_deltaT_TWCorr_EneRatioCorr_L;
    std::map<double,TF1*> fitFunc_deltaT_TWCorr_EneRatioCorr_R;
    std::map<double,TF1*> fitFunc_deltaT_EneRatioREF_TW_Corr_L;
    std::map<double,TF1*> fitFunc_deltaT_EneRatioREF_TW_Corr_R;
    std::map<double,TF1*> fitFunc_deltaT_TW_t1fineCorr_L;
    std::map<double,TF1*> fitFunc_deltaT_TW_t1fineCorr_R;
    std::map<double,TF1*> fitFunc_deltaT_EneRatioREF_t1fineCorr_L;
    std::map<double,TF1*> fitFunc_deltaT_EneRatioREF_t1fineCorr_R;
    std::map<double,TF1*> fitFunc_deltaT_EneRatio_t1fineCorr_L;
    std::map<double,TF1*> fitFunc_deltaT_EneRatio_t1fineCorr_R;


    for(auto stepLabel : stepLabels)
    {
        float Vov = map_Vovs[stepLabel];      
        float vth1 = map_ths[stepLabel];
	    
	for(int iBar = 0; iBar < 16; ++iBar){

            bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;        
            if (!barFound) continue;

	    std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
        
	    int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );        
	    if( !ranges["L-R"][index1] ) continue;
        
	    int nEnergyBins = ranges["L-R"][index1]->size()-1;

	    for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
	    {
	        //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
		double  index2( 10000000*iEnergyBin+index1 );
            
                std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));


              //---------------------------------------------------------------------
              // --- TH1F Delta T with Energy Ratio + t1Fine DUT correction
              //---------------------------------------------------------------------

              // --- Channel L:
              fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2] = DrawAndFitHistogram(
                                             h1_deltaT_EneRatioREF_t1fine_Corr_L[index2],
                                             Form("c_h1_deltaT_EneRatioREF_t1fine_Corr_L_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{L} w/ ERatio+t1fine^{DUT} corr [ps]",
                                             plotDir,"Loop5/TW_EneREF/DeltaT_EneRatioREF_t1fine_Corr",
                                             iBar,Vov,int(vth1),"L");


              // --- Channel R:
              fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2] = DrawAndFitHistogram(
                                             h1_deltaT_EneRatioREF_t1fine_Corr_R[index2],
                                             Form("c_h1_deltaT_EneRatioREF_t1fine_Corr_R_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{R} w/ ERatio+t1fine^{DUT} corr [ps]",
                                             plotDir,"Loop5/TW_EneREF/DeltaT_EneRatioREF_t1fine_Corr",
                                             iBar,Vov,int(vth1),"R");


            //-------------------------------------------------------------------------------------
            // --- TH2F + TProfile: Delta T with Energy Ratio + t1fine DUT corr vs Mean t1fine REF
            //-------------------------------------------------------------------------------------

            // --- Channel L:  
            DrawTH2WithProfile(
                            h2_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2],
                            p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2],
                            Form("c_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF_%s", labelLR_energyBin.c_str()),
                            "<t1fine^{REF}>","#Delta T_{L} w/ ERatio+t1fine^{DUT} corr [ps]",
                            plotDir,"Loop5/TW_EneREF/DeltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF_scatter",
                            iBar,Vov,vth1,"L");


            // --- Channel R:  
            DrawTH2WithProfile(
                            h2_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2],
                            p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2],
                            Form("c_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF_%s", labelLR_energyBin.c_str()),
                            "<t1fine^{REF}>","#Delta T_{R} w/ ERatio+t1fine^{DUT} corr [ps]",
                            plotDir,"Loop5/TW_EneREF/DeltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF_scatter",
                            iBar,Vov,vth1,"R");


            //-----------------------------------------------------------------------------------
            // --- TProfile Delta T with Energy Ratio + t1fine DUT correction vs Mean t1fine REF
            //-----------------------------------------------------------------------------------

            // --- Channel L:
            fitFunc_deltaT_EneRatio_t1fineCorr_L[index2] = DrawAndFitProfile(
                                                p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2],
                                                CTRMeans_ERatioREF_L,CTRSigmas_ERatioREF_L,ranges,false,
                                                index1,index2,
                                                "pol3",
                                                Form("c_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF_%s", labelLR_energyBin.c_str()),
                                                "<t1fine^{REF}>","#Delta T_{L} w/ ERatio+t1fine^{DUT} corr [ps]",
                                                plotDir,"Loop5/TW_EneREF/DeltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF",
                                                iBar, Vov, int(vth1),"L");

            // --- Channel R:
            fitFunc_deltaT_EneRatio_t1fineCorr_R[index2] = DrawAndFitProfile(
                                                p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2],
                                                CTRMeans_ERatioREF_R,CTRSigmas_ERatioREF_R,ranges,false,
                                                index1,index2,
                                                "pol3",
                                                Form("c_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF_%s", labelLR_energyBin.c_str()),
                                                "<t1fine^{REF}>","#Delta T_{R} w/ ERatio+t1fine^{DUT} corr [ps]",
                                                plotDir,"Loop5/TW_EneRatio/DeltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF",
                                                iBar, Vov, int(vth1),"R");


	  }
      }
    }



  std::cout << std::endl;
  std::cout << "[Draw Loop 5] ---> Done! " << std::endl;
  std::cout << std::endl;




  //////////////////////////////
  //--- 6th loop over events  //
  /////////////////////////////


  for(auto mapIt : trees)
    {
      ModuleEventClass* anEvent = new ModuleEventClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);

      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
        {
          if( entry%100000 == 0 ){
            std::cout << ">>> AMPL TIME-WALK loop 6: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
          }
          mapIt.second -> GetEntry(entry);

          bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

          long long timeMean_ext = 0.5 * (anEvent->timeL_ext + anEvent->timeR_ext);
          auto keyL = std::make_pair(index2,"L");
          auto keyR = std::make_pair(index2,"R");
          long long deltaT_L_raw = (anEvent->timeL - timeMean_ext) - offsets_deltaT_raw[keyL];
          long long deltaT_R_raw = (anEvent->timeR - timeMean_ext) - offsets_deltaT_raw[keyR];

          float energyMeanREF = 0.5 *(anEvent->energyL_ext + anEvent->energyR_ext);

          float deltaT_L_raw_low = - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_L_raw_hig =  5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_R_raw_low = - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_R_raw_hig =  5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);

          // -- EneRatioREF correction (from delta T raw):
          float EneRatio_corr_L = fitFunc_deltaT_L_raw_ERatioREF[index2]->Eval(anEvent->energyL/energyMeanREF);
          float EneRatio_corr_R = fitFunc_deltaT_R_raw_ERatioREF[index2]->Eval(anEvent->energyR/energyMeanREF);

          // -- EneRatio REF + t1fine corr
          //float EneR_corr_t1fine_L = fitFunc_deltaT_EneRCorrL_t1fineL[index2]->Eval(anEvent->t1fineL);
          //float EneR_corr_t1fine_R = fitFunc_deltaT_EneRCorrR_t1fineR[index2]->Eval(anEvent->t1fineR);
	  float EneR_corr_t1fine_L = GetProfileBinCorrection(p1_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2],anEvent->t1fineL); //--> corr bin x bin
          float EneR_corr_t1fine_R = GetProfileBinCorrection(p1_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2],anEvent->t1fineR);

          // -- Mean t1fine REF
          float Mean_t1fineREF = 0.5*(anEvent->t1fineL_ext + anEvent->t1fineR_ext);
          //float Mean_t1fineREF_corr_L = fitFunc_deltaT_EneRatio_t1fineCorr_L[index2]->Eval(Mean_t1fineREF);
	  //float Mean_t1fineREF_corr_R = fitFunc_deltaT_EneRatio_t1fineCorr_R[index2]->Eval(Mean_t1fineREF);


	  float EneRatioREF_Mean_t1fineREF_corr_L = GetProfileBinCorrection(p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2],Mean_t1fineREF); //--> corr bin x bin
          float EneRatioREF_Mean_t1fineREF_corr_R = GetProfileBinCorrection(p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2],Mean_t1fineREF);


	  // --- Histograms:
          std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));


          // -- delta T histograms
          if( h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2] == NULL )
          {

            // -- EneRatioREF corr + t1fine DUT corr + Mean t1fine REF corr
            h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2] = new TH1F(Form("h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
            h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2] = new TH1F(Form("h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	  }


          if((deltaT_L_raw > deltaT_L_raw_low ) && (deltaT_L_raw < deltaT_L_raw_hig)){
          h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2] -> Fill (deltaT_L_raw - EneRatio_corr_L - EneR_corr_t1fine_L - EneRatioREF_Mean_t1fineREF_corr_L);
	  }


	  if((deltaT_R_raw > deltaT_R_raw_low ) && (deltaT_R_raw < deltaT_R_raw_hig)){
          h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2] -> Fill (deltaT_R_raw - EneRatio_corr_R - EneR_corr_t1fine_R - EneRatioREF_Mean_t1fineREF_corr_R );
	  }

	}
    }

  std::cout << std::endl;
  std::cout << "[Loop 6] ---> Done! " << std::endl;
  std::cout << std::endl;





  //////////////////////////////
  //--- Draw plots 6th loop   //
  /////////////////////////////


    std::map<double,TF1*> fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L;
    std::map<double,TF1*> fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R;

    for(auto stepLabel : stepLabels)
    {
            float Vov = map_Vovs[stepLabel];
            float vth1 = map_ths[stepLabel];

            for(int iBar = 0; iBar < 16; ++iBar){

                    bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
                    if (!barFound) continue;

                    std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));

                    int index1( (10000*int(Vov*100.)) + (100*vth1) + iBar );
                    if( !ranges["L-R"][index1] ) continue;

                    int nEnergyBins = ranges["L-R"][index1]->size()-1;

                    for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
                    {
                            //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
                            double  index2( 10000000*iEnergyBin+index1 );


                            if(!h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2]) continue;
                            std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));



              //---------------------------------------------------------------------
              // --- TH1F Delta T with Energy Ratio + t1Fine DUT + Mean t1fine REF correction
              //---------------------------------------------------------------------

              // --- Channel L:
              TF1* f_L = DrawAndFitHistogram(
       			      h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2],
			      Form("c_h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L_%s", labelLR_energyBin.c_str()),
			      "#Delta T_{L} w/ ERatio+t1fine^{DUT}+<t1fine^{REF}> [ps]",
			      plotDir,"Loop5_deltaT_EneRatioREF_t1fine_Corr",
			      iBar,Vov,int(vth1),"L",
			      kBlue, 3.);
	      delete f_L;


              // --- Channel R:
              TF1* f_R = DrawAndFitHistogram(
       			      h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2],
			      Form("c_h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R_%s", labelLR_energyBin.c_str()),
			      "#Delta T_{R} w/ ERatio+t1fine^{DUT}+<t1fine^{REF}> [ps]",
			      plotDir,"Loop5_deltaT_EneRatioREF_t1fine_Corr",
			      iBar,Vov,int(vth1),"R",
			      kBlue, 3.);
	      delete f_R;

          }
      }
    }
  std::cout << std::endl;
  std::cout << "[Draw Loop 6] ---> Done! " << std::endl;
  std::cout << std::endl;




      int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
