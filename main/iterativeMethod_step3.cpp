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
#include <cmath>

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
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir_iterative_step3");
  std::string refCorrectionFile = opts.GetOpt<std::string>("Input.refCorrectionFile");
  //std::string refCorrectionFile = opts.GetOpt<std::string>("Output.refCorrectionFile");
  //TFile* outCorrFile = new TFile(refCorrectionFile.c_str(),"RECREATE");
  system(Form("mkdir -p %s",plotDir.c_str()));
  
  // --- From Loop 1:
  system(Form("mkdir -p %s/Loop1/tot_histo/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop1/energy_DUT_bar/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop1/energy_REF_bar/",plotDir.c_str()));

  // --- From Loop 5:
  system(Form("mkdir -p %s/Loop5/Optional/ToT/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5/DeltaT_Raw/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5/DeltaT_Raw_Raw/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5/Optional/EnergyRatioDUT_LR/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5/Optional/DeltaT_vs_ToT/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5/Optional/DeltaT_EneBin/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5/DeltaT_Mean_Raw/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5/Sigma_DeltaT_Raw_comparison/",plotDir.c_str()));

  // --- From Loop 6:
  system(Form("mkdir -p %s/Loop6/DeltaT_raw_vs_Energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/DeltaT_raw_vs_Energy_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/DeltaT_raw_vs_MeanEnergyREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/EneDUT_vs_MeanEneREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/DeltaT_Raw_wOffset/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/MeanfitFunc_deltaT_Raw/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/MeanProfile_deltaT_Raw/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/MeanProfileComparison_deltaT_Raw/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/TW_func_comparison/Side_L/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/TW_func_comparison/Side_R/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/Chi2_fitFunc_deltaT_vs_Ene/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/Slope_fitFunc_deltaT_vs_Ene/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/Fit_Function_comparison_differentevth/Side_L/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/Fit_Function_comparison_differentevth/Side_R/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6/Slope_fit_Pol1_deltaT_vs_Ene/",plotDir.c_str()));

  // --- From Loop 7:
  system(Form("mkdir -p %s/Loop7/TW_EneDUT/DeltaT_TW_corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop7/TW_EneDUT/DeltaT_TWCorr_vs_Energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop7/TW_EneDUT/DeltaT_TWCorr_vs_t1fine/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop7/TW_EneDUT/DeltaT_TWCorr_vs_t1fine_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop7/TW_EneDUT/Optional/DeltaT_TWCorr_vs_MeanEneREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop7/t1fineL_vs_t1fineR_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop7/Sigma_deltaT_EneCorr_vs_bar/",plotDir.c_str()));

  // --- From Loop 8:
  system(Form("mkdir -p %s/Loop8/TW_EneDUT/DeltaT_TW_t1fine_Corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop8/TW_EneDUT/DeltaT_TW_t1fine_Corr_vs_Meant1fineREF/",plotDir.c_str()));
  
  // --- From Summary plots:
  system(Form("mkdir -p %s/SummaryPlots/Sigma_DeltatT_vs_Bar/",plotDir.c_str()));
  system(Form("mkdir -p %s/SummaryPlots/Sigma_DeltatT_vs_vth/",plotDir.c_str()));
  
  // --- From Loop 9:
  //system(Form("mkdir -p %s/Loop9/DeltaT_TW_t1fine_Meant1fineREF_Corr/",plotDir.c_str()));
 

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
  std::string step1FileName= opts.GetOpt<std::string>("Input.step1FileName");
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

  // --- Load of the Amplitude walk correction functions (from the iterativeMethod_step2) of the REF bar:
  std::map<int,TF1*> fitFunc_deltaT_EneLCorr_L_ref = LoadTF1Map(refCorrectionFile,"f_twCorr_REF_L");
  std::map<int,TF1*> fitFunc_deltaT_EneRCorr_R_ref = LoadTF1Map(refCorrectionFile,"f_twCorr_REF_R");

  // --- define histograms
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameSingleChannel_step3_iterative");
  TFile* outFile = TFile::Open(outFileName.c_str(),"RECREATE");
  outFile->cd();

  // --- Plots Loop 5:
  std::map<double,TH1F*> h1_deltaT_R_raw;
  std::map<double,TH1F*> h1_deltaT_L_raw;
  std::map<double,TH1F*> h1_deltaT_R_raw_raw;
  std::map<double,TH1F*> h1_deltaT_L_raw_raw;
  std::map<double,TH1F*> h1_energyRatioDUT;
  std::map<double,TH1F*> h1_ToT_L_TW;
  std::map<double,TH1F*> h1_ToT_R_TW;
  std::map<double,TH2F*> h2_deltaT_vs_ToT_L;
  std::map<double,TH2F*> h2_deltaT_vs_ToT_R;
  std::map<double,TH1F*> h1_deltaT_Mean_raw;

  // --- Plots Loop 6:
  std::map<double,TProfile*> p1_deltaT_R_raw_vs_energyR;
  std::map<double,TProfile*> p1_deltaT_L_raw_vs_energyL;
  std::map<double,TH2F*> h2_deltaT_L_raw_vs_energyL;
  std::map<double,TH2F*> h2_deltaT_R_raw_vs_energyR;
  std::map<double,TH2F*> h2_deltaT_L_raw_vs_MeanEnergyREF;
  std::map<double,TH2F*> h2_deltaT_R_raw_vs_MeanEnergyREF;
  std::map<double,TProfile*> p1_deltaT_L_raw_vs_MeanEnergyREF;
  std::map<double,TProfile*> p1_deltaT_R_raw_vs_MeanEnergyREF;
  std::map<double,TH2F*> h2_EneDUT_L_vs_MeanEneREF;
  std::map<double,TH2F*> h2_EneDUT_R_vs_MeanEneREF;
  std::map<double,TH1F*> h1_deltaT_L_raw_wOffset;
  std::map<double,TH1F*> h1_deltaT_R_raw_wOffset;
  std::map< std::pair<float,int>, TGraphErrors* > g_Chi2_fitFunc_deltaT_vs_Ene_L;
  std::map< std::pair<float,int>, TGraphErrors* > g_Chi2_fitFunc_deltaT_vs_Ene_R;
  std::map< std::pair<float,int>, TGraphErrors* > g_Slope_fitFunc_deltaT_vs_Ene_L;
  std::map< std::pair<float,int>, TGraphErrors* > g_Slope_fitFunc_deltaT_vs_Ene_R;
  std::map< std::pair<float,int>, TGraphErrors* > g_Slope_fit_Pol1_deltaT_vs_Ene_L;
  std::map< std::pair<float,int>, TGraphErrors* > g_Slope_fit_Pol1_deltaT_vs_Ene_R;
  std::map< std::pair<float,int>, TGraphErrors* > g_MeanfitFunc_deltaT_Raw_L;
  std::map< std::pair<float,int>, TGraphErrors* > g_MeanfitFunc_deltaT_Raw_R;
  std::map< std::pair<float,int>, TGraphErrors* > g_mean_profile_deltaT_Raw_L;
  std::map< std::pair<float,int>, TGraphErrors* > g_mean_profile_deltaT_Raw_R;
  std::map< std::pair<float,int>, TGraphErrors* > g_mean_profileComparison_deltaT_Raw_L;
  std::map< std::pair<float,int>, TGraphErrors* > g_mean_profileComparison_deltaT_Raw_R;
  
  // --- Plots Loop 7:
  std::map<double,TH1F*> h1_deltaT_L_TWcorr;
  std::map<double,TH1F*> h1_deltaT_R_TWcorr;
  std::map<double,TProfile*> p1_deltaT_TWCorr_vs_energy_L;
  std::map<double,TProfile*> p1_deltaT_TWCorr_vs_energy_R;
  std::map<double,TH2F*> h2_deltaT_TWCorr_vs_energy_L;
  std::map<double,TH2F*> h2_deltaT_TWCorr_vs_energy_R;
  std::map<double,TProfile*> p1_deltaT_TWCorrL_vs_t1fineL;
  std::map<double,TProfile*> p1_deltaT_TWCorrR_vs_t1fineR;
  std::map<double,TH2F*> h2_deltaT_TWCorrL_vs_t1fineL;
  std::map<double,TH2F*> h2_deltaT_TWCorrR_vs_t1fineR;
  std::map<double,TH2F*> h2_t1fineL_vs_t1fineR;
  std::map<double,TProfile*> p1_t1fineL_vs_t1fineR;
  std::map<double,TProfile*> p1_deltaT_TWCorr_L_vs_MeanEneREF;
  std::map<double,TProfile*> p1_deltaT_TWCorr_R_vs_MeanEneREF;
  std::map<double,TH2F*> h2_deltaT_TWCorr_L_vs_MeanEneREF;
  std::map<double,TH2F*> h2_deltaT_TWCorr_R_vs_MeanEneREF;
  std::map< std::pair<float,int>, TGraphErrors* > g_Sigma_deltaT_EneCorr_L;
  std::map< std::pair<float,int>, TGraphErrors* > g_Sigma_deltaT_EneCorr_R;  
 
  // --- Plots Loop 8:
  std::map<double,TH1F*> h1_deltaT_TW_t1fine_Corr_L;
  std::map<double,TH1F*> h1_deltaT_TW_t1fine_Corr_R;
  std::map<double,TProfile*> p1_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF;
  std::map<double,TProfile*> p1_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF;
  std::map<double,TH2F*> h2_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF;
  std::map<double,TH2F*> h2_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF;

/*
  // --- Plots Loop 9:
  std::map<double,TH1F*> h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L;
  std::map<double,TH1F*> h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R;
*/

  //Delta T histograms per energy bin:
  std::map<double, bool> energyBinsDefined;
  std::map<double, std::array<float,91>> ene_BIN;
  const int n_BinEne = 90;
  std::map<double,TH1F*> h1_deltaT_L_ext_ENE; 


  std::map<std::string, std::map<int, std::vector<float>*> > ranges; //ranges[LRlabel][index]
  std::map<std::string, std::map<int,TF1*> > f_landau_ref;
  std::map<std::string, std::map<int,std::vector<float>*>> ranges_ref;
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
  TH1F* histo_L;
  TH1F* histo_R;
  TProfile* prof;
  TH2F* h2;
//=====================================================================================





  //////////////////////////////
  //--- Draw plots loop 1    //
  /////////////////////////////
  
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
      float vth = map_ths[stepLabel];
      std::string VovLabel(Form("Vov%.2f",Vov));
      std::string thLabel(Form("th%02.0f",vth));

      // --------------------------------------------
      // --- external bar mean energy plot:
      // --------------------------------------------
      c = new TCanvas(Form("c_energy_external_Vov%.2f_th%02.0f",Vov,vth), Form("c_energy_external_Vov%.2f_th%02.0f",Vov,vth));
      gPad -> SetLogy();
      histo = (TH1F*)( inFile->Get(Form("h1_energy_external_barL-R_Vov%.2f_th%02.0f", Vov, vth)));
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
	  c -> Print(Form("%s/Loop1/energy_REF_bar/c_energy_external__Vov%.2f_th%02.0f.png",plotDir.c_str(), Vov, vth));
	  c -> Print(Form("%s/Loop1/energy_REF_bar/c_energy_external__Vov%.2f_th%02.0f.pdf",plotDir.c_str(), Vov, vth));
	  delete c;
	  delete ftemp;
	}


      // --------------------------------------------
      // --- REF reference bar single channel energy plot:
      // --------------------------------------------

    c = new TCanvas(Form("c_energy_ref_L_Vov%.2f_th%02.0f",Vov,vth),Form("c_energy_ref_L_Vov%.2f_th%02.0f",Vov,vth));
    gPad->SetLogy();
    histo_L = (TH1F*)( inFile->Get(Form("h1_energy_ref_L_Vov%.2f_th%02.0f", Vov, vth)));

    if(histo_L)
    {
    histo_L->SetTitle(";Energy^{REF}_{L} [a.u.];entries");
    histo_L->SetLineColor(kRed);
    histo_L->SetLineWidth(2);
    histo_L->Draw();
    int index_ref = (10000*int(Vov*100.)) + (100*vth);
    ranges_ref["L"][index_ref] = new std::vector<float>;

    // ---------------------------------------------------
    // --- LANDau fit
    // ---------------------------------------------------

    float max = histo_L->GetBinCenter(histo_L->GetMaximumBin());
    f_landau_ref["L"][index_ref] = new TF1(Form("f_landau_ref_L_Vov%.2f_vth_%02.0f",Vov,vth),"[0]*TMath::Landau(x,[1],[2])",0,1000.);
    float xmin = max * 0.65;
    float xmax = std::min(max*2.5, 940.);
    f_landau_ref["L"][index_ref]->SetRange(xmin,xmax);
    f_landau_ref["L"][index_ref]->SetParameters(histo_L->Integral(histo_L->GetMaximumBin(), histo_L->GetNbinsX())/10,max,0.1*max);
    histo_L->Fit(f_landau_ref["L"][index_ref],"QRS");

    if(f_landau_ref["L"][index_ref]->GetParameter(1) > 0)
    {
        xmin = f_landau_ref["L"][index_ref]->GetParameter(1) - 2.0 * std::abs(f_landau_ref["L"][index_ref]->GetParameter(2));
        xmax = std::min(f_landau_ref["L"][index_ref]->GetParameter(1)*2.0,940.);
        f_landau_ref["L"][index_ref]->SetRange(xmin,xmax);
        histo_L->Fit(f_landau_ref["L"][index_ref],"QRS");
    }

    f_landau_ref["L"][index_ref]->SetLineColor(kBlack);
    f_landau_ref["L"][index_ref]->SetLineWidth(2);
    f_landau_ref["L"][index_ref]->Draw("same");

    // ---------------------------------------------------
    // --- SAVE ENERGY RANGE
    // ---------------------------------------------------

    float eneMin = f_landau_ref["L"][index_ref]->GetParameter(1)- 2.0 * std::abs(f_landau_ref["L"][index_ref]->GetParameter(2));
    float eneMax = std::min(f_landau_ref["L"][index_ref]->GetParameter(1)*2.0,940.);
    ranges_ref["L"][index_ref]->push_back(eneMin);
    ranges_ref["L"][index_ref]->push_back(eneMax);

    // draw cuts
    for(auto range : (*ranges_ref["L"][index_ref]))
    {
        TLine* line = new TLine(range,0.,range,histo_L->GetMaximum());
        line->SetLineWidth(2);
        line->SetLineStyle(7);
        line->Draw("same");
    }

    c->Print(Form("%s/Loop1/energy_REF_bar/c_energy_ref_L_Vov%.2f_th%02.0f.png",plotDir.c_str(), Vov, vth));
    c->Print(Form("%s/Loop1/energy_REF_bar/c_energy_ref_L_Vov%.2f_th%02.0f.pdf",plotDir.c_str(), Vov, vth));
    delete c;
}



// --------------------------------------------
// --- REF reference bar single channel energy plot:
// --------------------------------------------

c = new TCanvas(Form("c_energy_ref_R_Vov%.2f_th%02.0f",Vov,vth),
                Form("c_energy_ref_R_Vov%.2f_th%02.0f",Vov,vth));

gPad->SetLogy();

histo_R = (TH1F*)( inFile->Get(Form("h1_energy_ref_R_Vov%.2f_th%02.0f", Vov, vth)));

if(histo_R)
{
    histo_R->SetTitle(";Energy^{REF}_{R} [a.u.];entries");
    histo_R->SetLineColor(kRed);
    histo_R->SetLineWidth(2);
    histo_R->Draw();

    int index_ref = (10000*int(Vov*100.)) + (100*vth);
    ranges_ref["R"][index_ref] = new std::vector<float>;

    // ---------------------------------------------------
    // --- LANDau fit
    // ---------------------------------------------------

    float max = histo_R->GetBinCenter(histo_R->GetMaximumBin());
    f_landau_ref["R"][index_ref] = new TF1(Form("f_landau_ref_R_Vov%.2f_vth_%02.0f",Vov,vth),"[0]*TMath::Landau(x,[1],[2])",0,1000.);

    float xmin = max * 0.65;
    float xmax = std::min(max*2.5, 940.);

    f_landau_ref["R"][index_ref]->SetRange(xmin,xmax);
    f_landau_ref["R"][index_ref]->SetParameters(histo_R->Integral(histo_R->GetMaximumBin(), histo_R->GetNbinsX())/10,max,0.1*max);
    histo_R->Fit(f_landau_ref["R"][index_ref],"QRS");

    if(f_landau_ref["R"][index_ref]->GetParameter(1) > 0)
    {
        xmin = f_landau_ref["R"][index_ref]->GetParameter(1)- 2.0 * std::abs(f_landau_ref["R"][index_ref]->GetParameter(2));
        xmax = std::min(f_landau_ref["R"][index_ref]->GetParameter(1)*2.0,940.);
        f_landau_ref["R"][index_ref]->SetRange(xmin,xmax);
        histo_R->Fit(f_landau_ref["R"][index_ref],"QRS");
    }

    f_landau_ref["R"][index_ref]->SetLineColor(kBlack);
    f_landau_ref["R"][index_ref]->SetLineWidth(2);
    f_landau_ref["R"][index_ref]->Draw("same");

    // ---------------------------------------------------
    // --- SAVE ENERGY RANGE
    // ---------------------------------------------------
    float eneMin =f_landau_ref["R"][index_ref]->GetParameter(1)- 2.0 * std::abs(f_landau_ref["R"][index_ref]->GetParameter(2));
    float eneMax = std::min(f_landau_ref["R"][index_ref]->GetParameter(1)*2.0,940.);
    ranges_ref["R"][index_ref]->push_back(eneMin);
    ranges_ref["R"][index_ref]->push_back(eneMax);

    // draw cuts

    for(auto range : (*ranges_ref["R"][index_ref]))
    {
        TLine* line = new TLine(range,0.,range,histo_R->GetMaximum());
        line->SetLineWidth(2);
        line->SetLineStyle(7);
        line->Draw("same");
    }

    c->Print(Form("%s/Loop1/energy_REF_bar/c_energy_ref_R_Vov%.2f_th%02.0f.png",plotDir.c_str(), Vov, vth));
    c->Print(Form("%s/Loop1/energy_REF_bar/c_energy_ref_R_Vov%.2f_th%02.0f.pdf",plotDir.c_str(), Vov, vth));
    delete c;
}



      //--------------------------------------------------------
      // --- loop over bars
      for(int iBar = 0; iBar < 16; ++iBar) {

	bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
	if (!barFound) continue;

	int index( (10000*int(Vov*100.)) + (100*vth) + iBar );

	// -- loop over L, R, LR
	for(auto LRLabel : LRLabels ) {

	  //label histo
	  std::string label(Form("bar%02d%s_%s",iBar,LRLabel.c_str(),stepLabel.c_str()));

	  latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d%s}{V_{OV} = %.2f V, th. = %d DAC}",iBar,LRLabel.c_str(),Vov,int(vth)));
	  if (LRLabel == "L-R") {
	    latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth)));
	  }
	  latex -> SetNDC();
	  latex -> SetTextFont(42);
	  latex -> SetTextSize(0.04);
	  latex -> SetTextColor(kRed);

	  // -- draw ToT histograms only for L, R
	  if (LRLabel == "R" || LRLabel == "L")
	    {
	      
	
	      if (optional_plots){	    
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
	      }//optional_plots
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
	    
	    float max = histo->GetBinCenter(histo->GetMaximumBin());
	    //histo->GetXaxis()->SetRangeUser(0,1024);
	    histo->GetXaxis()->SetRangeUser(minE[std::make_pair(iBar, Vov)], 950); // minE to avoid fitting noise

	    f_gaus[index] = new TF1(Form("fit_energy_bar%02d%s_Vov%.2f_vth_%02.0f",iBar,LRLabel.c_str(),Vov,vth), "gaus", max-50, max+50);
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

	    f_landau[index] = new TF1(Form("f_landau_bar%02d%s_Vov%.2f_vth_%02.0f", iBar,LRLabel.c_str(),Vov,vth),"[0]*TMath::Landau(x,[1],[2])", 0,1000.);
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

	    if ( LRLabel=="L-R" && int(vth)==10) std::cout << iBar << "  " << Vov  << "  " << ranges[LRLabel][index] ->at(0) <<std::endl;


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
	  c -> Print(Form("%s/Loop1/energy_DUT_bar/c_energy__%s.png",plotDir.c_str(),label.c_str()));
	  c -> Print(Form("%s/Loop1/energy_DUT_bar/c_energy__%s.pdf",plotDir.c_str(),label.c_str()));
	  delete c;
	  delete latex;


	}// end loop over L, R, L-R labels


      }// -- end loop over bars

    } // -- end loop over stepLabels

  // ---  end 1st plots

  std::cout << std::endl;
  std::cout << "[DRAW Loop 1] ---> Done! " << std::endl;
  std::cout << std::endl;


  std::cout << std::endl;
  std::cout << "====================================================================================================" << std::endl;
  std::cout << " START OF THE COMPLETE ANALYSIS OF THE DUT MODULE ==> USE OF THE CALIBRATED REF BAR FROM THE STEP 2 " << std::endl;
  std::cout << "====================================================================================================" << std::endl;
  std::cout << std::endl;

 
  
// std::map<int,TF1*> fitFunc_deltaT_EneLCorr_L_ref = LoadTF1Map(refCorrectionFile,"f_twCorr_REF_L");
// std::map<int,TF1*> fitFunc_deltaT_EneRCorr_R_ref = LoadTF1Map(refCorrectionFile,"f_twCorr_REF_R");
 
// outFile->cd();




  //////////////////////////////
  //--- 5th loop over events  //
  /////////////////////////////

  std::map<int,std::map<int,bool> > accept;

  for(auto mapIt : trees)
    {
      ModuleEventWithRefClass* anEvent = new ModuleEventWithRefClass();

      mapIt.second -> SetBranchAddress("event",&anEvent);


      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
        {
          if( entry%100000 == 0 ) {
            std::cout << ">>> 5nd loop: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
          }

          mapIt.second -> GetEntry(entry);

          bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
          int index_ref = (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth);
          accept[index1][entry] = false;

          if(!ranges["L-R"][index1] ) continue;

          //--- SELECTION ON THE MEAN ENERGY REF BAR USED FOR EXTERNAL TIME REFERENCE:
          //if((0.5*(anEvent->energyL_ref + anEvent->energyR_ref)) > 300 || (0.5*(anEvent->energyL_ref + anEvent->energyR_ref)) < 230 ) continue;

	  if (anEvent->vth == 40 ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          if( energyBinAverage < 1 ) continue;

          accept[index1][entry] = true;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );

	  long long TW_corr_REF_L = fitFunc_deltaT_EneLCorr_L_ref[index_ref]->Eval(anEvent->energyL_ref);
          long long TW_corr_REF_R = fitFunc_deltaT_EneRCorr_R_ref[index_ref]->Eval(anEvent->energyR_ref);
	  
	  long long timeL_ref_corr = anEvent->timeL_ref - TW_corr_REF_L;
          long long timeR_ref_corr = anEvent->timeR_ref - TW_corr_REF_R;
          
	  long long timeMean_ext_corr = 0.5 * (timeL_ref_corr + timeR_ref_corr);
          long long timeMean_ext_raw = 0.5 * (anEvent->timeL_ref + anEvent->timeR_ref);

          long long deltaT_L_raw = (anEvent->timeL - timeMean_ext_corr);
          long long deltaT_R_raw = (anEvent->timeR - timeMean_ext_corr);
          
	  //long long deltaT_Mean_raw = (0.5*(anEvent->timeL + anEvent->timeR) - timeMean_ext_corr);

          // --- Histograms:
          std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth,energyBinAverage));

          if( h1_deltaT_L_raw[index2] == NULL )
          {
              h1_deltaT_L_raw_raw[index2] = new TH1F(Form("h1_deltaT_L_raw_raw_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
              h1_deltaT_R_raw_raw[index2] = new TH1F(Form("h1_deltaT_R_raw_raw_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
              h1_deltaT_L_raw[index2] = new TH1F(Form("h1_deltaT_L_raw_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
              h1_deltaT_R_raw[index2] = new TH1F(Form("h1_deltaT_R_raw_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	      h1_energyRatioDUT[index2] = new TH1F(Form("h1_energyRatioDUT_%s",labelLR_energyBin.c_str()),"",1000,0.,5.);
            //  h1_deltaT_Mean_raw[index2] = new TH1F(Form("h1_deltaT_Mean_raw_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
          }

          h1_deltaT_L_raw_raw[index2] -> Fill(anEvent->timeL - timeMean_ext_raw);
	  h1_deltaT_R_raw_raw[index2] -> Fill(anEvent->timeR - timeMean_ext_raw);
	  h1_deltaT_L_raw[index2]->Fill(deltaT_L_raw);
          h1_deltaT_R_raw[index2]->Fill(deltaT_R_raw);
          //h1_deltaT_Mean_raw[index2] -> Fill(deltaT_Mean_raw);

          h1_energyRatioDUT[index2] -> Fill( anEvent->energyR / anEvent->energyL );

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

        } // end loop over entries
    }
  std::cout << std::endl;
  std::cout << " [Loop 5] ---> Done! " << std::endl;
  std::cout << std::endl;




  //////////////////////////////
  //--- Draw plots loop 5    //
  /////////////////////////////

  std::map<double,TF1*> fitFunc_deltaT_Raw_L;
  std::map<double,TF1*> fitFunc_deltaT_Raw_R;
  std::map<double,TF1*> fitFunc_deltaT_Raw_Mean;

  std::map<double,float> CTRMeans_L;
  std::map<double,float> CTRSigmas_L;
  std::map<double,float> CTRMeans_R;
  std::map<double,float> CTRSigmas_R;
  ComputeCTRFromHistoMap(h1_deltaT_L_raw, CTRMeans_L, CTRSigmas_L);
  ComputeCTRFromHistoMap(h1_deltaT_R_raw, CTRMeans_R, CTRSigmas_R);

  std::map<std::pair<double,std::string>, double > offsets_deltaT_raw;
  std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> > mean_fitFunc_deltaT_Raw;
  std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> > Sigma_deltaT_Raw;
  std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> > Sigma_deltaT_Raw_Raw;


for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth = map_ths[stepLabel];

      for(int iBar = 0; iBar < 16; ++iBar)
        {
          bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
          if (!barFound) continue;

          std::string labelLR(Form("bar%02d_%s",iBar,stepLabel.c_str()));

          int index1( (10000*int(Vov*100.)) + (100*vth) + iBar );

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
                            TF1* f = DrawAndFitHistogram(
                                               h1_energyRatioDUT[index2],
                                               Form("c_energyRatio_%s", labelLR_energyBin.c_str()),
                                               "energy_{right} / energy_{left}",
                                               plotDir,"Loop5/Optional/EnergyRatioDUT_LR",
                                               iBar,Vov,int(vth));
                            delete f;



              //---------------------------------------------------------------------
              // --- Delta T Raw RAW (DUT - REF)
              //---------------------------------------------------------------------

              // --- Channel L 
              TF1* f_L_raw_raw  = DrawAndFitHistogram(
                                             h1_deltaT_L_raw_raw[index2],
                                             Form("c_deltaT_L_raw_raw_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{L} Raw RAW[ps]",
                                             plotDir,"Loop5/DeltaT_Raw_Raw",
                                             iBar,Vov,int(vth),"L",
                                             kBlack,3.);

              // --- Channel R              
              TF1* f_R_raw_raw = DrawAndFitHistogram(
                                             h1_deltaT_R_raw_raw[index2],
                                             Form("c_deltaT_R_raw_raw_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{R} Raw Raw [ps]",
                                             plotDir,"Loop5/DeltaT_Raw_Raw",
                                             iBar,Vov,int(vth),"R",
                                             kBlack,3.);

	      
              if (f_L_raw_raw && f_R_raw_raw) {

              Sigma_deltaT_Raw_Raw[std::make_tuple(iBar,Vov,vth,"L")]= std::make_pair(f_L_raw_raw->GetParameter(2), f_L_raw_raw->GetParError(2));
              Sigma_deltaT_Raw_Raw[std::make_tuple(iBar,Vov,vth,"R")]= std::make_pair(f_R_raw_raw->GetParameter(2), f_R_raw_raw->GetParError(2));
              }

	      delete f_L_raw_raw;
	      delete f_R_raw_raw;



              //---------------------------------------------------------------------
              // --- Delta T Raw (DUT - REF)
              //---------------------------------------------------------------------

              // --- Channel L 
              fitFunc_deltaT_Raw_L[index2] = DrawAndFitHistogram(
                                             h1_deltaT_L_raw[index2],
                                             Form("c_deltaT_L_raw_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{L} Raw [ps]",
                                             plotDir,"Loop5/DeltaT_Raw",
                                             iBar,Vov,int(vth),"L",
                                             kBlack,3.);

              // --- Channel R              
              fitFunc_deltaT_Raw_R[index2] = DrawAndFitHistogram(
                                             h1_deltaT_R_raw[index2],
                                             Form("c_deltaT_R_raw_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{R} Raw [ps]",
                                             plotDir,"Loop5/DeltaT_Raw",
                                             iBar,Vov,int(vth),"R",
                                             kBlack,3.);


               // --- Fill offsets_deltaT_raw map:
               offsets_deltaT_raw[std::pair(index2,"L")] = fitFunc_deltaT_Raw_L[index2]->GetParameter(1);
               offsets_deltaT_raw[std::pair(index2,"R")] = fitFunc_deltaT_Raw_R[index2]->GetParameter(1);


               // --- Fill the mean Delta T raw map
              if (fitFunc_deltaT_Raw_L[index2] && fitFunc_deltaT_Raw_R[index2]) {
              mean_fitFunc_deltaT_Raw[std::make_tuple(iBar,Vov,vth,"L")] = std::make_pair(fitFunc_deltaT_Raw_L[index2]->GetParameter(1), fitFunc_deltaT_Raw_L[index2]->GetParameter(2));
              mean_fitFunc_deltaT_Raw[std::make_tuple(iBar,Vov,vth,"R")] = std::make_pair(fitFunc_deltaT_Raw_R[index2]->GetParameter(1), fitFunc_deltaT_Raw_R[index2]->GetParameter(2));

              Sigma_deltaT_Raw[std::make_tuple(iBar,Vov,vth,"L")]= std::make_pair(fitFunc_deltaT_Raw_L[index2]->GetParameter(2), fitFunc_deltaT_Raw_L[index2]->GetParError(2));
              Sigma_deltaT_Raw[std::make_tuple(iBar,Vov,vth,"R")]= std::make_pair(fitFunc_deltaT_Raw_R[index2]->GetParameter(2), fitFunc_deltaT_Raw_R[index2]->GetParError(2));
	      
	      }



/*
              //---------------------------------------------------------------------
              // --- Delta T Mean Raw (<DUT> - <REF>)
              //---------------------------------------------------------------------

              fitFunc_deltaT_Raw_Mean[index2] = DrawAndFitHistogram(
                                             h1_deltaT_Mean_raw[index2],
                                             Form("c_deltaT_Mean_raw_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{Mean} Raw [ps]",
                                             plotDir,"Loop5/DeltaT_Mean_Raw",
                                             iBar,Vov,int(vth),"LR",
                                             kBlack,3.);
*/



              //---------------------------------------------------------------------
              // --- Additional plots ToT
              //---------------------------------------------------------------------

               if(optional_plots){


               // -- draw ToT channels side L
               DrawHistogram(h1_ToT_L_TW[index2],
                             Form("c_ToT_L_TW_%s", labelLR_energyBin.c_str()),
                             "ToT [ns]",
                             plotDir,"Loop5/Optional/ToT",
                             Form("ToT_L_TW_%s", labelLR_energyBin.c_str()),
                             BuildLatexSplit(iBar, Vov, vth, h1_ToT_L_TW[index2], "L"),
                             -1.5,1.5,true);

               // -- draw ToT channels side R
               DrawHistogram(h1_ToT_R_TW[index2],
                             Form("c_ToT_R_TW_%s", labelLR_energyBin.c_str()),
                             "ToT [ns]",
                             plotDir,"Loop5/Optional/ToT",
                             Form("ToT_R_TW_%s", labelLR_energyBin.c_str()),
                             BuildLatexSplit(iBar, Vov, vth, h1_ToT_R_TW[index2], "R"),
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
                            plotDir,"Loop5/Optional/DeltaT_vs_ToT",
                            iBar,Vov,vth,"L");

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_vs_ToT_R[index2],
                            nullptr,
                            Form("c_deltaT_vs_ToT_R_%s", labelLR_energyBin.c_str()),
                            "ToT_{R} [ps]","#Delta T [ps]",
                            plotDir,"Loop5/Optional/DeltaT_vs_ToT",
                            iBar,Vov,vth,"R");

               }// optional_plots


            } // --- end loop over energy bins

        } // --- end loop ober bars
   
    } // --- end loop over stepLabels
  std::cout << std::endl;
  std::cout << "[DRAW Loop 5] ---> Done! " << std::endl;
  std::cout << std::endl;

	  

//--------------------------------------------------------------------------
// --- Plot [sigma (Delta T) raw] and [sigma (Delta T) raw raw] vs bar index 
//--------------------------------------------------------------------------
DrawSigmaComparisonVsBar(Sigma_deltaT_Raw_Raw,Sigma_deltaT_Raw,barList,stepLabels,map_Vovs,map_ths,plotDir + "/Loop5/Sigma_DeltaT_Raw_comparison","SigmaComparison_REFCorrection");








  //////////////////////////////
  //--- 6rd loop over events  //
  /////////////////////////////

  for(auto mapIt : trees)
    {
      ModuleEventWithRefClass* anEvent = new ModuleEventWithRefClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);

      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
        {

          if( entry%100000 == 0 ) {
            std::cout << ">>> AMPLITUDE-WALK STUDIES -> Loop 6: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
          }

          mapIt.second -> GetEntry(entry);

          bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
          int index_ref = (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth);
          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );

          long long TW_corr_REF_L = fitFunc_deltaT_EneLCorr_L_ref[index_ref]->Eval(anEvent->energyL_ref);
          long long TW_corr_REF_R = fitFunc_deltaT_EneRCorr_R_ref[index_ref]->Eval(anEvent->energyR_ref);
          long long timeL_ref_corr = anEvent->timeL_ref - TW_corr_REF_L;
          long long timeR_ref_corr = anEvent->timeR_ref - TW_corr_REF_R;
          long long timeMean_ext_corr = 0.5 * (timeL_ref_corr + timeR_ref_corr);
          long long timeMean_ext_raw = 0.5 * (anEvent->timeL_ref + anEvent->timeR_ref);

          auto keyL = std::make_pair(index2,"L");
          auto keyR = std::make_pair(index2,"R");
          long long deltaT_L_raw = (anEvent->timeL - timeMean_ext_corr) - offsets_deltaT_raw[keyL];
          long long deltaT_R_raw = (anEvent->timeR - timeMean_ext_corr) - offsets_deltaT_raw[keyR];

          float energyMeanREF = 0.5 *(anEvent->energyL_ref + anEvent->energyR_ref);

          float deltaT_L_raw_low = - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_L_raw_hig =  5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_R_raw_low = - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_R_raw_hig =  5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);

          //float deltaT_L_raw_low = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          //float deltaT_L_raw_hig = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          //float deltaT_R_raw_low = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          //float deltaT_R_raw_hig = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);


          // --- Histograms:
          std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth,energyBinAverage));

          if( p1_deltaT_L_raw_vs_energyL[index2] == NULL )
            {

              h1_deltaT_L_raw_wOffset[index2] = new TH1F(Form("h1_deltaT_L_raw_wOffset_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
              h1_deltaT_R_raw_wOffset[index2] = new TH1F(Form("h1_deltaT_R_raw_wOffset_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);


              p1_deltaT_L_raw_vs_energyL[index2] = new TProfile(Form("p1_deltaT_L_raw_vs_energyL_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800.);
              //p1_deltaT_L_raw_vs_energyL[index2]->SetErrorOption("s"); //std Dev
              p1_deltaT_R_raw_vs_energyR[index2] = new TProfile(Form("p1_deltaT_R_raw_vs_energyR_%s",labelLR_energyBin.c_str()),"",50,ranges["R"][index1]->at(0),800.);
              //p1_deltaT_R_raw_vs_energyR[index2]->SetErrorOption("s"); //Std Dev

              h2_deltaT_L_raw_vs_energyL[index2] = new TH2F(Form("h2_deltaT_L_raw_vs_energyL_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800., 2000, -12000., 12000.);
              h2_deltaT_R_raw_vs_energyR[index2] = new TH2F(Form("h2_deltaT_R_raw_vs_energyR_%s",labelLR_energyBin.c_str()),"",50,ranges["R"][index1]->at(0),800., 2000, -12000., 12000.);

              if(optional_plots){
              h2_deltaT_L_raw_vs_MeanEnergyREF[index2] = new TH2F(Form("h2_deltaT_L_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),"",50,200.,800., 2000, -12000., 12000.);
              h2_deltaT_R_raw_vs_MeanEnergyREF[index2] = new TH2F(Form("h2_deltaT_R_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),"",50,200.,800., 2000, -12000., 12000.);
              p1_deltaT_L_raw_vs_MeanEnergyREF[index2] = new TProfile(Form("p1_deltaT_L_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
              p1_deltaT_R_raw_vs_MeanEnergyREF[index2] = new TProfile(Form("p1_deltaT_R_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
              }

              h2_EneDUT_L_vs_MeanEneREF[index2] = new TH2F(Form("h2_EneDUT_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",80,ranges["L"][index1]->at(0),800.,80,200.,800.);
              h2_EneDUT_R_vs_MeanEneREF[index2] = new TH2F(Form("h2_EneDUT_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",80,ranges["R"][index1]->at(0),800.,80,200.,800.);

            }


          if((deltaT_L_raw > deltaT_L_raw_low ) && (deltaT_L_raw < deltaT_L_raw_hig))
          {
             
  	     h1_deltaT_L_raw_wOffset[index2]->Fill(deltaT_L_raw);
             p1_deltaT_L_raw_vs_energyL[index2] -> Fill( anEvent->energyL, deltaT_L_raw );
             h2_deltaT_L_raw_vs_energyL[index2] -> Fill( anEvent->energyL, deltaT_L_raw );
             h2_EneDUT_L_vs_MeanEneREF[index2] -> Fill (anEvent->energyL,energyMeanREF);
             
	     if(optional_plots)
	     {
             p1_deltaT_L_raw_vs_MeanEnergyREF[index2] -> Fill( energyMeanREF, deltaT_L_raw);
             h2_deltaT_L_raw_vs_MeanEnergyREF[index2] -> Fill( energyMeanREF, deltaT_L_raw);
             }//optional_plots
          
	  }

          if((deltaT_R_raw > deltaT_R_raw_low ) && (deltaT_R_raw < deltaT_R_raw_hig))
          {
             h1_deltaT_R_raw_wOffset[index2]->Fill(deltaT_R_raw);
             p1_deltaT_R_raw_vs_energyR[index2] -> Fill( anEvent->energyR,deltaT_R_raw );
             h2_deltaT_R_raw_vs_energyR[index2] -> Fill( anEvent->energyR,deltaT_R_raw );
             h2_EneDUT_R_vs_MeanEneREF[index2] -> Fill(anEvent->energyR,energyMeanREF);
             
	     if(optional_plots)
	     {
             p1_deltaT_R_raw_vs_MeanEnergyREF[index2] -> Fill( energyMeanREF, deltaT_R_raw);
             h2_deltaT_R_raw_vs_MeanEnergyREF[index2] -> Fill( energyMeanREF, deltaT_R_raw);
             }//optional_plots
          
	  }


        }
    }


  std::cout << std::endl;
  std::cout << " [Loop 6] ---> Done! " << std::endl;
  std::cout << std::endl;
     	  
                                                   


  //////////////////////////////
  //--- Draw plots loop 6     //
  //////////////////////////////

  std::map<double,TF1*> fitFunc_deltaT_EneLCorr_L;
  std::map<double,TF1*> fitFunc_deltaT_EneRCorr_R;
  std::map<double,TF1*> fitFunc_deltaT_Raw_wOffset_L;
  std::map<double,TF1*> fitFunc_deltaT_Raw_wOffset_R;

  std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> > Chi2_fitFunc_deltaT_vs_Ene;

  double E_low = 400.;
  double E_high = 500;
  std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> > Slope_fitFunc_deltaT_vs_Ene;
  std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> > mean_fitFunc_deltaT_Raw_offset;
  std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> > mean_profile_deltaT_Raw;
  std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> > mean_profileComparison_deltaT_Raw;
  std::map<double,float> CTRMeans_DeltaT_Raw_Offset_L;
  std::map<double,float> CTRSigmas_DeltaT_Raw_Offset_L;
  std::map<double,float> CTRMeans_DeltaT_Raw_Offset_R;
  std::map<double,float> CTRSigmas_DeltaT_Raw_Offset_R;
  ComputeCTRFromHistoMap(h1_deltaT_L_raw_wOffset, CTRMeans_DeltaT_Raw_Offset_L, CTRSigmas_DeltaT_Raw_Offset_L);
  ComputeCTRFromHistoMap(h1_deltaT_R_raw_wOffset, CTRMeans_DeltaT_Raw_Offset_R, CTRSigmas_DeltaT_Raw_Offset_R);

  for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth = map_ths[stepLabel];

      for(int iBar = 0; iBar < 16; ++iBar){


        bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
        if (!barFound) continue;

        std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));

        int index1( (10000*int(Vov*100.)) + (100*vth) + iBar );
        if( !ranges["L-R"][index1] ) continue;

        int nEnergyBins = ranges["L-R"][index1]->size()-1;

        for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
          {
            //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
            double  index2( 10000000*iEnergyBin+index1 );


            if(!p1_deltaT_L_raw_vs_energyL[index2]) continue;

            std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));


              //---------------------------------------------------------------------
              // --- Delta T Raw (DUT - REF) with offset applied
              //---------------------------------------------------------------------

              // --- Channel L 
              fitFunc_deltaT_Raw_wOffset_L[index2] = DrawAndFitHistogram(
                                             h1_deltaT_L_raw_wOffset[index2],
                                             Form("c_deltaT_L_raw_wOffset_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{L} Raw w/ Offset [ps]",
                                             plotDir,"Loop6/DeltaT_Raw_wOffset",
                                             iBar,Vov,int(vth),"L");


              // --- Channel R                    
              fitFunc_deltaT_Raw_wOffset_R[index2] = DrawAndFitHistogram(
                                             h1_deltaT_R_raw_wOffset[index2],
                                             Form("c_deltaT_R_raw_wOffset_%s", labelLR_energyBin.c_str()),
                                             "#Delta T_{R} Raw w/ Offset [ps]",
                                             plotDir,"Loop6/DeltaT_Raw_wOffset",
                                             iBar,Vov,int(vth),"R");


              // --- Fill the mean Delta T raw with offset map
              if (fitFunc_deltaT_Raw_wOffset_L[index2] && fitFunc_deltaT_Raw_wOffset_R[index2]) {
              mean_fitFunc_deltaT_Raw_offset[std::make_tuple(iBar,Vov,vth,"L")] = std::make_pair(fitFunc_deltaT_Raw_wOffset_L[index2]->GetParameter(1), fitFunc_deltaT_Raw_wOffset_L[index2]->GetParameter(2));
              mean_fitFunc_deltaT_Raw_offset[std::make_tuple(iBar,Vov,vth,"R")] = std::make_pair(fitFunc_deltaT_Raw_wOffset_R[index2]->GetParameter(1), fitFunc_deltaT_Raw_wOffset_R[index2]->GetParameter(2));
              }


            //---------------------------------------------------------------------
            // --- Delta T Raw vs Energy DUT
            //---------------------------------------------------------------------

            // --- Channel L:
            fitFunc_deltaT_EneLCorr_L[index2] = DrawAndFitProfile(
                                                p1_deltaT_L_raw_vs_energyL[index2],
                                                CTRMeans_DeltaT_Raw_Offset_L,CTRSigmas_DeltaT_Raw_Offset_L,ranges,false,
                                                index1,index2,
                                                "pol3",
                                                Form("c_deltaT_L_raw_vs_energyL_%s", labelLR_energyBin.c_str()),
                                                "Energy_{L}^{DUT} [a.u.]","#Delta T_{L} raw [ps]",
                                                plotDir,"Loop6/DeltaT_raw_vs_Energy",
                                                iBar, Vov, int(vth),"L",
                                                true, ranges["L"][index1]->at(0),800);


            // --- Channel R:
            fitFunc_deltaT_EneRCorr_R[index2] = DrawAndFitProfile(
                                                p1_deltaT_R_raw_vs_energyR[index2],
                                                CTRMeans_DeltaT_Raw_Offset_R,CTRSigmas_DeltaT_Raw_Offset_R,ranges,false,
                                                index1,index2,
                                                "pol3",
                                                Form("c_deltaT_R_raw_vs_energyR_%s", labelLR_energyBin.c_str()),
                                                "Energy_{R}^{DUT} [a.u.]","#Delta T_{R} raw [ps]",
                                                plotDir,"Loop6/DeltaT_raw_vs_Energy",
                                                iBar, Vov, int(vth),"R",
                                                true, ranges["R"][index1]->at(0),800.);


            if(optional_plots){
            // --- Calculating the aritmetic mean of the points of the TProfile in a certain range 
            auto ProfileMean_L = ComputeProfileMeanInRange(p1_deltaT_L_raw_vs_energyL[index2],320.,400.);
            auto ProfileMean_R = ComputeProfileMeanInRange(p1_deltaT_R_raw_vs_energyR[index2],320.,400.);

            mean_profile_deltaT_Raw[std::make_tuple(iBar,Vov,vth,"L")] = ProfileMean_L;
            mean_profile_deltaT_Raw[std::make_tuple(iBar,Vov,vth,"R")] = ProfileMean_R;


            // --- Calculating the difference between the simple and statistically weighted arithmetic mean for a TProfile:
            auto ProfileMeanComparison_L = ComputeProfileMeansComparison(p1_deltaT_L_raw_vs_energyL[index2],320.,400.);
            auto ProfileMeanComparison_R = ComputeProfileMeansComparison(p1_deltaT_R_raw_vs_energyR[index2],320.,400.);

            mean_profileComparison_deltaT_Raw[std::make_tuple(iBar,Vov,vth,"L")] = std::make_pair(std::get<2>(ProfileMeanComparison_L),1);
            mean_profileComparison_deltaT_Raw[std::make_tuple(iBar,Vov,vth,"R")] = std::make_pair(std::get<2>(ProfileMeanComparison_R),1);


            Chi2_fitFunc_deltaT_vs_Ene[std::make_tuple(iBar,Vov,vth,"L")] = std::make_pair((fitFunc_deltaT_EneLCorr_L[index2]->GetChisquare()/fitFunc_deltaT_EneLCorr_L[index2]->GetNDF()),0.001);
            Chi2_fitFunc_deltaT_vs_Ene[std::make_tuple(iBar,Vov,vth,"R")] = std::make_pair((fitFunc_deltaT_EneRCorr_R[index2]->GetChisquare()/fitFunc_deltaT_EneRCorr_R[index2]->GetNDF()),0.001);

           }// optional_plots



            // --- Fill TProfile fit slope map:
            Slope_fitFunc_deltaT_vs_Ene[std::make_tuple(iBar,Vov,vth,"L")] = std::make_pair( (fitFunc_deltaT_EneLCorr_L[index2]->Eval(E_high) - fitFunc_deltaT_EneLCorr_L[index2]->Eval(E_low)) / (E_high - E_low), 0.001);
            Slope_fitFunc_deltaT_vs_Ene[std::make_tuple(iBar,Vov,vth,"R")] = std::make_pair((fitFunc_deltaT_EneRCorr_R[index2]->Eval(E_high) - fitFunc_deltaT_EneRCorr_R[index2]->Eval(E_low)) / (E_high - E_low), 0.001);



            //-------------------------------------------------------------------------------------
            // --- Delta T Raw vs Eniergy DUT --> TH2F + TProfile
            //-------------------------------------------------------------------------------------

            // --- Channel L:
            DrawTH2WithProfile(
                            h2_deltaT_L_raw_vs_energyL[index2],
                            p1_deltaT_L_raw_vs_energyL[index2],
                            Form("c_deltaT_L_raw_vs_energyL_scatter_%s", labelLR_energyBin.c_str()),
                            "Energy^{DUT}_{L}","#Delta T_{L} [ps]",
                            plotDir,"Loop6/DeltaT_raw_vs_Energy_scatter",
                            iBar,Vov,vth,"L");

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_R_raw_vs_energyR[index2],
                            p1_deltaT_R_raw_vs_energyR[index2],
                            Form("c_deltaT_R_raw_vs_energyR_scatter_%s", labelLR_energyBin.c_str()),
                            "Energy^{DUT}_{R}","#Delta T_{R} [ps]",
                            plotDir,"Loop6/DeltaT_raw_vs_Energy_scatter",
                            iBar,Vov,vth,"R");




           if(optional_plots) {
            //-------------------------------------------------------------------------------------
            // --- Delta T Raw vs Mean Energy REF (reference bar for external time) TH2F + TProfile
            //-------------------------------------------------------------------------------------
 
            // --- Channel L:
            DrawTH2WithProfile(
                            h2_deltaT_L_raw_vs_MeanEnergyREF[index2],
                            p1_deltaT_L_raw_vs_MeanEnergyREF[index2],
                            Form("c_deltaT_raw_vs_MeanEnergyREF_L_%s", labelLR_energyBin.c_str()),
                            "<E_{REF}>","#Delta T_{L} [ps]",
                            plotDir,"Loop6/DeltaT_raw_vs_MeanEnergyREF",
                            iBar,Vov,vth,"L");

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_R_raw_vs_MeanEnergyREF[index2],
                            p1_deltaT_R_raw_vs_MeanEnergyREF[index2],
                            Form("c_deltaT_raw_vs_MeanEnergyREF_R_%s", labelLR_energyBin.c_str()),
                            "<E_{REF}>","#Delta T_{R} [ps]",
                            plotDir,"Loop6/DeltaT_raw_vs_MeanEnergyREF",
                            iBar,Vov,vth,"R");


            }//optional_plots

            //-------------------------------------------------------------------------------------
            // --- Energy DUT vs Mean Energy REF --> TH2F
            //-------------------------------------------------------------------------------------

            // --- Channel L:
            DrawTH2WithProfile(
                            h2_EneDUT_L_vs_MeanEneREF[index2],
                            nullptr,
                            Form("c_h2_EneDUT_L_vs_MeanEneREF_%s", labelLR_energyBin.c_str()),
                            "Energy^{DUT}_{L} [a.u.]","<E^{REF}> [a.u.]",
                            plotDir,"Loop6/EneDUT_vs_MeanEneREF",
                            iBar,Vov,vth,"L");

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_EneDUT_R_vs_MeanEneREF[index2],
                            nullptr,
                            Form("c_h2_EneDUT_R_vs_MeanEneREF_%s", labelLR_energyBin.c_str()),
                            "Energy^{DUT}_{R} [a.u.]","<E^{REF}> [a.u.]",
                            plotDir,"Loop6/EneDUT_vs_MeanEneREF",
                            iBar,Vov,vth,"R");


          }
      }
    }
  std::cout << std::endl;
  std::cout << "[DRAW Loop 6] ---> Done! " << std::endl;
  std::cout << std::endl;




  //----------------------------------------------------------
  //  TGraphErrors: Mean Delta T Raw from gaussian fit vs bar
  //----------------------------------------------------------
  BuildGraphsErrors(mean_fitFunc_deltaT_Raw,g_MeanfitFunc_deltaT_Raw_L,g_MeanfitFunc_deltaT_Raw_R);
  PlotGraphsErrors(g_MeanfitFunc_deltaT_Raw_L,g_MeanfitFunc_deltaT_Raw_R,"Mean_DeltaT_Raw_TGraphError","#mu_{#DeltaT} [ps] ",plotDir + "Loop6/MeanfitFunc_deltaT_Raw", 0., 3000.);


  //----------------------------------------------------------
  //  TGraph: slope TProfile fit vs bar
  //----------------------------------------------------------
  BuildGraphsErrors(Slope_fitFunc_deltaT_vs_Ene,g_Slope_fitFunc_deltaT_vs_Ene_L,g_Slope_fitFunc_deltaT_vs_Ene_R);
  PlotGraphsErrors(g_Slope_fitFunc_deltaT_vs_Ene_L,g_Slope_fitFunc_deltaT_vs_Ene_R,"Slope_fitFunc_deltaT_vs_Ene","Fit function slope [ps/ADC]",plotDir + "Loop6/Slope_fitFunc_deltaT_vs_Ene", -1.5, 0.5);

    SaveProfileComparisonToTxt(Slope_fitFunc_deltaT_vs_Ene,"/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/TProfile_fit_slope/Limited_FitRange");



if(optional_plots){

  //----------------------------------------------------------
  //  TGraphErrors: Chi2/NDF TProfile fit vs bar
  //----------------------------------------------------------
  BuildGraphsErrors(Chi2_fitFunc_deltaT_vs_Ene,g_Chi2_fitFunc_deltaT_vs_Ene_L,g_Chi2_fitFunc_deltaT_vs_Ene_R);
  PlotGraphsErrors(g_Chi2_fitFunc_deltaT_vs_Ene_L,g_Chi2_fitFunc_deltaT_vs_Ene_R,"Chi2_fitFunc_deltaT_vs_Ene","#chi^{2}/NDF",plotDir + "Loop6/Chi2_fitFunc_deltaT_vs_Ene", 0., 15.);

  SaveProfileComparisonToTxt(Chi2_fitFunc_deltaT_vs_Ene,"/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/TProfile_fit_pol2_vs_pol3/chi2/pol3");


  //---------------------------------------------------------------
  //  TGraphErrors: Aritmetic Mean Delta T Raw from TProfile vs bar
  //--------------------------------------------------------------  
  BuildGraphsErrors(mean_profile_deltaT_Raw,g_mean_profile_deltaT_Raw_L,g_mean_profile_deltaT_Raw_R);
  PlotGraphsErrors(g_mean_profile_deltaT_Raw_L,g_mean_profile_deltaT_Raw_R,"Mean_DeltaT_Raw_TProfile","#mu_{#DeltaT} [ps] ",plotDir + "Loop6/MeanProfile_deltaT_Raw", 0., 3000.);


  //-------------------------------------------------------------------------------------------------
  // --- TGraph difference between (weithed mean) - (aritmetic mean) of Delta T raw form TProfile
  //-------------------------------------------------------------------------------------------------
  BuildGraphsErrors(mean_profileComparison_deltaT_Raw,g_mean_profileComparison_deltaT_Raw_L,g_mean_profileComparison_deltaT_Raw_R);
  PlotGraphsErrors(g_mean_profileComparison_deltaT_Raw_L,g_mean_profileComparison_deltaT_Raw_R,"Difference_Mean_DeltaT_Raw_TProfile","#mu_{#DeltaT}^{Weighted}-#mu_{#DeltaT}^{Simple} [ps] ",plotDir + "Loop6/MeanProfileComparison_deltaT_Raw", -100., 100.);


  // --- Save values in txt files:
  SaveProfileComparisonToTxt(
    mean_profileComparison_deltaT_Raw,"/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/Mean_DeltaT_vs_LO/DeltaT_simple_weighted");


  SaveProfileComparisonToTxt(
    mean_profile_deltaT_Raw,"/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/Mean_DeltaT_vs_LO/mean_profile_deltaT_Raw");


  SaveProfileComparisonToTxt(
    mean_fitFunc_deltaT_Raw_offset,"/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/Mean_DeltaT_vs_LO/mean_fitFunc_deltaT_Raw_offset");


  //---------------------------------------------------------------------
  // --- Plot Amplitude walk correction function (T1F) for comparison
  //---------------------------------------------------------------------
  DrawOverlayFits_CorrFunction(fitFunc_deltaT_EneLCorr_L, plotDir + "/Loop6/TW_func_comparison/Side_L","L");
  DrawOverlayFits_CorrFunction(fitFunc_deltaT_EneRCorr_R, plotDir + "/Loop6/TW_func_comparison/Side_R","R");

}//optional_plots


  //-----------------------------------------------------------------------------------
  // --- Plot Amplitude walk correction function (T1F) comparison same bar different vth
  //-----------------------------------------------------------------------------------
  DrawOverlayFits_SameBar_DifferentVth(fitFunc_deltaT_EneLCorr_L,plotDir + "/Loop6/Fit_Function_comparison_differentevth/Side_L","L",-500., 500.);
  DrawOverlayFits_SameBar_DifferentVth(fitFunc_deltaT_EneRCorr_R,plotDir + "/Loop6/Fit_Function_comparison_differentevth/Side_R","R",-500., 500.);

   






  //////////////////////////////
  //--- 7th loop over events  //
  /////////////////////////////

  for(auto mapIt : trees)
    {
      ModuleEventWithRefClass* anEvent = new ModuleEventWithRefClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);

      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
        {
          if( entry%100000 == 0 ){
            std::cout << ">>> TIME-WALK loop 7: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
          }
          mapIt.second -> GetEntry(entry);

          bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
          int index_ref = (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth);
          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
          //double indexBarID( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) );

	  long long TW_corr_REF_L = fitFunc_deltaT_EneLCorr_L_ref[index_ref]->Eval(anEvent->energyL_ref);
          long long TW_corr_REF_R = fitFunc_deltaT_EneRCorr_R_ref[index_ref]->Eval(anEvent->energyR_ref);
          long long timeL_ref_corr = anEvent->timeL_ref - TW_corr_REF_L;
          long long timeR_ref_corr = anEvent->timeR_ref - TW_corr_REF_R;
          long long timeMean_ext_corr = 0.5 * (timeL_ref_corr + timeR_ref_corr);
          long long timeMean_ext_raw = 0.5 * (anEvent->timeL_ref + anEvent->timeR_ref);

          auto keyL = std::make_pair(index2,"L");
          auto keyR = std::make_pair(index2,"R");
          long long deltaT_L_raw = (anEvent->timeL - timeMean_ext_corr) - offsets_deltaT_raw[keyL];
          long long deltaT_R_raw = (anEvent->timeR - timeMean_ext_corr) - offsets_deltaT_raw[keyR];

          //float t1fineMean = 0.5 * ( anEvent->t1fineR + anEvent->t1fineL );
          float energyMeanREF = 0.5 *(anEvent->energyL_ref + anEvent->energyR_ref);

          float deltaT_L_raw_low = - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_L_raw_hig =   5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_R_raw_low = - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_R_raw_hig =   5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          //float deltaT_L_raw_low = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          //float deltaT_L_raw_hig = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          //float deltaT_R_raw_low = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          //float deltaT_R_raw_hig = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);


          // -- TW correction term --> From fitFunc_deltaT_EneLCorr_L/R:
          float TW_corr_L = fitFunc_deltaT_EneLCorr_L[index2]->Eval(anEvent->energyL);
          float TW_corr_R = fitFunc_deltaT_EneRCorr_R[index2]->Eval(anEvent->energyR);


	  // --- Histograms:
          std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth,energyBinAverage));

          if( h1_deltaT_L_TWcorr[index2] == NULL )
          {

            h1_deltaT_L_TWcorr[index2] = new TH1F(Form("h1_deltaT_L_TWcorr_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
            h1_deltaT_R_TWcorr[index2] = new TH1F(Form("h1_deltaT_R_TWcorr_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);

            h2_deltaT_TWCorr_vs_energy_L[index2] = new TH2F(Form("h2_deltaT_TWCorr_vs_energy_L_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800., 2000, -12000., 12000.);
            h2_deltaT_TWCorr_vs_energy_R[index2] = new TH2F(Form("h2_deltaT_TWCorr_vs_energy_R_%s",labelLR_energyBin.c_str()),"",50,ranges["R"][index1]->at(0),800., 2000, -12000., 12000.);

            p1_deltaT_TWCorr_vs_energy_L[index2] = new TProfile(Form("p1_deltaT_TWCorr_vs_energy_L_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800.);    
            p1_deltaT_TWCorr_vs_energy_R[index2] = new TProfile(Form("p1_deltaT_TWCorr_vs_energy_R_%s",labelLR_energyBin.c_str()),"",50,ranges["R"][index1]->at(0),800.);

            p1_deltaT_TWCorrL_vs_t1fineL[index2] = new TProfile(Form("p1_deltaT_TWCorrL_vs_t1fineL_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
            p1_deltaT_TWCorrR_vs_t1fineR[index2] = new TProfile(Form("p1_deltaT_TWCorrR_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0,1000.);

            h2_deltaT_TWCorrL_vs_t1fineL[index2] = new TH2F(Form("h2_deltaT_TWCorrL_vs_t1fineL_%s",labelLR_energyBin.c_str()),"",50,0.,1000., 2000, -12000., 12000.);
            h2_deltaT_TWCorrR_vs_t1fineR[index2] = new TH2F(Form("h2_deltaT_TWCorrR_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0.,1000., 2000, -12000., 12000.);

	    p1_t1fineL_vs_t1fineR[index2] = new TProfile(Form("p1_t1fineL_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0.,1000.);
            h2_t1fineL_vs_t1fineR[index2]= new TH2F(Form("h2_t1fineL_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0.,1000.,50,0.,1000.);



            if(optional_plots)
	    {

            p1_deltaT_TWCorr_L_vs_MeanEneREF[index2] = new TProfile(Form("p1_deltaT_TWCorr_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
            p1_deltaT_TWCorr_R_vs_MeanEneREF[index2] = new TProfile(Form("p1_deltaT_TWCorr_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
            h2_deltaT_TWCorr_L_vs_MeanEneREF[index2] = new TH2F(Form("h2_deltaT_TWCorr_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.,2000, -12000., 12000.);
            h2_deltaT_TWCorr_R_vs_MeanEneREF[index2] = new TH2F(Form("h2_deltaT_TWCorr_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.,2000, -12000., 12000.);
            
	    }//optional_plots

            

	  }


          if((deltaT_L_raw > deltaT_L_raw_low ) && (deltaT_L_raw < deltaT_L_raw_hig)){

          h1_deltaT_L_TWcorr[index2] -> Fill(deltaT_L_raw - TW_corr_L);
          h2_deltaT_TWCorr_vs_energy_L[index2] -> Fill( anEvent->energyL, deltaT_L_raw - TW_corr_L );
          p1_deltaT_TWCorr_vs_energy_L[index2] -> Fill( anEvent->energyL, deltaT_L_raw - TW_corr_L );

          p1_deltaT_TWCorrL_vs_t1fineL[index2] -> Fill(anEvent->t1fineL, deltaT_L_raw - TW_corr_L);
          h2_deltaT_TWCorrL_vs_t1fineL[index2] -> Fill(anEvent->t1fineL, deltaT_L_raw - TW_corr_L);

          p1_t1fineL_vs_t1fineR[index2] -> Fill ( anEvent->t1fineL, anEvent->t1fineR );
          h2_t1fineL_vs_t1fineR[index2] -> Fill(anEvent->t1fineL, anEvent->t1fineR);

          if (optional_plots)
	  {
	  p1_deltaT_TWCorr_L_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_L_raw - TW_corr_L);
          h2_deltaT_TWCorr_L_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_L_raw - TW_corr_L);
	  }//optional_plots

          }


          if((deltaT_R_raw > deltaT_R_raw_low ) && (deltaT_R_raw < deltaT_R_raw_hig)){

          h1_deltaT_R_TWcorr[index2] -> Fill(deltaT_R_raw - TW_corr_R);
          p1_deltaT_TWCorr_vs_energy_R[index2] -> Fill( anEvent->energyR, deltaT_R_raw - TW_corr_R );
          h2_deltaT_TWCorr_vs_energy_R[index2] -> Fill( anEvent->energyR, deltaT_R_raw - TW_corr_R );

          p1_deltaT_TWCorrR_vs_t1fineR[index2] -> Fill(anEvent->t1fineR, deltaT_R_raw - TW_corr_R);
          h2_deltaT_TWCorrR_vs_t1fineR[index2] -> Fill(anEvent->t1fineR, deltaT_R_raw - TW_corr_R);

          if (optional_plots)
	  {
          p1_deltaT_TWCorr_R_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_R_raw - TW_corr_R);
          h2_deltaT_TWCorr_R_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_R_raw - TW_corr_R);
          }//optional_plots

          }

        }
    }

  std::cout << std::endl;
  std::cout << "[Loop 7] ---> Done! " << std::endl;
  std::cout << std::endl;





  //////////////////////////////
  //--- Draw plots loop 7     //
  //////////////////////////////

  //std::map<double,TF1*> fitFunc_deltaT_TWcorr_L;
  //std::map<double,TF1*> fitFunc_deltaT_TWcorr_R;
  std::map<double,TF1*> fitFunc_deltaT_TWCorrL_t1fineL;
  std::map<double,TF1*> fitFunc_deltaT_TWCorrR_t1fineR;

  std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> > Sigma_deltaT_EneCorr;

  std::map<double,float> CTRMeans_TW_L;
  std::map<double,float> CTRSigmas_TW_L;
  std::map<double,float> CTRMeans_TW_R;
  std::map<double,float> CTRSigmas_TW_R;
  ComputeCTRFromHistoMap(h1_deltaT_L_TWcorr, CTRMeans_TW_L, CTRSigmas_TW_L);
  ComputeCTRFromHistoMap(h1_deltaT_R_TWcorr, CTRMeans_TW_R, CTRSigmas_TW_R);


  for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth = map_ths[stepLabel];

      for(int iBar = 0; iBar < 16; ++iBar)
        {
          bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
          if (!barFound) continue;

          std::string labelLR(Form("bar%02d_%s",iBar,stepLabel.c_str()));

          int index1( (10000*int(Vov*100.)) + (100*vth) + iBar );

          if( !ranges["L-R"][index1] ) continue;

          int nEnergyBins = ranges["L-R"][index1]->size()-1;

          for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
            {
              //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
              double index2( 10000000*iEnergyBin+index1 );


              if (!h1_deltaT_L_TWcorr[index2]) continue;
              std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));

              //---------------------------------------------------------------------
              // --- Delta T with TW correction
              //---------------------------------------------------------------------

              // --- Channel L:
              TF1* f_L = DrawAndFitHistogram(
                              h1_deltaT_L_TWcorr[index2],
                              Form("c_h1_deltaT_L_TWcorr_%s", labelLR_energyBin.c_str()),
                              "#Delta T_{L} w/ TW corr [ps]",
                              plotDir,"Loop7/TW_EneDUT/DeltaT_TW_corr",
                              iBar,Vov,int(vth),"L");


              // --- Channel R:
              TF1* f_R = DrawAndFitHistogram(
                              h1_deltaT_R_TWcorr[index2],
                              Form("c_h1_deltaT_R_TWcorr_%s", labelLR_energyBin.c_str()),
                              "#Delta T_{R} w/ TW corr [ps]",
                              plotDir,"Loop7/TW_EneDUT/DeltaT_TW_corr",
                              iBar,Vov,int(vth),"R");


              // --- Fill the mean Delta T raw with offset map

              if (f_L && f_R) {

              Sigma_deltaT_EneCorr[std::make_tuple(iBar,Vov,vth,"L")]= std::make_pair(f_L->GetParameter(2), f_L->GetParError(2));
              Sigma_deltaT_EneCorr[std::make_tuple(iBar,Vov,vth,"R")]= std::make_pair(f_R->GetParameter(2), f_R->GetParError(2));
              }

              delete f_L;
              delete f_R;

            //-------------------------------------------------------------------------------------
            // --- TH2F + TProfile: Delta T with TW correction vs Energy DUT
            //-------------------------------------------------------------------------------------

            // --- Channel L:
            DrawTH2WithProfile(
                            h2_deltaT_TWCorr_vs_energy_L[index2],
                            p1_deltaT_TWCorr_vs_energy_L[index2],
                            Form("c_deltaT_TWCorr_L_vs_energyL_%s", labelLR_energyBin.c_str()),
                            "Energy_{L}^{DUT} [a.u.]","#Delta T_{L} w/ TW corr [ps]",
                            plotDir,"Loop7/TW_EneDUT/DeltaT_TWCorr_vs_Energy",
                            iBar,Vov,vth,"L");

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_TWCorr_vs_energy_R[index2],
                            p1_deltaT_TWCorr_vs_energy_R[index2],
                            Form("c_deltaT_TWCorr_R_vs_energyR_%s", labelLR_energyBin.c_str()),
                            "Energy_{R}^{DUT} [a.u.]","#Delta T_{R} w/ TW corr [ps]",
                            plotDir,"Loop7/TW_EneDUT/DeltaT_TWCorr_vs_Energy",
                            iBar,Vov,vth,"R");




            //---------------------------------------------------------------------
            // --- TProfile Delta T with TW correction vs t1Fine
            //---------------------------------------------------------------------

            // --- Channel L:
            fitFunc_deltaT_TWCorrL_t1fineL[index2] = DrawAndFitProfile(
                                                p1_deltaT_TWCorrL_vs_t1fineL[index2],
                                                CTRMeans_TW_L,CTRSigmas_TW_L,ranges,false,
                                                index1,index2,
                                                "pol3",
                                                Form("c_deltaT_TWCorr_L_vs_t1fineL_%s", labelLR_energyBin.c_str()),
                                                "t1fine_{L}^{DUT}","#Delta T_{L} w/ TW corr [ps]",
                                                plotDir,"Loop7/TW_EneDUT/DeltaT_TWCorr_vs_t1fine",
                                                iBar, Vov, int(vth),"L");


            // --- Channel R:
            fitFunc_deltaT_TWCorrR_t1fineR[index2] = DrawAndFitProfile(
                                                p1_deltaT_TWCorrR_vs_t1fineR[index2],
                                                CTRMeans_TW_R,CTRSigmas_TW_R,ranges,false,
                                                index1,index2,
                                                "pol3",
                                                Form("c_deltaT_TWCorr_R_vs_t1fineR_%s", labelLR_energyBin.c_str()),
                                                "t1fine_{R}^{DUT}","#Delta T_{R} w/ TW corr [ps]",
                                                plotDir,"Loop7/TW_EneDUT/DeltaT_TWCorr_vs_t1fine",
                                                iBar, Vov, int(vth),"R");
            //-------------------------------------------------------------------------------------
            // --- TH2F + TProfile: Delta T with TW correction vs t1Fine DUT
            //-------------------------------------------------------------------------------------

            // --- Channel L:
            DrawTH2WithProfile(
                            h2_deltaT_TWCorrL_vs_t1fineL[index2],
                            p1_deltaT_TWCorrL_vs_t1fineL[index2],
                            Form("c_deltaT_TWCorrL_vs_t1fineL_scatter_%s", labelLR_energyBin.c_str()),
                            "t1fine_{L}^{DUT}","#Delta T_{L} w/ TW corr [ps]",
                            plotDir,"Loop7/TW_EneDUT/DeltaT_TWCorr_vs_t1fine_scatter",
                            iBar,Vov,vth,"L");


            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_TWCorrR_vs_t1fineR[index2],
                            p1_deltaT_TWCorrR_vs_t1fineR[index2],
                            Form("c_deltaT_TWCorrR_vs_t1fineR_scatter_%s", labelLR_energyBin.c_str()),
                            "t1fine_{R}^{DUT}","#Delta T_{R} w/ TW corr [ps]",
                            plotDir,"Loop7/TW_EneDUT/DeltaT_TWCorr_vs_t1fine_scatter",
                            iBar,Vov,vth,"R");


           if(optional_plots){
            //-------------------------------------------------------------------------------------
            // --- TH2F + TProfile: Delta T with TW correction vs vs Mean Energy REF
            //-------------------------------------------------------------------------------------

            // --- Channel L:
            DrawTH2WithProfile(
                            h2_deltaT_TWCorr_L_vs_MeanEneREF[index2],
                            p1_deltaT_TWCorr_L_vs_MeanEneREF[index2],
                            Form("c_deltaT_TWCorr_L_vs_MeanEneREF_%s", labelLR_energyBin.c_str()),
                            "<E^{REF}>","#Delta T_{L} w/ TW corr [ps]",
                            plotDir,"Loop7/TW_EneDUT/Optional/DeltaT_TWCorr_vs_MeanEneREF",
                            iBar,Vov,vth,"L");

            // --- Channel R:
            DrawTH2WithProfile(
                            h2_deltaT_TWCorr_R_vs_MeanEneREF[index2],
                            p1_deltaT_TWCorr_R_vs_MeanEneREF[index2],
                            Form("c_deltaT_TWCorr_R_vs_MeanEneREF_%s", labelLR_energyBin.c_str()),
                            "<E^{REF}>","#Delta T_{R} w/ TW corr [ps]",
                            plotDir,"Loop7/TW_EneDUT/Optional/DeltaT_TWCorr_vs_MeanEneREF",
                            iBar,Vov,vth,"R");


            }// optional_plots

            //-------------------------------------------------------------------------------------
            // --- TH2F + TProfile: f1Fine ch L vs t1Fine ch R (DUT module)
            //-------------------------------------------------------------------------------------

            DrawTH2WithProfile(
                            h2_t1fineL_vs_t1fineR[index2],
                            p1_t1fineL_vs_t1fineR[index2],
                            Form("c_h2_t1fineL_vs_t1fineR_scatter_%s", labelLR_energyBin.c_str()),
                            "t1fine_{L}^{DUT}","t1fine_{R}^{DUT}",
                            plotDir,"Loop7/t1fineL_vs_t1fineR_scatter",
                            iBar,Vov,vth,"L_R");

            }
        }

    }


  std::cout << std::endl;
  std::cout << "[DRAW Loop 7] ---> Done! " << std::endl;
  std::cout << std::endl;


if(optional_plots){

  //----------------------------------------------------------
  //  TGraphErrors: Chi2/NDF TProfile fit vs bar
  //----------------------------------------------------------
  BuildGraphsErrors(Sigma_deltaT_EneCorr,g_Sigma_deltaT_EneCorr_L,g_Sigma_deltaT_EneCorr_R);
  PlotGraphsErrors(g_Sigma_deltaT_EneCorr_L,g_Sigma_deltaT_EneCorr_R,"Sigma_deltaT_EneCorr","#sigma_{#Delta T} [ps]",plotDir + "Loop7/Sigma_deltaT_EneCorr_vs_bar", 60.,85.);

  SaveProfileComparisonToTxt(Sigma_deltaT_EneCorr,"/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/TProfile_fit_pol2_vs_pol3/Sigma/pol2");

}//optional_plots







  //////////////////////////////
  //--- 8th loop over events  //
  /////////////////////////////
  for(auto mapIt : trees)
    {
      
      ModuleEventWithRefClass* anEvent = new ModuleEventWithRefClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);

      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
        {

          if( entry%100000 == 0 ) {
            std::cout << ">>> AMPLITUDE WALK CORRECTION: Loop 8: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
          }

          mapIt.second -> GetEntry(entry);

          bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
          int index_ref = (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth);
          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );

          long long TW_corr_REF_L = fitFunc_deltaT_EneLCorr_L_ref[index_ref]->Eval(anEvent->energyL_ref);
          long long TW_corr_REF_R = fitFunc_deltaT_EneRCorr_R_ref[index_ref]->Eval(anEvent->energyR_ref);
          long long timeL_ref_corr = anEvent->timeL_ref - TW_corr_REF_L;
          long long timeR_ref_corr = anEvent->timeR_ref - TW_corr_REF_R;
          long long timeMean_ext_corr = 0.5 * (timeL_ref_corr + timeR_ref_corr);
          long long timeMean_ext_raw = 0.5 * (anEvent->timeL_ref + anEvent->timeR_ref);

          auto keyL = std::make_pair(index2,"L");
          auto keyR = std::make_pair(index2,"R");
          long long deltaT_L_raw = (anEvent->timeL - timeMean_ext_corr) - offsets_deltaT_raw[keyL];
          long long deltaT_R_raw = (anEvent->timeR - timeMean_ext_corr) - offsets_deltaT_raw[keyR];

	  float energyMeanREF = 0.5 *(anEvent->energyL_ref + anEvent->energyR_ref);

          float deltaT_L_raw_low = - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_L_raw_hig =  5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_R_raw_low = - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_R_raw_hig =  5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          //float deltaT_L_raw_low = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          //float deltaT_L_raw_hig = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          //float deltaT_R_raw_low = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          //float deltaT_R_raw_hig = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);

	  // -- TW corr only fitFunc(energy):
          float TW_corr_L = fitFunc_deltaT_EneLCorr_L[index2]->Eval(anEvent->energyL);
          float TW_corr_R = fitFunc_deltaT_EneRCorr_R[index2]->Eval(anEvent->energyR);
	  
	  // -- TW + t1fine corr
	  //float TW_corr_t1fine_L = fitFunc_deltaT_TWCorrL_t1fineL[index2]->Eval(anEvent->t1fineL);
	  //float TW_corr_t1fine_R = fitFunc_deltaT_TWCorrR_t1fineR[index2]->Eval(anEvent->t1fineR);

	  float TW_corr_t1fine_L = GetProfileBinCorrection(p1_deltaT_TWCorrL_vs_t1fineL[index2],anEvent->t1fineL); //--> corr bin x bin
          float TW_corr_t1fine_R = GetProfileBinCorrection(p1_deltaT_TWCorrR_vs_t1fineR[index2],anEvent->t1fineR);

	  // -- Mean t1fine REF
	  float Mean_t1fineREF = 0.5*(anEvent->t1fineL_ref + anEvent->t1fineR_ref);


	  // --- Histograms:
          std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth,energyBinAverage));

          if( h1_deltaT_TW_t1fine_Corr_L[index2] == NULL )
          {	  
	  // -- TW + t1fine corr
	  h1_deltaT_TW_t1fine_Corr_L[index2] = new TH1F(Form("h1_deltaT_TW_t1fine_Corr_L_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
          h1_deltaT_TW_t1fine_Corr_R[index2] = new TH1F(Form("h1_deltaT_TW_t1fine_Corr_R_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);

	  // -- Mean t1fine REF
	  p1_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2] = new TProfile(Form("p1_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
          h2_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2] = new TH2F(Form("h2_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0.,1000.,2000,-12000.,12000.);
          p1_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2] = new TProfile(Form("p1_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
          h2_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2] = new TH2F(Form("h2_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0.,1000.,2000,-12000.,12000.);
	  }
	


	  // --- Fill Histograms:
	  if((deltaT_L_raw > deltaT_L_raw_low ) && (deltaT_L_raw < deltaT_L_raw_hig)){

	  h1_deltaT_TW_t1fine_Corr_L[index2] -> Fill(deltaT_L_raw - TW_corr_L - TW_corr_t1fine_L);
	  p1_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_L_raw - TW_corr_L - TW_corr_t1fine_L);
	  h2_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_L_raw - TW_corr_L - TW_corr_t1fine_L);
	  }


	  if((deltaT_R_raw > deltaT_R_raw_low ) && (deltaT_R_raw < deltaT_R_raw_hig)){

          h1_deltaT_TW_t1fine_Corr_R[index2] -> Fill(deltaT_R_raw - TW_corr_R -TW_corr_t1fine_R);
  	  p1_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_R_raw - TW_corr_R - TW_corr_t1fine_R);
          h2_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_R_raw - TW_corr_R - TW_corr_t1fine_R);
	  }

        }
    }

  std::cout << std::endl;
  std::cout << "[Loop 8] ---> Done! " << std::endl;
  std::cout << std::endl;







  //////////////////////////////
  //--- Draw plots 8th loop   //
  /////////////////////////////

    //std::map<double,TF1*> fitFunc_deltaT_TW_t1fineCorr_L;
    //std::map<double,TF1*> fitFunc_deltaT_TW_t1fineCorr_R;
    std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> > Sigma_deltaT_Ene_phase_Corr;

    for(auto stepLabel : stepLabels)
    {
        float Vov = map_Vovs[stepLabel];      
        float vth = map_ths[stepLabel];
	    
	for(int iBar = 0; iBar < 16; ++iBar){

            bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;        
            if (!barFound) continue;

	    std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));
        
	    int index1( (10000*int(Vov*100.)) + (100*vth) + iBar );        
	    if( !ranges["L-R"][index1] ) continue;
        
	    int nEnergyBins = ranges["L-R"][index1]->size()-1;

	    for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
	    {
	        //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
		double  index2( 10000000*iEnergyBin+index1 );
            
                std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));


              //---------------------------------------------------------------------
              // --- TH1F Delta T with TW + t1Fine DUT correction
              //---------------------------------------------------------------------

              // --- Channel L:
              TF1* f_L = DrawAndFitHistogram(
    			      h1_deltaT_TW_t1fine_Corr_L[index2],
			      Form("c_h1_deltaT_TW_t1fine_Corr_L_%s", labelLR_energyBin.c_str()),
			      "#Delta T_{L} w/ TW+t1fine^{DUT} corr [ps]",
			      plotDir,"Loop8/TW_EneDUT/DeltaT_TW_t1fine_Corr",
			      iBar,Vov,int(vth),"L");


              // --- Channel R:
              TF1* f_R = DrawAndFitHistogram(
       			      h1_deltaT_TW_t1fine_Corr_R[index2],
			      Form("c_h1_deltaT_TW_t1fine_Corr_R_%s", labelLR_energyBin.c_str()),
			      "#Delta T_{R} w/ TW+t1fine^{DUT} corr [ps]",
			      plotDir,"Loop8/TW_EneDUT/DeltaT_TW_t1fine_Corr",
			      iBar,Vov,int(vth),"R");

              if (f_L && f_R) {

              Sigma_deltaT_Ene_phase_Corr[std::make_tuple(iBar,Vov,vth,"L")]= std::make_pair(f_L->GetParameter(2), f_L->GetParError(2));
              Sigma_deltaT_Ene_phase_Corr[std::make_tuple(iBar,Vov,vth,"R")]= std::make_pair(f_R->GetParameter(2), f_R->GetParError(2));
              }

              delete f_L;
              delete f_R;



            //-------------------------------------------------------------------------------------
            // --- TH2F + TProfile: Delta T with TW + t1fine DUT corr vs Mean t1fine REF
            //-------------------------------------------------------------------------------------

            // --- Channel L:  
            DrawTH2WithProfile(
                            h2_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2],
                            p1_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2],
                            Form("c_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF_%s", labelLR_energyBin.c_str()),
                            "<t1fine^{REF}>","#Delta T_{L} w/ TW+t1fine^{DUT} corr [ps]",
                            plotDir,"Loop8/TW_EneDUT/DeltaT_TW_t1fine_Corr_vs_Meant1fineREF",
                            iBar,Vov,vth,"L");


            // --- Channel R:  
            DrawTH2WithProfile(
                            h2_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2],
                            p1_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2],
                            Form("c_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF_%s", labelLR_energyBin.c_str()),
                            "<t1fine^{REF}>","#Delta T_{R} w/ TW+t1fine^{DUT} corr [ps]",
                            plotDir,"Loop8/TW_EneDUT/DeltaT_TW_t1fine_Corr_vs_Meant1fineREF",
                            iBar,Vov,vth,"R");

	  }
      }
    }



  std::cout << std::endl;
  std::cout << "[Draw Loop 8] ---> Done! " << std::endl;
  std::cout << std::endl;

//------------------------------------------------------------------------------------------------------
// --- Sigma Delta T (Raw, Ene corr, Ene + pahse corr) values vs bar -> individual plot per side L and R
//------------------------------------------------------------------------------------------------------    
DrawSigmaEvolutionVsBar(Sigma_deltaT_Raw,Sigma_deltaT_EneCorr,Sigma_deltaT_Ene_phase_Corr,barList,stepLabels,map_Vovs,map_ths,plotDir + "/SummaryPlots/Sigma_DeltatT_vs_Bar","SigmaEvolution");



//------------------------------------------------------------------------------------------------------
// --- Sigma Delta T (Raw, Ene corr, Ene + pahse corr) values same DUT channels vs vth -> individual plot per side L and R
//------------------------------------------------------------------------------------------------------
DrawSigmaEvolutionVsThreshold(Sigma_deltaT_Raw,Sigma_deltaT_EneCorr,Sigma_deltaT_Ene_phase_Corr,barList,stepLabels,map_Vovs,map_ths,plotDir + "/SummaryPlots/Sigma_DeltatT_vs_vth","SigmaVsThreshold");

/*


  //////////////////////////////
  //--- 9th loop over events  //
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

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );
          int index_ref = (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth);
          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth) + anEvent->barID );

          long long timeMean_ext = 0.5 * (anEvent->timeL_ref + anEvent->timeR_ref);
          //auto keyL = std::make_pair(index2,"L");
          //auto keyR = std::make_pair(index2,"R");
          //long long deltaT_L_raw = (anEvent->timeL - timeMean_ext) - offsets_deltaT_raw[keyL];
          //long long deltaT_R_raw = (anEvent->timeR - timeMean_ext) - offsets_deltaT_raw[keyR];
          long long deltaT_L_raw = (anEvent->timeL - timeMean_ext);
          long long deltaT_R_raw = (anEvent->timeR - timeMean_ext);

          float energyMeanREF = 0.5 *(anEvent->energyL_ref + anEvent->energyR_ref);

          //float deltaT_L_raw_low = - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          //float deltaT_L_raw_hig =  5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          //float deltaT_R_raw_low = - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          //float deltaT_R_raw_hig =  5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_L_raw_low = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_L_raw_hig = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_R_raw_low = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_R_raw_hig = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);

          // -- TW corr only fitFunc(energy):
          float TW_corr_L = fitFunc_deltaT_EneLCorr_L[index2]->Eval(anEvent->energyL);
          float TW_corr_R = fitFunc_deltaT_EneRCorr_R[index2]->Eval(anEvent->energyR);
          
	  // -- TW + t1fine corr
          float TW_corr_t1fine_L = fitFunc_deltaT_TWCorrL_t1fineL[index2]->Eval(anEvent->t1fineL);
          float TW_corr_t1fine_R = fitFunc_deltaT_TWCorrR_t1fineR[index2]->Eval(anEvent->t1fineR);

          // -- Mean t1fine REF
          float Mean_t1fineREF = 0.5*(anEvent->t1fineL_ref + anEvent->t1fineR_ref);
          
	  //float Mean_t1fineREF_corr_L = fitFunc_deltaT_EneRatio_t1fineCorr_L[index2]->Eval(Mean_t1fineREF);
	  //float Mean_t1fineREF_corr_R = fitFunc_deltaT_EneRatio_t1fineCorr_R[index2]->Eval(Mean_t1fineREF);
	  float TW_Mean_t1fineREF_corr_L = GetProfileBinCorrection(p1_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2],Mean_t1fineREF); //--> corr bin x bin
          float TW_Mean_t1fineREF_corr_R = GetProfileBinCorrection(p1_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2],Mean_t1fineREF);


	  // --- Histograms:
          std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth,energyBinAverage));

          
	  if( h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2] == NULL )
          {
	    
	    // -- TW corr + t1fine DUT corr + Mean t1fine REF corr
	    h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2] = new TH1F(Form("h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
            h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2] = new TH1F(Form("h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	  }


          if((deltaT_L_raw > deltaT_L_raw_low ) && (deltaT_L_raw < deltaT_L_raw_hig)){
	  
		  h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2] -> Fill (deltaT_L_raw - TW_corr_L - TW_corr_t1fine_L - TW_Mean_t1fineREF_corr_L);
	  }


	  if((deltaT_R_raw > deltaT_R_raw_low ) && (deltaT_R_raw < deltaT_R_raw_hig)){
	  
		  h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2] -> Fill (deltaT_R_raw - TW_corr_R - TW_corr_t1fine_R - TW_Mean_t1fineREF_corr_R);
	  }

	}
    }

  std::cout << std::endl;
  std::cout << "[Loop 9] ---> Done! " << std::endl;
  std::cout << std::endl;


*/

/*


  //////////////////////////////
  //--- Draw plots 9th loop   //
  /////////////////////////////


    //std::map<double,TF1*> fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L;
    //std::map<double,TF1*> fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R;

    for(auto stepLabel : stepLabels)
    {
            float Vov = map_Vovs[stepLabel];
            float vth = map_ths[stepLabel];

            for(int iBar = 0; iBar < 16; ++iBar){

                    bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
                    if (!barFound) continue;

                    std::string labelLR(Form("bar%02dL-R_%s",iBar,stepLabel.c_str()));

                    int index1( (10000*int(Vov*100.)) + (100*vth) + iBar );
                    if( !ranges["L-R"][index1] ) continue;

                    int nEnergyBins = ranges["L-R"][index1]->size()-1;

                    for(int iEnergyBin = 1; iEnergyBin <= nEnergyBins; ++iEnergyBin)
                    {
                            //if (ranges["L-R"][index1]->at(iEnergyBin)<0) continue;
                            double  index2( 10000000*iEnergyBin+index1 );


                            if(!h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2]) continue;
                            std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));

              //---------------------------------------------------------------------
              // --- TH1F Delta T with TW + t1Fine DUT + Mean t1fine REF correction
              //---------------------------------------------------------------------
              
	      // --- Channel L:
              TF1* f_L = DrawAndFitHistogram(
			      h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2],
			      Form("c_h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L_%s", labelLR_energyBin.c_str()),
			      "#Delta T_{L} w/ TW+t1fine^{DUT}+<t1fine^{REF}> [ps]",
			      plotDir,"Loop6/DeltaT_TW_t1fine_Meant1fineREF_Corr",
			      iBar,Vov,int(vth),"L",
			      kBlue,3.);

	      delete f_L;

              // --- Channel R:
              TF1* f_R= DrawAndFitHistogram(
       			      h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2],
			      Form("c_h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R_%s", labelLR_energyBin.c_str()),
			      "#Delta T_{R} w/ TW+t1fine^{DUT}+<t1fine^{REF}> [ps]",
			      plotDir,"Loop6/DeltaT_TW_t1fine_Meant1fineREF_Corr",
			      iBar,Vov,int(vth),"R",
			      kBlue,3.);
	      delete f_R;

          }
      }
    }
  std::cout << std::endl;
  std::cout << "[Draw Loop 9] ---> Done! " << std::endl;
  std::cout << std::endl;


*/

      int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
