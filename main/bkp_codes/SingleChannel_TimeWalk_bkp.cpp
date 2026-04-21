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
  std::string plotDir = opts.GetOpt<std::string>("Output.plotDir_TW");
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop1_tot_histo/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop1_energy_TW/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop1_MeanTimeREF/",plotDir.c_str()));
 
  // --- From Loop 2:
  system(Form("mkdir -p %s/Loop2_ToT_binned/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_deltaT_Raw/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_energyRatio_LR/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_deltaT_vs_ToT/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_EnergyRatio_REF_TH1/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_deltaT_EneBin/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_DeltaPhase_REF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_DeltaTimeLR_REF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_REF_global/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_NumberEventsBar/",plotDir.c_str()));

  // --- From Loop 3:
  system(Form("mkdir -p %s/Loop3_deltaT_raw_vs_Energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_deltaT_raw_vs_Energy_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_DeltaT_raw_vs_EnergyRatioREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_DeltaT_raw_vs_EnergyRatioREF_scatter/",plotDir.c_str())); 
  system(Form("mkdir -p %s/Loop3_deltaT_raw_vs_MeanEnergyREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_EneDUT_vs_MeanEneREF/",plotDir.c_str()));
  
  // --- From Loop 4:  
  system(Form("mkdir -p %s/Loop4_deltaT_TW_corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_TWCorr_vs_Energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_TWCorr_vs_t1fine/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_TWCorr_vs_t1fine_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_TWCorr_vs_energyRatioREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_TWCorr_vs_energyRatioREF_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_EneRatioREF_corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_EneRatioREFCorr_vs_energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_EneRatioREFCorr_vs_energy_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_EneRatioREFCorr_vs_energyRatioREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_EneRatioREFCorr_vs_t1fine/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_EneRatioREFCorr_vs_t1fine_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Energy_highDeltaT/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_EneRatioREFCorr_vs_MeanEneREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_deltaT_TWCorr_vs_MeanEneREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop4_t1fineL_vs_t1fineR_scatter/",plotDir.c_str()));
  //system(Form("mkdir -p %s/BarID_highDeltaT/",plotDir.c_str()));
  
  //--Loop 5:
  system(Form("mkdir -p %s/Loop5_deltaT_TW_EneRatio_Corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5_deltaT_TW_EneRatio_Corr_vs_Energy/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5_deltaT_TW_EneRatio_Corr_vs_EneRatio/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5_deltaT_TW_EneRatio_Corr_vs_EneRatio_scatter/",plotDir.c_str())); 
 
  system(Form("mkdir -p %s/Loop5_deltaT_EneRatio_TW_Corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5_deltaT_EneRatio_TW_Corr_vs_Energy_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5_deltaT_EneRatio_TW_Corr_vs_EneRatioREF_scatter/",plotDir.c_str())); 

  system(Form("mkdir -p %s/Loop5_deltaT_TW_t1fine_Corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5_deltaT_EneRatioREF_t1fine_Corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5_deltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF_scatter/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5_deltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop5_deltaT_TW_t1fine_Corr_vs_Meant1fineREF/",plotDir.c_str()));

  // -- Loop 6:
  system(Form("mkdir -p %s/Loop6_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop6_deltaT_TW_t1fine_Meant1fineREF_Corr/",plotDir.c_str()));
 
 
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
  std::map<double,TH1F*> h1_energyRatio;
  std::map<double,TH2F*> h2_deltaT_vs_ToT_L;
  std::map<double,TH2F*> h2_deltaT_vs_ToT_R;
  std::map<double,TH1F*> h1_DeltaPhase_REF;
  std::map<double,TH1F*> h1_DeltaTimeLR_REF;
  TH1F* h1_DeltaTimeLR_REF_global;
  TH1F* h1_DeltaPhaseLR_REF_global;
  TH1F* h1_NumberEventsBar;

  //Loop 3:
  std::map<double,TProfile*> p1_deltaT_R_raw_vs_energyR;
  std::map<double,TProfile*> p1_deltaT_L_raw_vs_energyL;
  std::map<double,TH2F*> h2_deltaT_L_raw_vs_energyL;
  std::map<double,TH2F*> h2_deltaT_R_raw_vs_energyR;
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

  //Loop 4:
  std::map<double,TH1F*> h1_deltaT_L_TWcorr;
  std::map<double,TH1F*> h1_deltaT_R_TWcorr;
  std::map<double,TProfile*> p1_deltaT_TWCorr_vs_energy_L;
  std::map<double,TProfile*> p1_deltaT_TWCorr_vs_energy_R;
  std::map<double,TH2F*> h2_deltaT_TWCorr_vs_energy_L;
  std::map<double,TH2F*> h2_deltaT_TWCorr_vs_energy_R;
  std::map<double,TProfile*> p1_deltaT_TWCorr_L_vs_energyRatioREF;
  std::map<double,TProfile*> p1_deltaT_TWCorr_R_vs_energyRatioREF;
  std::map<double,TH2F*> h2_deltaT_TWCorr_L_vs_energyRatioREF;
  std::map<double,TH2F*> h2_deltaT_TWCorr_R_vs_energyRatioREF;
  std::map<double,TProfile*> p1_deltaT_TWCorrL_vs_t1fineL;
  std::map<double,TProfile*> p1_deltaT_TWCorrR_vs_t1fineR;
  std::map<double,TH2F*> h2_deltaT_TWCorrL_vs_t1fineL;
  std::map<double,TH2F*> h2_deltaT_TWCorrR_vs_t1fineR;
  std::map<double,TH1F*> h1_deltaT_L_EneRatioREFcorr;
  std::map<double,TH1F*> h1_deltaT_R_EneRatioREFcorr;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_L_vs_energy;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_R_vs_energy;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_L_vs_energy;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_R_vs_energy;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_L_vs_energyRatioREF;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_R_vs_energyRatioREF;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_L_vs_energyRatioREF;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_R_vs_energyRatioREF;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_L_vs_t1fineL;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_R_vs_t1fineR;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_L_vs_t1fineL;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_R_vs_t1fineR;
  
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_L_vs_MeanEneREF;
  std::map<double,TProfile*> p1_deltaT_EneRatioREFCorr_R_vs_MeanEneREF;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_L_vs_MeanEneREF;
  std::map<double,TH2F*> h2_deltaT_EneRatioREFCorr_R_vs_MeanEneREF;
  std::map<double,TProfile*> p1_deltaT_TWCorr_L_vs_MeanEneREF;
  std::map<double,TProfile*> p1_deltaT_TWCorr_R_vs_MeanEneREF;
  std::map<double,TH2F*> h2_deltaT_TWCorr_L_vs_MeanEneREF;
  std::map<double,TH2F*> h2_deltaT_TWCorr_R_vs_MeanEneREF;

  std::map<double,TH2F*> h2_t1fineL_vs_t1fineR; 
  std::map<double,TProfile*> p1_t1fineL_vs_t1fineR;
 
  // -- Test high deltaT:
  std::map<double,TH1F*> h1_energyL_high_deltaT;
  std::map<double,TH1F*> h1_energyR_high_deltaT;
  std::map<double,TH2F*> h2_energy_vs_deltaT_high_L;
  std::map<double,TH2F*> h2_energy_vs_deltaT_high_R;
  std::map<double,TH1F*> h1_bar_high_deltaT;

  //Loop 5:
  std::map<double,TH1F*> h1_deltaT_TWCorr_EneRatioCorr_L;
  std::map<double,TH2F*> h2_deltaT_TWCorr_EneRatioCorr_vs_energy_L;
  std::map<double,TProfile*> p1_deltaT_TWCorr_EneRatioCorr_vs_energy_L;
  std::map<double,TProfile*> p1_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L;
  std::map<double,TH2F*> h2_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L;
  std::map<double,TH2F*> h2_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R;
  std::map<double,TH1F*> h1_deltaT_TWCorr_EneRatioCorr_R;
  std::map<double,TH2F*> h2_deltaT_TWCorr_EneRatioCorr_vs_energy_R;
  std::map<double,TProfile*> p1_deltaT_TWCorr_EneRatioCorr_vs_energy_R;
  std::map<double,TProfile*> p1_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R;

  std::map<double,TH1F*> h1_deltaT_EneRatioREF_TW_Corr_L;
  std::map<double,TH1F*> h1_deltaT_EneRatioREF_TW_Corr_R;
  std::map<double,TH2F*> h2_deltaT_EneRatioREF_TW_Corr_vs_energy_L;
  std::map<double,TH2F*> h2_deltaT_EneRatioREF_TW_Corr_vs_energy_R;
  std::map<double,TProfile*> p1_deltaT_EneRatioREF_TW_Corr_vs_energy_L;
  std::map<double,TProfile*> p1_deltaT_EneRatioREF_TW_Corr_vs_energy_R;
  std::map<double,TH2F*> h2_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_L;
  std::map<double,TH2F*> h2_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_R;
  std::map<double,TProfile*> p1_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_L;
  std::map<double,TProfile*> p1_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_R;
  
  std::map<double,TH1F*> h1_deltaT_TW_t1fine_Corr_L;
  std::map<double,TH1F*> h1_deltaT_TW_t1fine_Corr_R;
  std::map<double,TH1F*> h1_deltaT_EneRatioREF_t1fine_Corr_L;
  std::map<double,TH1F*> h1_deltaT_EneRatioREF_t1fine_Corr_R;
  std::map<double,TProfile*> p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF;
  std::map<double,TProfile*> p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF;
  std::map<double,TH2F*> h2_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF;
  std::map<double,TH2F*> h2_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF;
  std::map<double,TProfile*> p1_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF;
  std::map<double,TProfile*> p1_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF;
  std::map<double,TH2F*> h2_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF;
  std::map<double,TH2F*> h2_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF;


  //Loop 6:
  std::map<double,TH1F*> h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L;
  std::map<double,TH1F*> h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R;
  std::map<double,TH1F*> h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L;
  std::map<double,TH1F*> h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R;


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

	  //--SELECTION SLICE MEAN ENERGY REF BAR 7:
	  //if((0.5*(anEvent->energyL_ext + anEvent->energyR_ext)) > 320 || (0.5*(anEvent->energyL_ext + anEvent->energyR_ext)) < 240 ) continue;

	  int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

	  if( energyBinAverage < 1 ) continue;

	  accept[index1][entry] = true;

	  double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

	  long long timeMean_ext = 0.5 * (anEvent->timeL_ext + anEvent->timeR_ext);
          long long deltaT_L_raw = anEvent->timeL - timeMean_ext;
          long long deltaT_R_raw = anEvent->timeR - timeMean_ext;

	  float energyMeanREF = 0.5 *(anEvent->energyL_ext + anEvent->energyR_ext);

	  //REF bar only:
	  long long phaseBarREF_L = std::fmod(anEvent->timeL_ext - (6250 * anEvent->t1coarseL_ext), 6250.0); 
	  long long phaseBarREF_R = std::fmod(anEvent->timeR_ext - (6250 * anEvent->t1coarseR_ext), 6250.0);
          long long DeltaT_REF = anEvent->timeL_ext - anEvent->timeR_ext;

//-------------
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
//-------------

	  if( h1_deltaT_L_raw[index2] == NULL )
	    {
	      std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));

	      h1_deltaT_L_raw[index2] = new TH1F(Form("h1_deltaT_L_raw_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
              h1_deltaT_R_raw[index2] = new TH1F(Form("h1_deltaT_R_raw_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	      
	      h1_ToT_L_TW[index2] = new TH1F(Form("h1_ToT_L_TW_%s",labelLR_energyBin.c_str()),"",2000,-60,60.);
              h1_ToT_R_TW[index2] = new TH1F(Form("h1_ToT_R_TW_%s",labelLR_energyBin.c_str()),"",2000,-60,60.);
              h2_deltaT_vs_ToT_L[index2] = new TH2F(Form("h2_deltaT_vs_ToT_L_%s",labelLR_energyBin.c_str()),"",2000,-10, 20, 100, -12000., 12000.);
              h2_deltaT_vs_ToT_R[index2] = new TH2F(Form("h2_deltaT_vs_ToT_R_%s",labelLR_energyBin.c_str()),"",2000,-10, 20, 100, -12000., 12000.);

	    }

	  // -- raw Delta T histograms
	  h1_deltaT_L_raw[index2]->Fill(deltaT_L_raw);      
	  h1_deltaT_R_raw[index2]->Fill(deltaT_R_raw);
	    
          // -- TbT histograms
	  h1_ToT_L_TW[index2] -> Fill(anEvent->totL);
	  h1_ToT_L_TW[index2] -> Fill(anEvent->totR);
	  h2_deltaT_vs_ToT_L[index2] -> Fill( anEvent->totL,deltaT_L_raw);
	  h2_deltaT_vs_ToT_R[index2] -> Fill( anEvent->totR,deltaT_R_raw);


	  // -- energy Ratio histograms
	  if( h1_energyRatio[index2] == NULL )
	    {
	      std::string labelLR_energyBin(Form("bar%02dL-R_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
	      h1_energyRatio[index2] = new TH1F(Form("h1_energyRatio_%s",labelLR_energyBin.c_str()),"",1000,0.,5.);
              h1_energyRatioL_REF[index2] = new TH1F(Form("h1_energyRatioL_REF_%s",labelLR_energyBin.c_str()),"",500,0.,5.);
	      h1_energyRatioR_REF[index2] = new TH1F(Form("h1_energyRatioR_REF_%s",labelLR_energyBin.c_str()),"",500,0.,5.);
	    }
	 
	  h1_energyRatio[index2] -> Fill( anEvent->energyR / anEvent->energyL );
	  h1_energyRatioL_REF[index2] -> Fill( anEvent->energyL/energyMeanREF );
	  h1_energyRatioR_REF[index2] -> Fill( anEvent->energyR/energyMeanREF );


	  // -- Delta Phase REF bar histograms
          if( h1_DeltaPhase_REF[indexREF] == NULL )
            {
              std::string label_ext(Form("Vov%.2f_th%02d",anEvent->Vov,anEvent->vth1));
	      std::string label_ext_global(Form("Vov"));
              h1_DeltaPhase_REF[indexREF] = new TH1F(Form("h1_DeltaPhase_REF_%s",label_ext.c_str()),"",2000,-12000,12000.);
	      h1_DeltaTimeLR_REF[indexREF] = new TH1F(Form("h1_DeltaTimeLR_REF_%s",label_ext.c_str()),"",2000,-12000,12000.);
            }
          h1_DeltaPhase_REF[indexREF] -> Fill(phaseBarREF_L  - phaseBarREF_R);
          h1_DeltaTimeLR_REF[indexREF] ->Fill(DeltaT_REF);



	  // GLOBAL histogram (no indexREF dependence)
          if( h1_DeltaTimeLR_REF_global == NULL )
          {
		  h1_DeltaTimeLR_REF_global = new TH1F("h1_DeltaTimeLR_REF_global",";#DeltaT_{REF} = t_{L} - t_{R} [ps];Entries",2000,-12000,12000.);
		  h1_DeltaPhaseLR_REF_global = new TH1F("h1_DeltaPhaseLR_REF_global",";#Delta #phi_{REF} = #phi_{L} - #phi_{R} [ps];Entries",2000,-12000,12000.);
	  }
	  h1_DeltaTimeLR_REF_global->Fill(DeltaT_REF);
	  if(phaseBarREF_L!=0 && phaseBarREF_R!=0) {
	  h1_DeltaPhaseLR_REF_global->Fill(phaseBarREF_L  - phaseBarREF_R);

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

  std::map<double,float> CTRMeans_L;
  std::map<double,float> CTRSigmas_L;
  std::map<double,TF1*> fitFunc_deltaT_Raw_L;
  std::map<double,TF1*> fitFunc_deltaT_Raw_R;
  std::map<double,TF1*> fitFunc_energyRatio;
  TF1* fitFunc_DeltaPhaseLR_REF_global;
  TF1* fitFunc_DeltaTimeLR_REF_global;



  std::string label_1(Form("global"));


  if(h1_NumberEventsBar!=0) {
   c = new TCanvas(Form("c_h1_NumberEventsBar_%s",label_1.c_str()),Form("c_h1_NumberEventsBar_%s",label_1.c_str()));
   histo = h1_NumberEventsBar;
   //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
   //histo -> GetXaxis() -> SetRangeUser(-1000,1000);
   //histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
   histo -> SetTitle(Form(";DUT bar number;entries"));
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


   c -> Print(Form("%s/Loop2_NumberEventsBar/c_NumberEventsBar_%s.png",plotDir.c_str(),label_1.c_str()));
   c -> Print(Form("%s/Loop2_NumberEventsBar/c_NumberEventsBar_%s.pdf",plotDir.c_str(),label_1.c_str()));
   delete latex;
   delete c;
  }






  for(auto mapIt : h1_deltaT_L_raw)
    {
      double index = mapIt.first;

      FindSmallestInterval(vals,h1_deltaT_L_raw[index],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans_L[index] = mean;
      CTRSigmas_L[index] = effSigma;

    }


  std::map<double,float> CTRMeans_R;
  std::map<double,float> CTRSigmas_R;

  for(auto mapIt : h1_deltaT_R_raw)
    {
      double index = mapIt.first;

      FindSmallestInterval(vals,h1_deltaT_R_raw[index],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans_R[index] = mean;
      CTRSigmas_R[index] = effSigma;
    }

//----------------------------

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
		    //<< " | Xcenter=" << std::setw(8) << x_center
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

//-----------------------------
              std::string label(Form("global"));

	      //Delta Time ChL-ChR REF bar -- Global
              c = new TCanvas(Form("c_DeltaTimeLR_REF_global_%s",label.c_str()),Form("c_DeltaTimeLR_REF_global_%s",label.c_str()));
              histo = h1_DeltaTimeLR_REF_global;
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> GetXaxis() -> SetRangeUser(-1000,1000);
	      histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#Delta t = (t_{R} - t _{L});entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

	      fitFunc_DeltaTimeLR_REF_global = new TF1(Form("fitFunc_DeltaTimeLR_REF_global_%s",label.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
              histo -> Fit(fitFunc_DeltaTimeLR_REF_global,"QNRS");
              histo -> Fit(fitFunc_DeltaTimeLR_REF_global,"QSR+","",fitFunc_DeltaTimeLR_REF_global->GetParameter(1)-2.*fitFunc_DeltaTimeLR_REF_global->GetParameter(2),fitFunc_DeltaTimeLR_REF_global->GetParameter(1)+2.*fitFunc_DeltaTimeLR_REF_global->GetParameter(2));
              histo -> Fit(fitFunc_DeltaTimeLR_REF_global,"QSR+","",fitFunc_DeltaTimeLR_REF_global->GetParameter(1)-2.*fitFunc_DeltaTimeLR_REF_global->GetParameter(2),fitFunc_DeltaTimeLR_REF_global->GetParameter(1)+2.*fitFunc_DeltaTimeLR_REF_global->GetParameter(2));

              fitFunc_DeltaTimeLR_REF_global -> SetLineColor(kBlue);
              fitFunc_DeltaTimeLR_REF_global -> SetLineWidth(2);
              fitFunc_DeltaTimeLR_REF_global -> Draw("same");

	      latex = new TLatex(0.20,0.85,Form("#splitline{REF Module - Bar 8}""{Mean = %.3f, stdDev = %.3f}", histo->GetMean(), histo->GetRMS()));
              latex -> SetNDC();
	      latex->SetTextAlign(12);
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop2_REF_global/c_DeltaTimeLR_REF_global_%s.png",plotDir.c_str(),label.c_str()));
              c -> Print(Form("%s/Loop2_REF_global/c_DeltaTimeLR_REF_global_%s.pdf",plotDir.c_str(),label.c_str()));
              delete latex;
              delete c;



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

//---------------------------------------------------------




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
	      if (!h1_energyRatio[index2]) continue;

	      std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));

	      c = new TCanvas(Form("c_energyRatio_%s",labelLR_energyBin.c_str()),Form("c_energyRatio_%s",labelLR_energyBin.c_str()));
	      histo = h1_energyRatio[index2];
	      histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
	      histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
	      histo -> SetTitle(Form(";energy_{right} / energy_{left};entries"));
	      histo -> SetLineColor(kRed);
	      histo -> SetLineWidth(2);
	      histo -> Draw();
	      histo -> Write();

	      fitFunc_energyRatio[index2] = new TF1(Form("fitFunc_energyRatio_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
	      histo -> Fit(fitFunc_energyRatio[index2],"QNRS");
	      histo -> Fit(fitFunc_energyRatio[index2],"QSR+","",fitFunc_energyRatio[index2]->GetParameter(1)-2.*fitFunc_energyRatio[index2]->GetParameter(2),fitFunc_energyRatio[index2]->GetParameter(1)+2.*fitFunc_energyRatio[index2]->GetParameter(2));
	      histo -> Fit(fitFunc_energyRatio[index2],"QSR+","",fitFunc_energyRatio[index2]->GetParameter(1)-2.*fitFunc_energyRatio[index2]->GetParameter(2),fitFunc_energyRatio[index2]->GetParameter(1)+2.*fitFunc_energyRatio[index2]->GetParameter(2));

	      fitFunc_energyRatio[index2] -> SetLineColor(kBlack);
	      fitFunc_energyRatio[index2] -> SetLineWidth(2);
	      fitFunc_energyRatio[index2] -> Draw("same");

	      //FIXME
	      //fitFunc_energyRatio[index2] -> SetParameter(1,histo->GetMean());
	      //fitFunc_energyRatio[index2] -> SetParameter(2,histo->GetRMS());

	      latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kRed);
	      latex -> Draw("same");

	      c -> Print(Form("%s/Loop2_energyRatio_LR/c_energyRatio__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	      c -> Print(Form("%s/Loop2_energyRatio_LR/c_energyRatio__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	      delete latex;
	      delete c;




	      // -- raw Delta T - Ch L 
              if (!h1_deltaT_L_raw[index2]) continue;

	       c = new TCanvas(Form("c_deltaT_L_raw_%s",labelLR_energyBin.c_str()),Form("c_deltaT_L_raw_%s",labelLR_energyBin.c_str()));
               histo = h1_deltaT_L_raw[index2];
               histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
               histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
               histo -> SetTitle(Form(";#Delta T_{L} [ps];entries"));
               histo -> SetLineColor(kRed);
               histo -> SetLineWidth(2);
	       //c -> SetLogy();
               histo -> Draw();
	       histo -> Write();

	       fitFunc_deltaT_Raw_L[index2] = new TF1(Form("fitFunc_deltaT_Raw_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
               histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QNRS");
               histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QSR+","",fitFunc_deltaT_Raw_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2),fitFunc_deltaT_Raw_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2));
               histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QSR+","",fitFunc_deltaT_Raw_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2),fitFunc_deltaT_Raw_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2));

               fitFunc_deltaT_Raw_L[index2] -> SetLineColor(kBlack);
               fitFunc_deltaT_Raw_L[index2] -> SetLineWidth(2);
               fitFunc_deltaT_Raw_L[index2] -> Draw("same");

	       latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), histo->GetMean(), histo->GetRMS()));
               latex -> SetNDC();
               latex -> SetTextFont(42);
               latex -> SetTextSize(0.04);
               latex -> SetTextColor(kRed);
               latex -> Draw("same");

               c -> Print(Form("%s/Loop2_deltaT_Raw/deltaT_L_raw__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
               c -> Print(Form("%s/Loop2_deltaT_Raw/deltaT_L_raw__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
               delete c;
               delete latex;


	       // -- raw Delta T - Ch R
               if (!h1_deltaT_R_raw[index2]) continue;
	                
	       c = new TCanvas(Form("c_deltaT_R_ext_%s",labelLR_energyBin.c_str()),Form("c_deltaT_R_ext_%s",labelLR_energyBin.c_str()));
               histo = h1_deltaT_R_raw[index2];
               histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
               histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
               histo -> SetTitle(Form(";#Delta T_{R} [ps];entries"));
               histo -> SetLineColor(kRed);
               histo -> SetLineWidth(2);
               histo -> Draw();
               histo -> Write();

	       fitFunc_deltaT_Raw_R[index2] = new TF1(Form("fitFunc_deltaT_Raw_R_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
               histo -> Fit(fitFunc_deltaT_Raw_R[index2],"QNRS");
               histo -> Fit(fitFunc_deltaT_Raw_R[index2],"QSR+","",fitFunc_deltaT_Raw_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_Raw_R[index2]->GetParameter(2),fitFunc_deltaT_Raw_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_Raw_R[index2]->GetParameter(2));
               histo -> Fit(fitFunc_deltaT_Raw_R[index2],"QSR+","",fitFunc_deltaT_Raw_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_Raw_R[index2]->GetParameter(2),fitFunc_deltaT_Raw_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_Raw_R[index2]->GetParameter(2));

               fitFunc_deltaT_Raw_R[index2] -> SetLineColor(kBlack);
               fitFunc_deltaT_Raw_R[index2] -> SetLineWidth(2);
               fitFunc_deltaT_Raw_R[index2] -> Draw("same");

	       latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), histo->GetMean(), histo->GetRMS()));
               latex -> SetNDC();
               latex -> SetTextFont(42);
               latex -> SetTextSize(0.04);
               latex -> SetTextColor(kRed);
               latex -> Draw("same");

               c -> Print(Form("%s/Loop2_deltaT_Raw/deltaT_R_ext__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
               c -> Print(Form("%s/Loop2_deltaT_Raw/deltaT_R_ext__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
               delete c;
               delete latex;

/*
      std::cout << std::fixed << std::setprecision(3);
            std::cout
                    << " | Bar="  << std::setw(2)  << iBar
                    << " | vth=" << std::setw(5) << map_ths[stepLabel]
                    << " | Means=" << std::setw(8) << CTRMeans_TW_L[index2]
                    << " | Sigmas="    << std::setw(8) << CTRSigmas_TW_L[index2]
                    << std::endl;

*/


	       // -- draw ToT Ch L
               if (!h1_ToT_L_TW[index2]) continue;

	       c = new TCanvas(Form("c_ToT_L_TW_%s",labelLR_energyBin.c_str()),Form("c_ToT_L_TW_%s",labelLR_energyBin.c_str()));
	       histo = h1_ToT_L_TW[index2];
	       //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
	       histo -> GetXaxis() -> SetRangeUser(-1.5,1.5);
	       histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
	       histo -> SetTitle(Form(";(ToT)[ns];entries"));
               histo -> SetLineColor(kRed);
	       histo -> SetLineWidth(2);
	       c -> SetLogy();
	       histo -> Draw();
               histo -> Write();

	       latex = new TLatex(0.25,0.85,Form("#splitline{bar %02d}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), histo->GetMean(), histo->GetRMS()));
	       latex -> SetNDC();
	       latex -> SetTextFont(42);
	       latex -> SetTextSize(0.03);
	       latex -> SetTextColor(kRed);
	       latex -> Draw("same");

	       c -> Print(Form("%s/Loop2_ToT_binned/ToT_L_TW_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	       c -> Print(Form("%s/Loop2_ToT_binned/ToT_L_TW_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	       delete c;
	       delete latex;



	       // -- draw ToT Ch R
	       if (!h1_ToT_R_TW[index2]) continue;
	       c = new TCanvas(Form("c_ToT_R_TW_%s",labelLR_energyBin.c_str()),Form("c_ToT_R_TW_%s",labelLR_energyBin.c_str()));
               histo = h1_ToT_R_TW[index2];
               histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
               histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
               histo -> SetTitle(Form(";(ToT)[ns];entries"));
               histo -> SetLineColor(kRed);
               histo -> SetLineWidth(2);
               c -> SetLogy();
               histo -> Draw();
               histo -> Write();

               latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), histo->GetMean(), histo->GetRMS()));
               latex -> SetNDC();
               latex -> SetTextFont(42);
               latex -> SetTextSize(0.04);
               latex -> SetTextColor(kRed);
               latex -> Draw("same");

               c -> Print(Form("%s/Loop2_ToT_binned/ToT_R_TW_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
               c -> Print(Form("%s/Loop2_ToT_binned/ToT_R_TW_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
               delete c;
               delete latex;



	       // -- TH2F ToT vs delta T:
	       if(!h2_deltaT_vs_ToT_L[index2]) continue;

               c = new TCanvas(Form("c_deltaT_vs_ToT_L_%s",labelLR_energyBin.c_str()),Form("c_deltaT_vs_ToT_L_%s",labelLR_energyBin.c_str()));
               c -> SetGridy();

               h2 = h2_deltaT_vs_ToT_L[index2];
               //h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) -600., h2 -> GetMean(2)+600);
               h2 -> SetTitle(Form(";ToT Left;#Deltat Left [ps]"));
               h2 -> Draw("colz");

	       latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
               latex -> SetNDC();
               latex -> SetTextFont(42);
               latex -> SetTextSize(0.04);
               latex -> SetTextColor(kRed);
               latex -> Draw("same");

	       c -> Print(Form("%s/Loop2_deltaT_vs_ToT/c_deltaT_vs_ToT_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	       c -> Print(Form("%s/Loop2_deltaT_vs_ToT/c_deltaT_vs_ToT_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
               delete c;
               delete latex;



	    // -- Energy Ratio REF: [Ene DUT / Mean Ene REF (bar 7)]  
	    
	    //Ch L:
	    if(!h1_energyRatioL_REF[index2]) continue;  
	    
	    h1_energyRatioL_REF[index2];
	    c = new TCanvas(Form("c_energyRatioL_REF_%s",labelLR_energyBin.c_str()),Form("c_energyRatioL_REF_%s",labelLR_energyBin.c_str()));
            histo = h1_energyRatioL_REF[index2];
            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
            histo -> SetTitle(Form(";(E_{L}/E_{REF});entries"));
            histo -> SetLineColor(kRed);
            histo -> SetLineWidth(2);
            //c -> SetLogy();
            histo -> Draw();
            histo -> Write();

	    /*
            fitFunc_deltaT_Raw_L[index2] = new TF1(Form("fitFunc_deltaT_Raw_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
            histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QNRS");
            histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QSR+","",fitFunc_deltaT_Raw_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2),fitFunc_deltaT_Raw_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2));
            histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QSR+","",fitFunc_deltaT_Raw_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2),fitFunc_deltaT_Raw_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2));

            fitFunc_deltaT_Raw_L[index2] -> SetLineColor(kBlack);
            fitFunc_deltaT_Raw_L[index2] -> SetLineWidth(2);
            fitFunc_deltaT_Raw_L[index2] -> Draw("same");

	    */
            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), histo->GetMean(), histo->GetRMS()));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

            c -> Print(Form("%s/Loop2_EnergyRatio_REF_TH1/energyRatioL_REF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop2_EnergyRatio_REF_TH1/energyRatioL_REF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete c;
            delete latex;



	    //Ch R:
	    if(!h1_energyRatioR_REF[index2]) continue;
	    h1_energyRatioR_REF[index2];
            c = new TCanvas(Form("c_energyRatioR_REF_%s",labelLR_energyBin.c_str()),Form("c_energyRatioR_REF_%s",labelLR_energyBin.c_str()));
            histo = h1_energyRatioR_REF[index2];
            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
            histo -> SetTitle(Form(";(E_{R}/E_{REF});entries"));
            histo -> SetLineColor(kRed);
            histo -> SetLineWidth(2);
            //c -> SetLogy();
            histo -> Draw();
            histo -> Write();

            /*
            fitFunc_deltaT_Raw_L[index2] = new TF1(Form("fitFunc_deltaT_Raw_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
            histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QNRS");
            histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QSR+","",fitFunc_deltaT_Raw_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2),fitFunc_deltaT_Raw_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2));
            histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QSR+","",fitFunc_deltaT_Raw_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2),fitFunc_deltaT_Raw_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2));

            fitFunc_deltaT_Raw_L[index2] -> SetLineColor(kBlack);
            fitFunc_deltaT_Raw_L[index2] -> SetLineWidth(2);
            fitFunc_deltaT_Raw_L[index2] -> Draw("same");

            */
            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), histo->GetMean(), histo->GetRMS()));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

            c -> Print(Form("%s/Loop2_EnergyRatio_REF_TH1/energyRatioR_REF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop2_EnergyRatio_REF_TH1/energyRatioR_REF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete c;
            delete latex;


            } // --- end loop over energy bins

        } // --- end loop ober bars


      int indexREF( (10000*(Vov*100.)) + (100*vth1)); 
      std::string label_ext(Form("barExt_%s",stepLabel.c_str()));            
      
      if(!h1_DeltaTimeLR_REF[indexREF]) continue;
      c = new TCanvas(Form("c_DeltaTimeLR_REF_%s",label_ext.c_str()),Form("c_DeltaTimeLR_REF_%s",label_ext.c_str()));
      histo = h1_DeltaTimeLR_REF[indexREF];
      histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
      histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
      histo -> SetTitle(Form(";t_{L} - t_{R} [ps];entries"));
      histo -> SetLineColor(kRed);
      histo -> SetLineWidth(2);
      //c -> SetLogy();
      histo -> Draw();
      histo -> Write();

            /*
            fitFunc_deltaT_Raw_L[index2] = new TF1(Form("fitFunc_deltaT_Raw_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
            histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QNRS");
            histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QSR+","",fitFunc_deltaT_Raw_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2),fitFunc_deltaT_Raw_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2));
            histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QSR+","",fitFunc_deltaT_Raw_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2),fitFunc_deltaT_Raw_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2));

            fitFunc_deltaT_Raw_L[index2] -> SetLineColor(kBlack);
            fitFunc_deltaT_Raw_L[index2] -> SetLineWidth(2);
            fitFunc_deltaT_Raw_L[index2] -> Draw("same");

            */
            
      latex = new TLatex(0.40,0.85,Form("#splitline{bar 7 REF}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}", Vov, int(vth1), histo->GetMean(), histo->GetRMS()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");

      c -> Print(Form("%s/Loop2_DeltaTimeLR_REF/DeltaTimeLR_REF_%s.png",plotDir.c_str(),label_ext.c_str()));
      c -> Print(Form("%s/Loop2_DeltaTimeLR_REF/DeltaTimeLR_REF_%s.pdf",plotDir.c_str(),label_ext.c_str()));
      delete c;
      delete latex;



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

            /*
            fitFunc_deltaT_Raw_L[index2] = new TF1(Form("fitFunc_deltaT_Raw_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
            histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QNRS");
            histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QSR+","",fitFunc_deltaT_Raw_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2),fitFunc_deltaT_Raw_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2));
            histo -> Fit(fitFunc_deltaT_Raw_L[index2],"QSR+","",fitFunc_deltaT_Raw_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2),fitFunc_deltaT_Raw_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_Raw_L[index2]->GetParameter(2));

            fitFunc_deltaT_Raw_L[index2] -> SetLineColor(kBlack);
            fitFunc_deltaT_Raw_L[index2] -> SetLineWidth(2);
            fitFunc_deltaT_Raw_L[index2] -> Draw("same");

            */

      latex = new TLatex(0.40,0.85,Form("#splitline{bar 7 REF}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}", Vov, int(vth1), histo->GetMean(), histo->GetRMS()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
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

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );


          long long timeMean_ext = 0.5 * (anEvent->timeL_ext + anEvent->timeR_ext);
          long long deltaT_L_raw = anEvent->timeL - timeMean_ext;
          long long deltaT_R_raw = anEvent->timeR - timeMean_ext;

	  float energyMeanREF = 0.5 *(anEvent->energyL_ext + anEvent->energyR_ext);


	  float deltaT_L_raw_low = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
	  float deltaT_L_raw_hig = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);

	  float deltaT_R_raw_low = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_R_raw_hig = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);

          if( p1_deltaT_L_raw_vs_energyL[index2] == NULL )
            {
              std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
              
	      p1_deltaT_L_raw_vs_energyL[index2] = new TProfile(Form("p1_deltaT_L_raw_vs_energyL_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800.);
	      //p1_deltaT_L_raw_vs_energyL[index2]->SetErrorOption("s"); //std DEv
	      p1_deltaT_R_raw_vs_energyR[index2] = new TProfile(Form("p1_deltaT_R_raw_vs_energyR_%s",labelLR_energyBin.c_str()),"",50,ranges["R"][index1]->at(0),800.);
              //p1_deltaT_R_raw_vs_energyR[index2]->SetErrorOption("s"); //Std DEV
	      
	      h2_deltaT_L_raw_vs_energyL[index2] = new TH2F(Form("h2_deltaT_L_raw_vs_energyL_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800., 2000, -12000., 12000.);
	      h2_deltaT_R_raw_vs_energyR[index2] = new TH2F(Form("h2_deltaT_R_raw_vs_energyR_%s",labelLR_energyBin.c_str()),"",50,ranges["R"][index1]->at(0),800., 2000, -12000., 12000.);

	      h2_deltaT_L_raw_vs_MeanEnergyREF[index2] = new TH2F(Form("h2_deltaT_L_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),"",50,200.,800., 2000, -12000., 12000.);
              h2_deltaT_R_raw_vs_MeanEnergyREF[index2] = new TH2F(Form("h2_deltaT_R_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),"",50,200.,800., 2000, -12000., 12000.);
	      p1_deltaT_L_raw_vs_MeanEnergyREF[index2] = new TProfile(Form("p1_deltaT_L_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
	      p1_deltaT_R_raw_vs_MeanEnergyREF[index2] = new TProfile(Form("p1_deltaT_R_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);


	      h2_EneDUT_L_vs_MeanEneREF[index2] = new TH2F(Form("h2_EneDUT_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",80,ranges["L"][index1]->at(0),800.,80,200.,800.);
              h2_EneDUT_R_vs_MeanEneREF[index2] = new TH2F(Form("h2_EneDUT_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",80,ranges["R"][index1]->at(0),800.,80,200.,800.);

	    }


	  if(p1_deltaT_L_raw_vs_energyRatioREF[index2] == NULL ) 
            {

		    std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
      		    p1_deltaT_L_raw_vs_energyRatioREF[index2] = new TProfile(Form("p1_deltaT_L_raw_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5.);
		    h2_deltaT_L_raw_vs_energyRatioREF[index2] = new TH2F(Form("h2_deltaT_L_raw_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);
	    }	      

	  if(p1_deltaT_R_raw_vs_energyRatioREF[index2] == NULL )
            {
              std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
		    p1_deltaT_R_raw_vs_energyRatioREF[index2] = new TProfile(Form("p1_deltaT_R_raw_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5.);
		    h2_deltaT_R_raw_vs_energyRatioREF[index2] = new TH2F(Form("h2_deltaT_R_raw_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);
            }

	  if((deltaT_L_raw > deltaT_L_raw_low ) && (deltaT_L_raw < deltaT_L_raw_hig)){
          p1_deltaT_L_raw_vs_energyL[index2] -> Fill( anEvent->energyL, deltaT_L_raw );
	  h2_deltaT_L_raw_vs_energyL[index2] -> Fill( anEvent->energyL, deltaT_L_raw );

	  p1_deltaT_L_raw_vs_energyRatioREF[index2] -> Fill( (anEvent->energyL)/(energyMeanREF), deltaT_L_raw);
	  h2_deltaT_L_raw_vs_energyRatioREF[index2] -> Fill( (anEvent->energyL)/(energyMeanREF), deltaT_L_raw );
	  
	  p1_deltaT_L_raw_vs_MeanEnergyREF[index2] -> Fill( energyMeanREF, deltaT_L_raw);
          h2_deltaT_L_raw_vs_MeanEnergyREF[index2] -> Fill( energyMeanREF, deltaT_L_raw);
	  
	  h2_EneDUT_L_vs_MeanEneREF[index2] -> Fill (anEvent->energyL,energyMeanREF);
	  }
	  
	  if((deltaT_R_raw > deltaT_R_raw_low ) && (deltaT_R_raw < deltaT_R_raw_hig)){
	  p1_deltaT_R_raw_vs_energyR[index2] -> Fill( anEvent->energyR,deltaT_R_raw );
	  h2_deltaT_R_raw_vs_energyR[index2] -> Fill( anEvent->energyR,deltaT_R_raw );
	  
	  p1_deltaT_R_raw_vs_energyRatioREF[index2] -> Fill( (anEvent->energyR)/(energyMeanREF),deltaT_R_raw);
	  h2_deltaT_R_raw_vs_energyRatioREF[index2] -> Fill( (anEvent->energyR)/(energyMeanREF),deltaT_R_raw );
	  
          p1_deltaT_R_raw_vs_MeanEnergyREF[index2] -> Fill( energyMeanREF, deltaT_R_raw);
          h2_deltaT_R_raw_vs_MeanEnergyREF[index2] -> Fill( energyMeanREF, deltaT_R_raw);

	  h2_EneDUT_R_vs_MeanEneREF[index2] -> Fill(anEvent->energyR,energyMeanREF);
	  }



          }


      std::cout << std::endl;
    }





  //-----------------------
  //--- draw 3rd loop plots

  std::map<double,TF1*> fitFunc_deltaT_EneLCorr_L;
  std::map<double,TF1*> fitFunc_deltaT_EneRCorr_R;
  std::map<double,TF1*> fitFunc_deltaT_L_raw_ERatioREF;
  std::map<double,TF1*> fitFunc_deltaT_R_raw_ERatioREF;


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
            
	    
	    if(!p1_deltaT_L_raw_vs_energyL[index2]) continue;

            std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));


	    // -- Profile Delta T vs Energy
	    
	    //Ch L:
	    c = new TCanvas(Form("c_deltaT_L_raw_vs_energyL_%s",labelLR_energyBin.c_str()),Form("c_deltaT_L_raw_vs_energyL_%s",labelLR_energyBin.c_str()));

            prof = p1_deltaT_L_raw_vs_energyL[index2];
            prof -> SetTitle(Form(";Energy (Ch L);#Delta T_{L}[ps]"));
            prof -> GetYaxis() -> SetRangeUser(CTRMeans_L[index2]-5.*CTRSigmas_L[index2],CTRMeans_L[index2]+5.*CTRSigmas_L[index2]);
            prof -> Draw("");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

	    float fitXMin_L = ranges["L"][index1]->at(0);
	    float fitXMax_L = ranges["L"][index1]->at(1);

	    fitFunc_deltaT_EneLCorr_L[index2] = new TF1(Form("fitFunc_deltaT_EneLCorr_L_%s",labelLR_energyBin.c_str()),"pol3",fitXMin_L,fitXMax_L);
	    prof -> Fit(fitFunc_deltaT_EneLCorr_L[index2],"QRS+");
	    fitFunc_deltaT_EneLCorr_L[index2] -> SetLineColor(kRed);
	    fitFunc_deltaT_EneLCorr_L[index2] -> SetLineWidth(2);
	    fitFunc_deltaT_EneLCorr_L[index2] -> Draw("same");

	    c -> Print(Form("%s/Loop3_deltaT_raw_vs_Energy/c_deltaT_L_raw_vs_energyL__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop3_deltaT_raw_vs_Energy/c_deltaT_L_raw_vs_energyL__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;
	    
	    
	    /*	    
	    // -------------------------------------------------------------
	    if (iBar ==8 && vth1 ==10){

	    for (int bin = 1; bin <= prof->GetNbinsX(); ++bin){
            double x_center  = prof->GetBinCenter(bin);
            //double x_min     = prof->GetBinLowEdge(bin);
            //double x_max     = x_min + prof->GetBinWidth(bin);
            double y_mean    = prof->GetBinContent(bin);
            double y_err     = prof->GetBinError(bin);      // errore sulla media
            double entries   = prof->GetBinEntries(bin);

	    std::cout << std::fixed << std::setprecision(3);
	    std::cout
	    << "vth="     << std::setw(3)  << vth1
	    << " | bar="  << std::setw(2)  << iBar
	    << " | bin="  << std::setw(3)  << bin
	    << " | Xcenter=" << std::setw(8) << x_center
	    << " | entries=" << std::setw(6) << entries
	    << " | mean="    << std::setw(8) << y_mean
	    << " | err="     << std::setw(8) << y_err
	    << std::endl;
	    }	    
	    }
            */


	    //Ch R:
            c = new TCanvas(Form("c_deltaT_R_raw_vs_energyR_%s",labelLR_energyBin.c_str()),Form("c_deltaT_R_raw_vs_energyR_%s",labelLR_energyBin.c_str()));

	    prof = p1_deltaT_R_raw_vs_energyR[index2];
            prof -> SetTitle(Form(";Energy (Ch R);#Delta T_{R}[ps]"));
            prof -> GetYaxis() -> SetRangeUser(CTRMeans_R[index2]-5.*CTRSigmas_R[index2],CTRMeans_R[index2]+5.*CTRSigmas_R[index2]);
            prof -> Draw("");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

	    float fitXMin_R = ranges["R"][index1]->at(0);
            float fitXMax_R = ranges["R"][index1]->at(1);

            fitFunc_deltaT_EneRCorr_R[index2] = new TF1(Form("fitFunc_deltaT_EneRCorr_R_%s",labelLR_energyBin.c_str()),"pol3",fitXMin_R,fitXMax_R);
            prof -> Fit(fitFunc_deltaT_EneRCorr_R[index2],"QRS+");
            fitFunc_deltaT_EneRCorr_R[index2] -> SetLineColor(kRed);
            fitFunc_deltaT_EneRCorr_R[index2] -> SetLineWidth(2);
            fitFunc_deltaT_EneRCorr_R[index2] -> Draw("same");

            c -> Print(Form("%s/Loop3_deltaT_raw_vs_Energy/c_deltaT_R_raw_vs_energyR__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop3_deltaT_raw_vs_Energy/c_deltaT_R_raw_vs_energyR__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;




	    // -- profile + TH2F deltaT vs energy (No TW corr)

	    //Ch L:
            c = new TCanvas(Form("c_deltaT_L_raw_vs_energyL_scatter_%s",labelLR_energyBin.c_str()),Form("c_deltaT_L_raw_vs_energyL_scatter_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_L_raw_vs_energyL[index2];
            h2->SetTitle(";Energy (Ch L) [a.u.];#Delta T_{L} [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_L_raw_vs_energyL[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop3_deltaT_raw_vs_Energy_scatter/c_deltaT_L_raw_vs_energyL__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop3_deltaT_raw_vs_Energy_scatter/c_deltaT_L_raw_vs_energyL__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    //ChR
	    c = new TCanvas(Form("c_deltaT_R_raw_vs_energyR_scatter_%s",labelLR_energyBin.c_str()),Form("c_deltaT_R_raw_vs_energyR_scatter_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_R_raw_vs_energyR[index2];
            h2->SetTitle(";Energy (Ch R) [a.u.];#Delta T_{R} [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_R_raw_vs_energyR[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop3_deltaT_raw_vs_Energy_scatter/c_deltaT_R_raw_vs_energyR__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop3_deltaT_raw_vs_Energy_scatter/c_deltaT_R_raw_vs_energyR__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;



	    // -- DeltaT vs energy Ratio [Energy DUT]/[Mean Energy REF bar]
	    
	    // Ch L
            c = new TCanvas(Form("c_deltaT_L_raw_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_L_raw_vs_energyRatioREF_%s",labelLR_energyBin.c_str()));

            prof = p1_deltaT_L_raw_vs_energyRatioREF[index2];
            prof -> SetTitle(Form(";E_{L}/E_{REF};#Delta T_{L}[ps]"));
	    prof -> GetYaxis() -> SetRangeUser(CTRMeans_L[index2]-5.*CTRSigmas_L[index2],CTRMeans_L[index2]+5.*CTRSigmas_L[index2]);
            prof -> GetXaxis() -> SetRangeUser(0.,3.5);
            prof -> Draw("");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

            //float fitXMin_EneRatio_L = 0.3;
            //float fitXMax_EneRatio_L = 3.5;
	   
	    auto range_L = GetDynamicFitRange(prof, 0.05);
	    float fitXMin_EneRatio_L = range_L.first;
	    float fitXMax_EneRatio_L = range_L.second; 

            fitFunc_deltaT_L_raw_ERatioREF[index2] = new TF1(Form("fitFunc_deltaT_raw_L_ERatioREF_%s",labelLR_energyBin.c_str()),"pol3",fitXMin_EneRatio_L,fitXMax_EneRatio_L);
            prof -> Fit(fitFunc_deltaT_L_raw_ERatioREF[index2],"QRS+");
            fitFunc_deltaT_L_raw_ERatioREF[index2] -> SetLineColor(kRed);
            fitFunc_deltaT_L_raw_ERatioREF[index2] -> SetLineWidth(2);
            fitFunc_deltaT_L_raw_ERatioREF[index2] -> Draw("same");
	     
            c -> Print(Form("%s/Loop3_DeltaT_raw_vs_EnergyRatioREF/c_deltaT_L_raw_vs_energyRatioREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop3_DeltaT_raw_vs_EnergyRatioREF/c_deltaT_L_raw_vs_energyRatioREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    // Ch R
            c = new TCanvas(Form("c_deltaT_R_raw_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_R_raw_vs_energyRatioREF_%s",labelLR_energyBin.c_str()));

            prof = p1_deltaT_R_raw_vs_energyRatioREF[index2];
            prof -> SetTitle(Form(";E_{R}/E_{REF};#Delta T_{R}[ps]"));
            prof -> GetXaxis() -> SetRangeUser(0.,3.5);
	    prof -> GetYaxis() -> SetRangeUser(CTRMeans_R[index2]-5.*CTRSigmas_R[index2],CTRMeans_R[index2]+5.*CTRSigmas_R[index2]);
            prof -> Draw("");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

            //float fitXMin_EneRatio_R = 0.3;
            //float fitXMax_EneRatio_R = 3.5;

	    auto range_R = GetDynamicFitRange(prof, 0.05);
            float fitXMin_EneRatio_R = range_R.first;
            float fitXMax_EneRatio_R = range_R.second;

            fitFunc_deltaT_R_raw_ERatioREF[index2] = new TF1(Form("fitFunc_deltaT_R_raw_ERatioREF_%s",labelLR_energyBin.c_str()),"pol3",fitXMin_EneRatio_R,fitXMax_EneRatio_R);
            prof -> Fit(fitFunc_deltaT_R_raw_ERatioREF[index2],"QRS+");
            fitFunc_deltaT_R_raw_ERatioREF[index2] -> SetLineColor(kRed);
            fitFunc_deltaT_R_raw_ERatioREF[index2] -> SetLineWidth(2);
            fitFunc_deltaT_R_raw_ERatioREF[index2] -> Draw("same");

            c -> Print(Form("%s/Loop3_DeltaT_raw_vs_EnergyRatioREF/c_p1_deltaT_R_vs_energyRatioR_REF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop3_DeltaT_raw_vs_EnergyRatioREF/c_p1_deltaT_R_vs_energyRatioR_REF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;






	    // -- profile + TH2F deltaT vs energy Ratio [Energy DUT]/[Mean Energy REF bar]

            //Ch L:
            c = new TCanvas(Form("c_deltaT_L_raw_vs_energyRatioREF_scatter_%s",labelLR_energyBin.c_str()),Form("c_deltaT_L_raw_vs_energyRatioREF_scatter_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_L_raw_vs_energyRatioREF[index2];
            h2->SetTitle(";E_{L}/E_{REF};#Delta T_{L}[ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            //h2 -> GetYaxis() -> SetRangeUser(CTRMeans_L[index2]-5.*CTRSigmas_L[index2],CTRMeans_L[index2]+5.*CTRSigmas_L[index2]);
	    h2->GetXaxis()->SetRangeUser(0., 3.5);
	    h2->Draw("colz");

            prof = p1_deltaT_L_raw_vs_energyRatioREF[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop3_DeltaT_raw_vs_EnergyRatioREF_scatter/c_deltaT_L_raw_vs_energyRatioREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop3_DeltaT_raw_vs_EnergyRatioREF_scatter/c_deltaT_L_raw_vs_energyRatioREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    //Ch R:
            c = new TCanvas(Form("c_deltaT_R_raw_vs_energyRatioREF_scatter_%s",labelLR_energyBin.c_str()),Form("c_deltaT_R_raw_vs_energyRatioREF_scatter_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_R_raw_vs_energyRatioREF[index2];
            h2->SetTitle(";E_{R}/E_{REF};#Delta T_{R}[ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            //h2 -> GetYaxis() -> SetRangeUser(CTRMeans_R[index2]-5.*CTRSigmas_R[index2],CTRMeans_R[index2]+5.*CTRSigmas_R[index2]);
	    h2->GetXaxis()->SetRangeUser(0., 3.5);
	    h2->Draw("colz");

            prof = p1_deltaT_R_raw_vs_energyRatioREF[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop3_DeltaT_raw_vs_EnergyRatioREF_scatter/c_deltaT_R_raw_vs_energyRatioREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop3_DeltaT_raw_vs_EnergyRatioREF_scatter/c_deltaT_R_raw_vs_energyRatioREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;






            // -- profile + TH2F deltaT raw vs Mean Energy REF bar 7

            //Ch L:
            c = new TCanvas(Form("c_p1_deltaT_L_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),Form("c_p1_deltaT_L_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_L_raw_vs_MeanEnergyREF[index2];
            h2->SetTitle(";<E_{REF}>;#Delta T_{L}[ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            //h2 -> GetYaxis() -> SetRangeUser(CTRMeans_L[index2]-5.*CTRSigmas_L[index2],CTRMeans_L[index2]+5.*CTRSigmas_L[index2]);
            //h2->GetXaxis()->SetRangeUser(0., 3.5);
            h2->Draw("colz");

            prof = p1_deltaT_L_raw_vs_MeanEnergyREF[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop3_deltaT_raw_vs_MeanEnergyREF/c_deltaT_raw_vs_MeanEnergyREF_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop3_deltaT_raw_vs_MeanEnergyREF/c_deltaT_raw_vs_MeanEnergyREF_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;



            //Ch R:
            c = new TCanvas(Form("c_p1_deltaT_R_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()),Form("c_p1_deltaT_R_raw_vs_MeanEnergyREF_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_R_raw_vs_MeanEnergyREF[index2];
            h2->SetTitle(";<E_{REF}>;#Delta T_{R}[ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            //h2 -> GetYaxis() -> SetRangeUser(CTRMeans_L[index2]-5.*CTRSigmas_L[index2],CTRMeans_L[index2]+5.*CTRSigmas_L[index2]);
            //h2->GetXaxis()->SetRangeUser(0., 3.5);
            h2->Draw("colz");

            prof = p1_deltaT_R_raw_vs_MeanEnergyREF[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop3_deltaT_raw_vs_MeanEnergyREF/c_deltaT_raw_vs_MeanEnergyREF_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop3_deltaT_raw_vs_MeanEnergyREF/c_deltaT_raw_vs_MeanEnergyREF_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;





            // -- TH2F Energy Ch DUT vs Mean Energy REF bar 7

            //Ch L:
            c = new TCanvas(Form("c_h2_EneDUT_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),Form("c_h2_EneDUT_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_EneDUT_L_vs_MeanEneREF[index2];
            h2->SetTitle(";E_{L} [a.u.];<E_{REF}> [a.u.]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            //h2 -> GetYaxis() -> SetRangeUser(CTRMeans_L[index2]-5.*CTRSigmas_L[index2],CTRMeans_L[index2]+5.*CTRSigmas_L[index2]);
            //h2->GetXaxis()->SetRangeUser(0., 3.5);
            h2->Draw("colz");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop3_EneDUT_vs_MeanEneREF/c_EneDUT_L_vs_MeanEneREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop3_EneDUT_vs_MeanEneREF/c_EneDUT_L_vs_MeanEneREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


            //Ch R:
            c = new TCanvas(Form("c_h2_EneDUT_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),Form("c_h2_EneDUT_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_EneDUT_R_vs_MeanEneREF[index2];
            h2->SetTitle(";E_{R} [a.u.];<E_{REF}> [a.u.]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            //h2 -> GetYaxis() -> SetRangeUser(CTRMeans_L[index2]-5.*CTRSigmas_L[index2],CTRMeans_L[index2]+5.*CTRSigmas_L[index2]);
            //h2->GetXaxis()->SetRangeUser(0., 3.5);
            h2->Draw("colz");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop3_EneDUT_vs_MeanEneREF/c_EneDUT_R_vs_MeanEneREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop3_EneDUT_vs_MeanEneREF/c_EneDUT_R_vs_MeanEneREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


          }
      }
    }










  
  //------------------------
  //--- 4th loop over events
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

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );
          double indexBarID( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) );
	  
	  long long timeMean_ext = 0.5 * (anEvent->timeL_ext + anEvent->timeR_ext);
	  long long deltaT_L_raw = anEvent->timeL - timeMean_ext;
          long long deltaT_R_raw = anEvent->timeR - timeMean_ext;
          //float t1fineMean = 0.5 * ( anEvent->t1fineR + anEvent->t1fineL );
          float energyMeanREF = 0.5 *(anEvent->energyL_ext + anEvent->energyR_ext);

	  float deltaT_L_raw_low = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_L_raw_hig = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);

	  float deltaT_R_raw_low = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_R_raw_hig = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);

	  
	  // -- TW corr only fitFunc(energy):
	  float TW_corr_L = fitFunc_deltaT_EneLCorr_L[index2]->Eval(anEvent->energyL);
          float TW_corr_R = fitFunc_deltaT_EneRCorr_R[index2]->Eval(anEvent->energyR);

	  // -- TW corr fitFunc(energy) - fitFunc(MPV energy):
	  //float TW_corr = fitFunc_deltaT_EneLCorr_L[index2]->Eval(anEvent->energyL) - fitFunc_deltaT_EneLCorr_L[index2]->Eval(f_landau[index1]->GetParameter(1));
	  
	  // -- EneRatioREF correction (from delta T raw):
	  float EneRatio_corr_L = fitFunc_deltaT_L_raw_ERatioREF[index2]->Eval(anEvent->energyL/energyMeanREF);
          float EneRatio_corr_R = fitFunc_deltaT_R_raw_ERatioREF[index2]->Eval(anEvent->energyR/energyMeanREF);

	  std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));

	  // -- delta T histograms
          if( h1_deltaT_L_TWcorr[index2] == NULL )
	  {

	    // -- Delta T w/ TW corr  
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



	    // -- Delta T w/ TW corr vs energyRatioREF
	    p1_deltaT_TWCorr_L_vs_energyRatioREF[index2] = new TProfile(Form("p1_deltaT_TWCorr_L_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5.);
	    p1_deltaT_TWCorr_R_vs_energyRatioREF[index2] = new TProfile(Form("p1_deltaT_TWCorr_R_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5.);

	    h2_deltaT_TWCorr_L_vs_energyRatioREF[index2] = new TH2F(Form("h2_deltaT_TWCorr_L_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);
	    h2_deltaT_TWCorr_R_vs_energyRatioREF[index2] = new TH2F(Form("h2_deltaT_TWCorr_R_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);



	    // -- Delta T w/ EneRatioREF corr
	    h1_deltaT_L_EneRatioREFcorr[index2] = new TH1F(Form("h1_deltaT_L_EneRatioREFcorr_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
            h1_deltaT_R_EneRatioREFcorr[index2] = new TH1F(Form("h1_deltaT_R_EneRatioREFcorr_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);

	    p1_deltaT_EneRatioREFCorr_L_vs_energy[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_L_vs_energy_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800.);
	    p1_deltaT_EneRatioREFCorr_R_vs_energy[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_R_vs_energy_%s",labelLR_energyBin.c_str()),"",50,ranges["R"][index1]->at(0),800.);
	  
	    h2_deltaT_EneRatioREFCorr_L_vs_energy[index2] = new TH2F(Form("h2_deltaT_TWCorr_vs_energy_L_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800., 2000, -12000., 12000.);
	    h2_deltaT_EneRatioREFCorr_R_vs_energy[index2] = new TH2F(Form("h2_deltaT_TWCorr_vs_energy_R_%s",labelLR_energyBin.c_str()),"",50,ranges["R"][index1]->at(0),800., 2000, -12000., 12000.);
            
	    p1_deltaT_EneRatioREFCorr_L_vs_energyRatioREF[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_L_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5.);
            p1_deltaT_EneRatioREFCorr_R_vs_energyRatioREF[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_R_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5.);

            h2_deltaT_EneRatioREFCorr_L_vs_energyRatioREF[index2] = new TH2F(Form("h2_deltaT_EneRatioREFCorr_L_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);
            h2_deltaT_EneRatioREFCorr_R_vs_energyRatioREF[index2] = new TH2F(Form("h2_deltaT_EneRatioREFCorr_R_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);

            p1_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_L_vs_t1fineL_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
            p1_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_R_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
	    
	    h2_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2] = new TH2F(Form("h2_deltaT_EneRatioREFCorr_L_vs_t1fineL_%s",labelLR_energyBin.c_str()),"",50,0.,1000., 2000, -12000., 12000.);
	    h2_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2] = new TH2F(Form("h2_deltaT_EneRatioREFCorr_R_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0.,1000., 2000, -12000., 12000.);
	   
	  
	    // -- Delta corr (TW or EneRatio Corr) vs Mean Energy REF bar

            p1_deltaT_EneRatioREFCorr_L_vs_MeanEneREF[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
            p1_deltaT_EneRatioREFCorr_R_vs_MeanEneREF[index2] = new TProfile(Form("p1_deltaT_EneRatioREFCorr_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
            h2_deltaT_EneRatioREFCorr_L_vs_MeanEneREF[index2] = new TH2F(Form("h2_deltaT_EneRatioREFCorr_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.,2000, -12000., 12000.);
            h2_deltaT_EneRatioREFCorr_R_vs_MeanEneREF[index2] = new TH2F(Form("h2_deltaT_EneRatioREFCorr_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.,2000, -12000., 12000.);

	    p1_deltaT_TWCorr_L_vs_MeanEneREF[index2] = new TProfile(Form("p1_deltaT_TWCorr_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
            p1_deltaT_TWCorr_R_vs_MeanEneREF[index2] = new TProfile(Form("p1_deltaT_TWCorr_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.);
            h2_deltaT_TWCorr_L_vs_MeanEneREF[index2] = new TH2F(Form("h2_deltaT_TWCorr_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.,2000, -12000., 12000.);
            h2_deltaT_TWCorr_R_vs_MeanEneREF[index2] = new TH2F(Form("h2_deltaT_TWCorr_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),"",50,200.,800.,2000, -12000., 12000.);


	    //Scatter plot t1Fine_L vs t1Fine_R:
	    p1_t1fineL_vs_t1fineR[index2] = new TProfile(Form("p1_t1fineL_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0.,1000.);
            h2_t1fineL_vs_t1fineR[index2]= new TH2F(Form("h2_t1fineL_vs_t1fineR_%s",labelLR_energyBin.c_str()),"",50,0.,1000.,50,0.,1000.);

	  
	  }

	  
	  if((deltaT_L_raw > deltaT_L_raw_low ) && (deltaT_L_raw < deltaT_L_raw_hig)){
          // -- Delta T w/ TW corr  
	  h1_deltaT_L_TWcorr[index2] -> Fill(deltaT_L_raw - TW_corr_L);
	  h2_deltaT_TWCorr_vs_energy_L[index2] -> Fill( anEvent->energyL, deltaT_L_raw - TW_corr_L );
	  p1_deltaT_TWCorr_vs_energy_L[index2] -> Fill( anEvent->energyL, deltaT_L_raw - TW_corr_L );
	  
	  p1_deltaT_TWCorrL_vs_t1fineL[index2] -> Fill(anEvent->t1fineL, deltaT_L_raw - TW_corr_L);
          h2_deltaT_TWCorrL_vs_t1fineL[index2] -> Fill(anEvent->t1fineL, deltaT_L_raw - TW_corr_L);

	  p1_deltaT_TWCorr_L_vs_energyRatioREF[index2] -> Fill ( anEvent->energyL/energyMeanREF, deltaT_L_raw - TW_corr_L );
          h2_deltaT_TWCorr_L_vs_energyRatioREF[index2] -> Fill ( anEvent->energyL/energyMeanREF, deltaT_L_raw - TW_corr_L );

	  //Scatter plot t1Fine_L vs t1Fine_R:
	  p1_t1fineL_vs_t1fineR[index2] -> Fill ( anEvent->t1fineL, anEvent->t1fineR );
	  h2_t1fineL_vs_t1fineR[index2] -> Fill(anEvent->t1fineL, anEvent->t1fineR);


	  // -- Delta T w/ EneRatioREF corr
	  h1_deltaT_L_EneRatioREFcorr[index2] -> Fill ( deltaT_L_raw - EneRatio_corr_L );
	  p1_deltaT_EneRatioREFCorr_L_vs_energy[index2] -> Fill ( anEvent->energyL, deltaT_L_raw - EneRatio_corr_L );
	  h2_deltaT_EneRatioREFCorr_L_vs_energy[index2] -> Fill ( anEvent->energyL, deltaT_L_raw - EneRatio_corr_L );
	  p1_deltaT_EneRatioREFCorr_L_vs_energyRatioREF[index2] -> Fill ( anEvent->energyL/energyMeanREF, deltaT_L_raw - EneRatio_corr_L );
	  h2_deltaT_EneRatioREFCorr_L_vs_energyRatioREF[index2] -> Fill ( anEvent->energyL/energyMeanREF, deltaT_L_raw - EneRatio_corr_L );
	  p1_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2] -> Fill(anEvent->t1fineL, deltaT_L_raw - EneRatio_corr_L);
	  h2_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2] -> Fill(anEvent->t1fineL, deltaT_L_raw - EneRatio_corr_L);


	  p1_deltaT_TWCorr_L_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_L_raw - TW_corr_L);
	  h2_deltaT_TWCorr_L_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_L_raw - TW_corr_L);
          p1_deltaT_EneRatioREFCorr_L_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_L_raw - EneRatio_corr_L);
          h2_deltaT_EneRatioREFCorr_L_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_L_raw - EneRatio_corr_L);	  
	  }


	  if((deltaT_R_raw > deltaT_R_raw_low ) && (deltaT_R_raw < deltaT_R_raw_hig)){	  
          // -- Delta T w/ TW corr		  
	  h1_deltaT_R_TWcorr[index2] -> Fill(deltaT_R_raw - TW_corr_R);
	  p1_deltaT_TWCorr_vs_energy_R[index2] -> Fill( anEvent->energyR, deltaT_R_raw - TW_corr_R );
          h2_deltaT_TWCorr_vs_energy_R[index2] -> Fill( anEvent->energyR, deltaT_R_raw - TW_corr_R );
	  
	  p1_deltaT_TWCorrR_vs_t1fineR[index2] -> Fill(anEvent->t1fineR, deltaT_R_raw - TW_corr_R);
	  h2_deltaT_TWCorrR_vs_t1fineR[index2] -> Fill(anEvent->t1fineR, deltaT_R_raw - TW_corr_R);

	  p1_deltaT_TWCorr_R_vs_energyRatioREF[index2] -> Fill ( anEvent->energyR/energyMeanREF, deltaT_R_raw - TW_corr_R );
	  h2_deltaT_TWCorr_R_vs_energyRatioREF[index2] -> Fill ( anEvent->energyR/energyMeanREF, deltaT_R_raw - TW_corr_R );
	  
	  // -- Delta T w/ EneRatioREF corr
	  h1_deltaT_R_EneRatioREFcorr[index2] -> Fill ( deltaT_R_raw - EneRatio_corr_R );
	  p1_deltaT_EneRatioREFCorr_R_vs_energy[index2] -> Fill ( anEvent->energyR, deltaT_R_raw - EneRatio_corr_R );
	  h2_deltaT_EneRatioREFCorr_R_vs_energy[index2] -> Fill ( anEvent->energyR, deltaT_R_raw - EneRatio_corr_R );
	  p1_deltaT_EneRatioREFCorr_R_vs_energyRatioREF[index2] -> Fill ( anEvent->energyR/energyMeanREF, deltaT_R_raw - EneRatio_corr_R );
          h2_deltaT_EneRatioREFCorr_R_vs_energyRatioREF[index2] -> Fill ( anEvent->energyR/energyMeanREF, deltaT_R_raw - EneRatio_corr_R );
	  p1_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2] -> Fill(anEvent->t1fineR, deltaT_R_raw - EneRatio_corr_R);
	  h2_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2] -> Fill(anEvent->t1fineR, deltaT_R_raw - EneRatio_corr_R);
	  
          p1_deltaT_TWCorr_R_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_R_raw - TW_corr_R);
          h2_deltaT_TWCorr_R_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_R_raw - TW_corr_R);
          p1_deltaT_EneRatioREFCorr_R_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_R_raw - EneRatio_corr_R);
          h2_deltaT_EneRatioREFCorr_R_vs_MeanEneREF[index2] -> Fill(energyMeanREF,deltaT_R_raw - EneRatio_corr_R);	  
	  }


	  //----------------------------------------------------------------------
	  // -- Test high deltaT > 150
	  if( h2_energy_vs_deltaT_high_L[index2] == NULL ){

		  h2_energy_vs_deltaT_high_L[index2] = new TH2F(Form("h2_energy_vs_deltaT_high_L_%s",labelLR_energyBin.c_str()),"",250, 0., 3000., 90,ranges["L"][index1]->at(0),800.);
		  h2_energy_vs_deltaT_high_R[index2] = new TH2F(Form("h2_energy_vs_deltaT_high_R_%s",labelLR_energyBin.c_str()),"",250, 0., 3000., 90,ranges["L"][index1]->at(0),800.);
                  h1_energyL_high_deltaT[index2] = new TH1F(Form("h1_energyL_high_deltaT_%s",labelLR_energyBin.c_str()),"",512,0,1024);
	          h1_energyR_high_deltaT[index2] = new TH1F(Form("h1_energyR_high_deltaT_%s",labelLR_energyBin.c_str()),"",512,0,1024);
	  }

	  //if( h1_bar_high_deltaT[indexBarID] == NULL ){
          //        h1_bar_high_deltaT[indexBarID] = new TH1F(Form("h1_bar_high_deltaT_%s",labelLR_energyBin.c_str()),"",16,-0.5,15.5);
          //}

	  if((deltaT_L_raw - TW_corr_L) > 150.) {
		  h1_energyL_high_deltaT[index2] -> Fill(anEvent->energyL);
		  h1_energyR_high_deltaT[index2] -> Fill(anEvent->energyR);
		  h2_energy_vs_deltaT_high_L[index2] -> Fill (deltaT_L_raw - TW_corr_L, anEvent->energyL);
		  h2_energy_vs_deltaT_high_R[index2] -> Fill (deltaT_R_raw - TW_corr_R, anEvent->energyR);

		  //h1_bar_high_deltaT[indexBarID] -> Fill (anEvent->barID);
		  //std::cout << " | Delta T (> 150 ps): " << deltaT_L_ext - TW_corr_L << "Ene Ratio REF" << anEvent->energyL/energyMeanREF << " | Bar DUT: " << anEvent->barID << std::endl; 

	  }
	  //----------------------------------------------------------------------

	  
	}
    }




  //-----------------------
  //--- draw plots loop 4th

  std::map<double,TF1*> fitFunc_deltaT_TWcorr_L;
  std::map<double,TF1*> fitFunc_deltaT_TWcorr_R;
  std::map<double,TF1*> fitFunc_deltaT_L_TWCorr_EneRatioREF;
  std::map<double,TF1*> fitFunc_deltaT_R_TWCorr_EneRatioREF;
  std::map<double,TF1*> fitFunc_deltaT_EneRatioREFcorr_L;
  std::map<double,TF1*> fitFunc_deltaT_EneRatioREFcorr_R;
  std::map<double,TF1*> fitFunc_deltaT_EneRatioREF_TW_L;
  std::map<double,TF1*> fitFunc_deltaT_EneRatioREF_TW_R;
  std::map<double,TF1*> fitFunc_deltaT_EneRCorrL_t1fineL;
  std::map<double,TF1*> fitFunc_deltaT_EneRCorrR_t1fineR;
  std::map<double,TF1*> fitFunc_deltaT_TWCorrL_t1fineL;
  std::map<double,TF1*> fitFunc_deltaT_TWCorrR_t1fineR; 

  std::map<double,float> CTRMeans_TW_L;
  std::map<double,float> CTRSigmas_TW_L;

    for(auto mapIt : h1_deltaT_L_TWcorr)
    {
        double index = mapIt.first;

      FindSmallestInterval(vals,h1_deltaT_L_TWcorr[index],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans_TW_L[index] = mean;
      CTRSigmas_TW_L[index] = effSigma;
    }


    std::map<double,float> CTRMeans_TW_R;
    std::map<double,float> CTRSigmas_TW_R;

    for(auto mapIt : h1_deltaT_R_TWcorr)
    {
        double index = mapIt.first;

      FindSmallestInterval(vals,h1_deltaT_R_TWcorr[index],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans_TW_R[index] = mean;
      CTRSigmas_TW_R[index] = effSigma;
    }



    std::map<double,float> CTRMeans_ERatioREF_L;
    std::map<double,float> CTRSigmas_ERatioREF_L;

    for(auto mapIt : h1_deltaT_L_EneRatioREFcorr)
    {
        double index = mapIt.first;

      FindSmallestInterval(vals,h1_deltaT_L_EneRatioREFcorr[index],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans_ERatioREF_L[index] = mean;
      CTRSigmas_ERatioREF_L[index] = effSigma;
    }



     std::map<double,float> CTRMeans_ERatioREF_R;
    std::map<double,float> CTRSigmas_ERatioREF_R;

    for(auto mapIt : h1_deltaT_R_EneRatioREFcorr)
    {
        double index = mapIt.first;

      FindSmallestInterval(vals,h1_deltaT_R_EneRatioREFcorr[index],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans_ERatioREF_R[index] = mean;
      CTRSigmas_ERatioREF_R[index] = effSigma;
    }



  for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth1 = map_ths[stepLabel];

      /*
      //------------------------------------------------------------
      double indexBarID( (10000*int(Vov*100.)) + (100*vth1) );

      // Bar ID high deltaT:
      if (!h1_bar_high_deltaT[indexBarID]) continue;
      std::string labelBarID(Form("%s",stepLabel.c_str()));

      c = new TCanvas(Form("c_h1_bar_high_deltaT_%s",labelBarID.c_str()),Form("c_h1_bar_high_deltaT_%s",labelBarID.c_str()));
      histo = h1_bar_high_deltaT[indexBarID];
      histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
      histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
      histo -> SetTitle(Form(";Bar ID;entries"));
      histo -> SetLineColor(kBlack);
      histo -> SetLineWidth(2);
      //c -> SetLogy();
      histo -> Draw();
      histo -> Write();

      latex = new TLatex(0.40,0.85,Form("#splitline{#Delta T > 150 ps}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}", Vov, int(vth1), histo->GetMean(), histo->GetRMS()));
      latex -> SetNDC();
      latex -> SetTextFont(42);
      latex -> SetTextSize(0.04);
      latex -> SetTextColor(kRed);
      latex -> Draw("same");

      c -> Print(Form("%s/BarID_highDeltaT/BarID_high_deltaT__%s.png",plotDir.c_str(),labelBarID.c_str()));
      c -> Print(Form("%s/BarID_highDeltaT/BarID_high_deltaT__%s.pdf",plotDir.c_str(),labelBarID.c_str()));
      delete c;
      delete latex;
      //-------------------------------------------------------------
      */


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


              if (!h1_deltaT_L_TWcorr[index2]) continue;
              std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));

	      // -- Delta T histograms with Time-Walk correction

	      // Ch L:  
              c = new TCanvas(Form("c_h1_deltaT_L_TWcorr_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_L_TWcorr_%s",labelLR_energyBin.c_str()));
              histo = h1_deltaT_L_TWcorr[index2];
              histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#Delta T_{L} w/ TW [ps];entries"));
              histo -> SetLineColor(kBlue);
              histo -> SetLineWidth(2);
              //c -> SetLogy();
              histo -> Draw();
              histo -> Write();

	      fitFunc_deltaT_TWcorr_L[index2] = new TF1(Form("fitFunc_deltaT_TWcorr_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
              histo -> Fit(fitFunc_deltaT_TWcorr_L[index2],"QNRS");
              histo -> Fit(fitFunc_deltaT_TWcorr_L[index2],"QSR+","",fitFunc_deltaT_TWcorr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_TWcorr_L[index2]->GetParameter(2),fitFunc_deltaT_TWcorr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_TWcorr_L[index2]->GetParameter(2));
              histo -> Fit(fitFunc_deltaT_TWcorr_L[index2],"QSR+","",fitFunc_deltaT_TWcorr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_TWcorr_L[index2]->GetParameter(2),fitFunc_deltaT_TWcorr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_TWcorr_L[index2]->GetParameter(2));

              fitFunc_deltaT_TWcorr_L[index2] -> SetLineColor(kBlack);
              fitFunc_deltaT_TWcorr_L[index2] -> SetLineWidth(2);
              fitFunc_deltaT_TWcorr_L[index2] -> Draw("same");

              latex = new TLatex(0.40,0.85,Form("#splitline{w/ TW bar %02d L}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_TWcorr_L[index2]->GetParameter(1), fitFunc_deltaT_TWcorr_L[index2]->GetParameter(2)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop4_deltaT_TW_corr/deltaT_TWCorr_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop4_deltaT_TW_corr/deltaT_TWCorr_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete c;
              delete latex;
	    
	   
	      //Ch R:
	      c = new TCanvas(Form("c_h1_deltaT_R_TWcorr_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_R_TWcorr_%s",labelLR_energyBin.c_str()));
              histo = h1_deltaT_R_TWcorr[index2];
              histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#Delta T_{R} w/ TW [ps];entries"));
              histo -> SetLineColor(kBlue);
              histo -> SetLineWidth(2);
              //c -> SetLogy();
              histo -> Draw();
              histo -> Write();

              fitFunc_deltaT_TWcorr_R[index2] = new TF1(Form("fitFunc_deltaT_TWcorr_R_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
              histo -> Fit(fitFunc_deltaT_TWcorr_R[index2],"QNRS");
              histo -> Fit(fitFunc_deltaT_TWcorr_R[index2],"QSR+","",fitFunc_deltaT_TWcorr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_TWcorr_R[index2]->GetParameter(2),fitFunc_deltaT_TWcorr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_TWcorr_R[index2]->GetParameter(2));
              histo -> Fit(fitFunc_deltaT_TWcorr_R[index2],"QSR+","",fitFunc_deltaT_TWcorr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_TWcorr_R[index2]->GetParameter(2),fitFunc_deltaT_TWcorr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_TWcorr_R[index2]->GetParameter(2));

              fitFunc_deltaT_TWcorr_R[index2] -> SetLineColor(kBlack);
              fitFunc_deltaT_TWcorr_R[index2] -> SetLineWidth(2);
              fitFunc_deltaT_TWcorr_R[index2] -> Draw("same");

              latex = new TLatex(0.40,0.85,Form("#splitline{w/ TW bar %02d R}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_TWcorr_R[index2]->GetParameter(1), fitFunc_deltaT_TWcorr_R[index2]->GetParameter(2)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop4_deltaT_TW_corr/deltaT_TWCorr_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop4_deltaT_TW_corr/deltaT_TWCorr_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete c;
              delete latex;





	    // -- Profile + TH2F deltaT vs energy (with TW corr)

	    //Ch L:
	    c = new TCanvas(Form("c_deltaT_TWCorr_L_vs_energyL_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_L_vs_energyL_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

	    h2 = h2_deltaT_TWCorr_vs_energy_L[index2];
	    h2->SetTitle(";Energy (Ch L) [a.u.];#Delta T_{L} [ps]");
	    h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
	    h2->Draw("colz");

	    prof = p1_deltaT_TWCorr_vs_energy_L[index2];
	    prof->SetLineColor(kRed);
	    prof->SetLineWidth(3);
	    prof->Draw("same");     

	    latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	    latex->SetNDC();
	    latex->SetTextFont(42);
	    latex->SetTextSize(0.04);
	    latex->SetTextColor(kRed);
	    latex->Draw("same");

    	    c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_Energy/deltaT_TWCorr_vs_energyL__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
       	    c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_Energy/deltaT_TWCorr_vs_energyL__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	    delete latex;
	    delete c;


	    //Ch R:
            c = new TCanvas(Form("c_deltaT_TWCorr_R_vs_energyR_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_R_vs_energyR_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_TWCorr_vs_energy_R[index2];
            h2->SetTitle(";Energy (Ch R) [a.u.];#Delta T_{R} [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_TWCorr_vs_energy_R[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_Energy/deltaT_TWCorr_vs_energyR__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_Energy/deltaT_TWCorr_vs_energyR__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;



	    // -- Test high deltaT:

	    //Ch L:
	    c = new TCanvas(Form("c_energyL_high_deltaT_%s",labelLR_energyBin.c_str()),Form("c_energyL_high_deltaT_%s",labelLR_energyBin.c_str()));
            histo = h1_energyL_high_deltaT[index2];
            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
            histo -> SetTitle(Form(";energy [a.u.];entries"));
            histo -> SetLineColor(kRed);
            histo -> SetLineWidth(2);
            //c -> SetLogy();
            histo -> Draw();
            histo -> Write();

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), histo->GetMean(), histo->GetRMS()));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

            c -> Print(Form("%s/Energy_highDeltaT/energyL_high_deltaT__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Energy_highDeltaT/energyL_high_deltaT__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete c;
            delete latex;


	    // Ch R:
	    c = new TCanvas(Form("c_h1_energyR_high_deltaT_%s",labelLR_energyBin.c_str()),Form("c_h1_energyR_high_deltaT_%s",labelLR_energyBin.c_str()));
            histo = h1_energyR_high_deltaT[index2];
            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
            histo -> SetTitle(Form(";energy [a.u.];entries"));
            histo -> SetLineColor(kRed);
            histo -> SetLineWidth(2);
            //c -> SetLogy();
            histo -> Draw();
            histo -> Write();

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), histo->GetMean(), histo->GetRMS()));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

            c -> Print(Form("%s/Energy_highDeltaT/energyR_high_deltaT__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Energy_highDeltaT/energyR_high_deltaT__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete c;
            delete latex;



	    // -----
	    

	    // Ch L:
	    c = new TCanvas(Form("c_h2_energy_vs_deltaT_high_L_%s",labelLR_energyBin.c_str()),Form("c_h2_energy_vs_deltaT_high_L_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            // --- Draw TH2F first ---
            h2 = h2_energy_vs_deltaT_high_L[index2];
            h2->SetTitle(";#DeltaT (Ch L) [ps];Energy (Ch L) [a.u.]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            // --- Text ---
            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Energy_highDeltaT/c_h2_energy_vs_deltaT_high_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Energy_highDeltaT/c_h2_energy_vs_deltaT_high_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;



	    // Ch R:
            c = new TCanvas(Form("c_h2_energy_vs_deltaT_high_R_%s",labelLR_energyBin.c_str()),Form("c_h2_energy_vs_deltaT_high_R_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_energy_vs_deltaT_high_R[index2];
            h2->SetTitle(";#Delta T_{R} [ps];Energy (Ch L) [a.u.]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Energy_highDeltaT/c_h2_energy_vs_deltaT_high_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Energy_highDeltaT/c_h2_energy_vs_deltaT_high_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    // -- delta T TW corr vs t1fine

	    //Ch L:
            c = new TCanvas(Form("c_deltaT_TWCorr_L_vs_t1fineL_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_L_vs_t1fineL_%s",labelLR_energyBin.c_str()));
            prof = p1_deltaT_TWCorrL_vs_t1fineL[index2];
            prof -> SetTitle(Form(";t1fine_{L};#Delta T_{L}[ps]"));
	    prof -> GetYaxis() -> SetRangeUser(CTRMeans_TW_L[index2]-600.,CTRMeans_TW_L[index2]+600.);
            prof -> Draw("");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

	    
            //float fitXMin_L = ranges["L"][index1]->at(0);
            //float fitXMax_L = ranges["L"][index1]->at(1);

            fitFunc_deltaT_TWCorrL_t1fineL[index2] = new TF1(Form("fitFunc_deltaT_TWCorrL_t1fineL_%s",labelLR_energyBin.c_str()),"pol3",300,850);
            prof -> Fit(fitFunc_deltaT_TWCorrL_t1fineL[index2],"QRS+");
            fitFunc_deltaT_TWCorrL_t1fineL[index2] -> SetLineColor(kRed);
            fitFunc_deltaT_TWCorrL_t1fineL[index2] -> SetLineWidth(2);
            fitFunc_deltaT_TWCorrL_t1fineL[index2] -> Draw("same");
	    

            c -> Print(Form("%s/Loop4_deltaT_TWCorr_vs_t1fine/deltaT_TWCorrL_vs_t1fineL__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop4_deltaT_TWCorr_vs_t1fine/deltaT_TWCorrL_vs_t1fineL__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    //Ch R:
	    c = new TCanvas(Form("c_deltaT_TWCorrR_vs_t1fineR_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorrR_vs_t1fineR_%s",labelLR_energyBin.c_str()));
            prof = p1_deltaT_TWCorrR_vs_t1fineR[index2];
            prof -> SetTitle(Form(";t1fine_{R};#Delta T_{R}[ps]"));
            prof -> GetYaxis() -> SetRangeUser(CTRMeans_TW_R[index2]-600.,CTRMeans_TW_R[index2]+600.);
            prof -> Draw("");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

	    fitFunc_deltaT_TWCorrR_t1fineR[index2] = new TF1(Form("fitFunc_deltaT_TWCorrR_t1fineR_%s",labelLR_energyBin.c_str()),"pol3",300,850);
            prof -> Fit(fitFunc_deltaT_TWCorrR_t1fineR[index2],"QRS+");
            fitFunc_deltaT_TWCorrR_t1fineR[index2] -> SetLineColor(kRed);
            fitFunc_deltaT_TWCorrR_t1fineR[index2] -> SetLineWidth(2);
            fitFunc_deltaT_TWCorrR_t1fineR[index2] -> Draw("same");

	    c -> Print(Form("%s/Loop4_deltaT_TWCorr_vs_t1fine/deltaT_TWCorrR_vs_t1fineR__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop4_deltaT_TWCorr_vs_t1fine/deltaT_TWCorrR_vs_t1fineR__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;



	    // -- Profile + TH2F deltaT TW corr vs t1Fine

            //Ch L:
            c = new TCanvas(Form("c_deltaT_TWCorrL_vs_t1fineL_scatter_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorrL_vs_t1fineL_scatter_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_TWCorrL_vs_t1fineL[index2];
            h2->SetTitle(";t1fine_{L};#Delta T_{L} TW corr [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_TWCorrL_vs_t1fineL[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_t1fine_scatter/deltaT_TW_L_vs_t1fineL__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_t1fine_scatter/deltaT_TW_L_vs_t1fineL__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;



            //Ch R:
            c = new TCanvas(Form("c_deltaT_TWCorrR_vs_t1fineR_scatter_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorrR_vs_t1fineR_scatter_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_TWCorrR_vs_t1fineR[index2];
            h2->SetTitle(";t1fine_{R};#Delta T_{R} TW corr [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_TWCorrR_vs_t1fineR[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_t1fine_scatter/deltaT_TW_R_vs_t1fineR__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_t1fine_scatter/deltaT_TW_R_vs_t1fineR__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;



	    // -- DeltaT Corr vs energy Ratio [Ch DUT]/[mean REF bar]

	    //Ch L:
            c = new TCanvas(Form("c_deltaT_TWCorr_L_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_L_vs_energyRatioREF_%s",labelLR_energyBin.c_str()));
            prof = p1_deltaT_TWCorr_L_vs_energyRatioREF[index2];
            prof -> SetTitle(Form(";E_{L}/E_{REF};#Delta T_{L} w/ TW [ps]"));
            prof -> GetXaxis() -> SetRangeUser(0.,3.5);
	    prof -> GetYaxis() -> SetRangeUser(CTRMeans_TW_L[index2] - 5.*CTRSigmas_TW_L[index2],CTRMeans_TW_L[index2] + 5.*CTRSigmas_TW_L[index2]);
            prof -> Draw("");

            latex = new TLatex(0.40,0.85,Form("#splitline{w/ TW corr - bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

	    auto range_L = GetDynamicFitRange(prof, 0.05);
            float fitXMin_EneRatio_L = range_L.first;
            float fitXMax_EneRatio_L = range_L.second;


            fitFunc_deltaT_L_TWCorr_EneRatioREF[index2] = new TF1(Form("fitFunc_deltaT_L_TWCorr_EneRatioREF_%s",labelLR_energyBin.c_str()),"pol3",fitXMin_EneRatio_L,fitXMax_EneRatio_L);
            prof -> Fit(fitFunc_deltaT_L_TWCorr_EneRatioREF[index2],"QRS+");
            fitFunc_deltaT_L_TWCorr_EneRatioREF[index2] -> SetLineColor(kRed);
            fitFunc_deltaT_L_TWCorr_EneRatioREF[index2] -> SetLineWidth(2);
            fitFunc_deltaT_L_TWCorr_EneRatioREF[index2] -> Draw("same");

            c -> Print(Form("%s/Loop4_deltaT_TWCorr_vs_energyRatioREF/p1_deltaT_TWCorr_L_vs_energyRatioREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop4_deltaT_TWCorr_vs_energyRatioREF/p1_deltaT_TWCorr_L_vs_energyRatioREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    //Ch R:
            c = new TCanvas(Form("c_deltaT_TWCorr_R_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_R_vs_energyRatioREF_%s",labelLR_energyBin.c_str()));
            prof = p1_deltaT_TWCorr_R_vs_energyRatioREF[index2];
            prof -> SetTitle(Form(";E_{R}/E_{REF};#Delta T_{R} w/ TW [ps]"));
            prof -> GetXaxis() -> SetRangeUser(0.,3.5);
	    prof -> GetYaxis() -> SetRangeUser(CTRMeans_TW_R[index2] - 5.*CTRSigmas_TW_R[index2],CTRMeans_TW_R[index2] + 5.*CTRSigmas_TW_R[index2]);
            prof -> Draw("");

            latex = new TLatex(0.40,0.85,Form("#splitline{w/ TW corr - bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

	    auto range_R = GetDynamicFitRange(prof, 0.05);
            float fitXMin_EneRatio_R = range_R.first;
            float fitXMax_EneRatio_R = range_R.second;


            fitFunc_deltaT_R_TWCorr_EneRatioREF[index2] = new TF1(Form("fitFunc_deltaT_R_TWCorr_EneRatioREF_%s",labelLR_energyBin.c_str()),"pol3",fitXMin_EneRatio_R,fitXMax_EneRatio_R);
            prof -> Fit(fitFunc_deltaT_R_TWCorr_EneRatioREF[index2],"QRS+");
            fitFunc_deltaT_R_TWCorr_EneRatioREF[index2] -> SetLineColor(kRed);
            fitFunc_deltaT_R_TWCorr_EneRatioREF[index2] -> SetLineWidth(2);
            fitFunc_deltaT_R_TWCorr_EneRatioREF[index2] -> Draw("same");

            c -> Print(Form("%s/Loop4_deltaT_TWCorr_vs_energyRatioREF/p1_deltaT_TWCorr_R_vs_energyRatioREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop4_deltaT_TWCorr_vs_energyRatioREF/p1_deltaT_TWCorr_R_vs_energyRatioREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;




	    // -- Profile + TH2F deltaT vs energy Ratio REF (with TW corr)

            //Ch L:
            c = new TCanvas(Form("c_deltaT_TWCorr_L_vs_energyRatioREF_scatter_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_L_vs_energyRatioREF_scatter_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_TWCorr_L_vs_energyRatioREF[index2];
            h2->SetTitle(";E_{L}/E_{REF} [a.u.];#Delta T_{L} [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_TWCorr_L_vs_energyRatioREF[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_energyRatioREF_scatter/h2_deltaT_TWCorr_L_vs_energyRatioREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_energyRatioREF_scatter/h2_deltaT_TWCorr_L_vs_energyRatioREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    //Ch R:
            c = new TCanvas(Form("c_deltaT_TWCorr_R_vs_energyRatioREF_scatter_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_R_vs_energyRatioREF_scatter_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_TWCorr_R_vs_energyRatioREF[index2];
            h2->SetTitle(";E_{L}/E_{REF} [a.u.];#Delta T_{R} [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_TWCorr_R_vs_energyRatioREF[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_energyRatioREF_scatter/h2_deltaT_TWCorr_R_vs_energyRatioREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_energyRatioREF_scatter/h2_deltaT_TWCorr_R_vs_energyRatioREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;



	    // -- Delta T histograms with Energy RatioREF correction

              // Ch L:
              c = new TCanvas(Form("c_h1_deltaT_L_EneRatioREFcorr_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_L_EneRatioREFcorr_%s",labelLR_energyBin.c_str()));
              histo = h1_deltaT_L_EneRatioREFcorr[index2];
              histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#Delta T_{L} w/ ERatioREF [ps];entries"));
              histo -> SetLineColor(kGreen);
              histo -> SetLineWidth(2);
              //c -> SetLogy();
              histo -> Draw();
              histo -> Write();

              fitFunc_deltaT_EneRatioREFcorr_L[index2] = new TF1(Form("fitFunc_deltaT_EneRatioREFcorr_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
              histo -> Fit(fitFunc_deltaT_EneRatioREFcorr_L[index2],"QNRS");
              histo -> Fit(fitFunc_deltaT_EneRatioREFcorr_L[index2],"QSR+","",fitFunc_deltaT_EneRatioREFcorr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREFcorr_L[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREFcorr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREFcorr_L[index2]->GetParameter(2));
              histo -> Fit(fitFunc_deltaT_EneRatioREFcorr_L[index2],"QSR+","",fitFunc_deltaT_EneRatioREFcorr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREFcorr_L[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREFcorr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREFcorr_L[index2]->GetParameter(2));

              fitFunc_deltaT_EneRatioREFcorr_L[index2] -> SetLineColor(kBlack);
              fitFunc_deltaT_EneRatioREFcorr_L[index2] -> SetLineWidth(2);
              fitFunc_deltaT_EneRatioREFcorr_L[index2] -> Draw("same");

              latex = new TLatex(0.40,0.85,Form("#splitline{w/ ERatioREF bar %02d L}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_EneRatioREFcorr_L[index2]->GetParameter(1), fitFunc_deltaT_EneRatioREFcorr_L[index2]->GetParameter(2)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop4_deltaT_EneRatioREF_corr/deltaT_EneRatioREFCorr_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop4_deltaT_EneRatioREF_corr/deltaT_EneRatioREFCorr_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete c;
              delete latex;


	      // Ch R:
              c = new TCanvas(Form("c_h1_deltaT_R_EneRatioREFcorr_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_R_EneRatioREFcorr_%s",labelLR_energyBin.c_str()));
              histo = h1_deltaT_R_EneRatioREFcorr[index2];
              histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";#Delta T_{R} w/ ERatioREF [ps];entries"));
              histo -> SetLineColor(kGreen);
              histo -> SetLineWidth(2);
              //c -> SetLogy();
              histo -> Draw();
              histo -> Write();

              fitFunc_deltaT_EneRatioREFcorr_R[index2] = new TF1(Form("fitFunc_deltaT_EneRatioREFcorr_R_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
              histo -> Fit(fitFunc_deltaT_EneRatioREFcorr_R[index2],"QNRS");
              histo -> Fit(fitFunc_deltaT_EneRatioREFcorr_R[index2],"QSR+","",fitFunc_deltaT_EneRatioREFcorr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREFcorr_R[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREFcorr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREFcorr_R[index2]->GetParameter(2));
              histo -> Fit(fitFunc_deltaT_EneRatioREFcorr_R[index2],"QSR+","",fitFunc_deltaT_EneRatioREFcorr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREFcorr_R[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREFcorr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREFcorr_R[index2]->GetParameter(2));

              fitFunc_deltaT_EneRatioREFcorr_R[index2] -> SetLineColor(kBlack);
              fitFunc_deltaT_EneRatioREFcorr_R[index2] -> SetLineWidth(2);
              fitFunc_deltaT_EneRatioREFcorr_R[index2] -> Draw("same");

              latex = new TLatex(0.40,0.85,Form("#splitline{w/ ERatioREF bar %02d R}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_EneRatioREFcorr_R[index2]->GetParameter(1), fitFunc_deltaT_EneRatioREFcorr_R[index2]->GetParameter(2)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop4_deltaT_EneRatioREF_corr/deltaT_EneRatioREFCorr_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop4_deltaT_EneRatioREF_corr/deltaT_EneRatioREFCorr_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete c;
              delete latex;



	      // -- delta T Ene ratio REF corr vs energy

              //Ch L:
              c = new TCanvas(Form("c_p1_deltaT_EneRatioREFCorr_L_vs_energy_%s",labelLR_energyBin.c_str()),Form("c_p1_deltaT_EneRatioREFCorr_L_vs_energy_%s",labelLR_energyBin.c_str()));
              prof = p1_deltaT_EneRatioREFCorr_L_vs_energy[index2];
              prof -> SetTitle(Form(";E_{L};#Delta T_{L} w/ ERatioREF [ps]"));
              prof -> GetYaxis() -> SetRangeUser(CTRMeans_ERatioREF_L[index2]-5.*CTRSigmas_ERatioREF_L[index2],CTRMeans_ERatioREF_L[index2]+5.*CTRSigmas_ERatioREF_L[index2]);
              prof -> Draw("");

              latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");

	      float fitXMin_EneratioREF_TW_L = ranges["L"][index1]->at(0);
              float fitXMax_EneratioREF_TW_L = ranges["L"][index1]->at(1);

              fitFunc_deltaT_EneRatioREF_TW_L[index2] = new TF1(Form("fitFunc_deltaT_EneRatioREF_TW_L_%s",labelLR_energyBin.c_str()),"pol3",fitXMin_EneratioREF_TW_L,fitXMax_EneratioREF_TW_L);
              prof -> Fit(fitFunc_deltaT_EneRatioREF_TW_L[index2],"QRS+");
              fitFunc_deltaT_EneRatioREF_TW_L[index2] -> SetLineColor(kRed);
              fitFunc_deltaT_EneRatioREF_TW_L[index2] -> SetLineWidth(2);
              fitFunc_deltaT_EneRatioREF_TW_L[index2] -> Draw("same");

              c -> Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_energy/deltaT_EneRatioREFCorr_L_vs_energy__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_energy/deltaT_EneRatioREFCorr_L_vs_energy__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;


	      //Ch R:
              c = new TCanvas(Form("c_p1_deltaT_EneRatioREFCorr_R_vs_energy_%s",labelLR_energyBin.c_str()),Form("c_p1_deltaT_EneRatioREFCorr_R_vs_energy_%s",labelLR_energyBin.c_str()));
              prof = p1_deltaT_EneRatioREFCorr_R_vs_energy[index2];
              prof -> SetTitle(Form(";E_{R};#Delta T_{R} w/ ERatioREF [ps]"));
	      prof -> GetYaxis() -> SetRangeUser(CTRMeans_ERatioREF_R[index2]-5.*CTRSigmas_ERatioREF_R[index2],CTRMeans_ERatioREF_R[index2]+5.*CTRSigmas_ERatioREF_R[index2]);
              prof -> Draw("");

              latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kRed);
              latex -> Draw("same");


              float fitXMin_EneratioREF_TW_R = ranges["L"][index1]->at(0);
              float fitXMax_EneratioREF_TW_R = ranges["L"][index1]->at(1);

              fitFunc_deltaT_EneRatioREF_TW_R[index2] = new TF1(Form("fitFunc_deltaT_EneRatioREF_TW_R_%s",labelLR_energyBin.c_str()),"pol3",fitXMin_EneratioREF_TW_R,fitXMax_EneratioREF_TW_R);
              prof -> Fit(fitFunc_deltaT_EneRatioREF_TW_R[index2],"QRS+");
              fitFunc_deltaT_EneRatioREF_TW_R[index2] -> SetLineColor(kRed);
              fitFunc_deltaT_EneRatioREF_TW_R[index2] -> SetLineWidth(2);
              fitFunc_deltaT_EneRatioREF_TW_R[index2] -> Draw("same");


              c -> Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_energy/deltaT_EneRatioREFCorr_R_vs_energy__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_energy/deltaT_EneRatioREFCorr_R_vs_energy__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;



	      // -- Profile + TH2F deltaT with EneRatioREF corr vs energy

            //Ch L:
            c = new TCanvas(Form("c_deltaT_EneRatioREFCorr_L_vs_energy_scatter_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREFCorr_L_vs_energy_scatter%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_EneRatioREFCorr_L_vs_energy[index2];
            h2->SetTitle(";Energy (Ch L) [a.u.];#Delta T_{L} w/ ERatioREF [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_EneRatioREFCorr_L_vs_energy[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_energy_scatter/deltaT_EneRatioREFCorr_L_vs_energy__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_energy_scatter/deltaT_EneRatioREFCorr_L_vs_energy__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    //Ch R:
            c = new TCanvas(Form("c_deltaT_EneRatioREFCorr_R_vs_energy_scatter_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREFCorr_R_vs_energy_scatter%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_EneRatioREFCorr_R_vs_energy[index2];
            h2->SetTitle(";Energy (Ch R) [a.u.];#Delta T_{R} w/ ERatioREF [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_EneRatioREFCorr_R_vs_energy[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_energy_scatter/deltaT_EneRatioREFCorr_R_vs_energy__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_energy_scatter/deltaT_EneRatioREFCorr_R_vs_energy__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    // -- Profile + TH2F deltaT with EneRatioREF corr vs EneRatioREF

            //Ch L:
            c = new TCanvas(Form("c_deltaT_EneRatioREFCorr_L_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREFCorr_L_vs_energyRatioREF%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_EneRatioREFCorr_L_vs_energyRatioREF[index2];
            h2->SetTitle(";E_{L}/<E_{REF}> [a.u.];#Delta T_{L} w/ ERatioREF [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->GetXaxis()->SetRangeUser(0., 3.5);
	    h2->Draw("colz");

            prof = p1_deltaT_EneRatioREFCorr_L_vs_energyRatioREF[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_energyRatioREF/deltaT_EneRatioREFCorr_L_vs_energyRatioREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_energyRatioREF/deltaT_EneRatioREFCorr_L_vs_energyRatioREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    //Ch R:
            c = new TCanvas(Form("c_deltaT_EneRatioREFCorr_R_vs_energyRatioREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREFCorr_R_vs_energyRatioREF%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_EneRatioREFCorr_R_vs_energyRatioREF[index2];
            h2->SetTitle(";E_{R}/<E_{REF}> [a.u.];#Delta T_{R} w/ ERatioREF [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
	    h2->GetXaxis()->SetRangeUser(0., 3.5);
            h2->Draw("colz");

            prof = p1_deltaT_EneRatioREFCorr_R_vs_energyRatioREF[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_energyRatioREF/deltaT_EneRatioREFCorr_R_vs_energyRatioREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_energyRatioREF/deltaT_EneRatioREFCorr_R_vs_energyRatioREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;
	   



	    // -- delta T EneRatioREF corr vs t1fine

            //Ch L:
            c = new TCanvas(Form("c_deltaT_EneRatioREFCorr_L_vs_t1fineL_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREFCorr_L_vs_t1fineL_%s",labelLR_energyBin.c_str()));
            prof = p1_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2];
            prof -> SetTitle(Form(";t1fine_{L};#Delta T_{L} EneRatioREF [ps]"));
            prof -> GetYaxis() -> SetRangeUser(CTRMeans_ERatioREF_L[index2]-600.,CTRMeans_ERatioREF_L[index2]+600.);
            prof -> Draw("");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

            
           // float fitXMin_L = ranges["L"][index1]->at(0);
           // float fitXMax_L = ranges["L"][index1]->at(1);

            fitFunc_deltaT_EneRCorrL_t1fineL[index2] = new TF1(Form("fitFunc_deltaT_EneRCorrL_t1fineL_%s",labelLR_energyBin.c_str()),"pol3",300,850);
            prof -> Fit(fitFunc_deltaT_EneRCorrL_t1fineL[index2],"QRS+");
            fitFunc_deltaT_EneRCorrL_t1fineL[index2] -> SetLineColor(kRed);
            fitFunc_deltaT_EneRCorrL_t1fineL[index2] -> SetLineWidth(2);
            fitFunc_deltaT_EneRCorrL_t1fineL[index2] -> Draw("same");
            

            c -> Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_t1fine/deltaT_EneRatioREFCorr_L_vs_t1fineL__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_t1fine/deltaT_EneRatioREFCorr_L_vs_t1fineL__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;



            //Ch R:
            c = new TCanvas(Form("c_deltaT_EneRatioREFCorr_R_vs_t1fineR_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREFCorr_R_vs_t1fineR_%s",labelLR_energyBin.c_str()));
            prof = p1_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2];
            prof -> SetTitle(Form(";t1fine_{R};#Delta T_{R} EneRatioREF [ps]"));
	    prof -> GetYaxis() -> SetRangeUser(CTRMeans_ERatioREF_R[index2]-600.,CTRMeans_ERatioREF_R[index2]+600.);
            prof -> Draw("");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex -> SetNDC();
            latex -> SetTextFont(42);
            latex -> SetTextSize(0.04);
            latex -> SetTextColor(kRed);
            latex -> Draw("same");

            
           // float fitXMin_EneRatio_t1fineL = ranges["L"][index1]->at(0);
           // float fitXMin_EneRatio_t1fineL = ranges["L"][index1]->at(1);

            fitFunc_deltaT_EneRCorrR_t1fineR[index2] = new TF1(Form("fitFunc_deltaT_EneRCorrR_t1fineR_%s",labelLR_energyBin.c_str()),"pol3",300,850);
            prof -> Fit(fitFunc_deltaT_EneRCorrR_t1fineR[index2],"QRS+");
            fitFunc_deltaT_EneRCorrR_t1fineR[index2] -> SetLineColor(kRed);
            fitFunc_deltaT_EneRCorrR_t1fineR[index2] -> SetLineWidth(2);
            fitFunc_deltaT_EneRCorrR_t1fineR[index2] -> Draw("same");
            

            c -> Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_t1fine/deltaT_EneRatioREFCorr_R_vs_t1fineR__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c -> Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_t1fine/deltaT_EneRatioREFCorr_R_vs_t1fineR__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


            // -- Profile + TH2F deltaT EneRatioREF corr vs t1Fine

            //Ch L:
            c = new TCanvas(Form("c_deltaT_EneRatioREFCorr_L_vs_t1fineL_scatter_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREFCorr_L_vs_t1fineL_scatter_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2];
            h2->SetTitle(";t1fine_{L};#Delta T_{L} EneRatioREF corr [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_t1fine_scatter/deltaT_EneRatioREFCorr_L_vs_t1fineL__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_t1fine_scatter/deltaT_EneRatioREFCorr_L_vs_t1fineL__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


            //Ch R:
            c = new TCanvas(Form("c_deltaT_EneRatioREFCorr_R_vs_t1fineR_scatter_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREFCorr_R_vs_t1fineR_scatter_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2];
            h2->SetTitle(";t1fine_{R};#Delta T_{R} EneRatioREF corr [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_t1fine_scatter/deltaT_EneRatioREFCorr_R_vs_t1fineR__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_t1fine_scatter/deltaT_EneRatioREFCorr_R_vs_t1fineR__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;




            // -- Profile + TH2F deltaT TW corr vs Mean Energy REF bar

            //Ch L:
            c = new TCanvas(Form("c_deltaT_TWCorr_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_TWCorr_L_vs_MeanEneREF[index2];
            h2->SetTitle(";<E_{REF}>;#Delta T_{L} TW corr [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_TWCorr_L_vs_MeanEneREF[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_MeanEneREF/deltaT_TWCorr_L_vs_MeanEneREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_MeanEneREF/deltaT_TWCorr_L_vs_MeanEneREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


            //Ch R:
            c = new TCanvas(Form("c_deltaT_TWCorr_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_TWCorr_R_vs_MeanEneREF[index2];
            h2->SetTitle(";<E_{REF}>;#Delta T_{R} TW corr [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_TWCorr_R_vs_MeanEneREF[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_MeanEneREF/deltaT_TWCorr_R_vs_MeanEneREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_TWCorr_vs_MeanEneREF/deltaT_TWCorr_R_vs_MeanEneREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


            // -- Profile + TH2F deltaT EneRatioREF corr vs Mean Energy REF bar

            //Ch L:
            c = new TCanvas(Form("c_deltaT_EneRatioREFCorr_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREFCorr_L_vs_MeanEneREF_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_EneRatioREFCorr_L_vs_MeanEneREF[index2];
            h2->SetTitle(";<E_{REF}>;#Delta T_{L} EneRatioREF corr [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_EneRatioREFCorr_L_vs_MeanEneREF[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_MeanEneREF/deltaT_EneRatioREFCorr_L_vs_MeanEneREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_MeanEneREF/deltaT_EneRatioREFCorr_L_vs_MeanEneREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    //Ch R:
            c = new TCanvas(Form("c_deltaT_EneRatioREFCorr_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREFCorr_R_vs_MeanEneREF_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_deltaT_EneRatioREFCorr_R_vs_MeanEneREF[index2];
            h2->SetTitle(";<E_{REF}>;#Delta T_{R} EneRatioREF corr [ps]");
            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_deltaT_EneRatioREFCorr_R_vs_MeanEneREF[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_MeanEneREF/deltaT_EneRatioREFCorr_R_vs_MeanEneREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_deltaT_EneRatioREFCorr_vs_MeanEneREF/deltaT_EneRatioREFCorr_R_vs_MeanEneREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;





	    //Scatter plot t1FineL vs t1FineR:

            c = new TCanvas(Form("c_h2_t1fineL_vs_t1fineR_scatter_%s",labelLR_energyBin.c_str()),Form("c_h2_t1fineL_vs_t1fineR_scatter_%s",labelLR_energyBin.c_str()));
            c->SetGridy();

            h2 = h2_t1fineL_vs_t1fineR[index2];
            h2->SetTitle(";t1Fine_{DUT} L;t1Fine_{DUT} R");
            //h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
            h2->Draw("colz");

            prof = p1_t1fineL_vs_t1fineR[index2];
            prof->SetLineColor(kRed);
            prof->SetLineWidth(3);
            prof->Draw("same");

            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
            latex->SetNDC();
            latex->SetTextFont(42);
            latex->SetTextSize(0.04);
            latex->SetTextColor(kRed);
            latex->Draw("same");

            c->Print(Form("%s/Loop4_t1fineL_vs_t1fineR_scatter/h2_t1fineL_vs_t1fineR_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
            c->Print(Form("%s/Loop4_t1fineL_vs_t1fineR_scatter/h2_t1fineL_vs_t1fineR_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
            delete latex;
            delete c;


	    }
	}

    }




  






  //------------------------
  //--- 5nd loop over events
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
          long long deltaT_L_raw = anEvent->timeL - timeMean_ext;
          long long deltaT_R_raw = anEvent->timeR - timeMean_ext;

	  float energyMeanREF = 0.5 *(anEvent->energyL_ext + anEvent->energyR_ext);


          float deltaT_L_raw_low = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_L_raw_hig = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);

          float deltaT_R_raw_low = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_R_raw_hig = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);

	  // -- TW corr only fitFunc(energy):
          float TW_corr_L = fitFunc_deltaT_EneLCorr_L[index2]->Eval(anEvent->energyL);
          float TW_corr_R = fitFunc_deltaT_EneRCorr_R[index2]->Eval(anEvent->energyR);

	  // -- EneRatioREF correction (from delta T raw):
          float EneRatio_corr_L = fitFunc_deltaT_L_raw_ERatioREF[index2]->Eval(anEvent->energyL/energyMeanREF);
          float EneRatio_corr_R = fitFunc_deltaT_R_raw_ERatioREF[index2]->Eval(anEvent->energyR/energyMeanREF);

	  // -- TW + EneRatioREF corr
	  float TW_corr_EneRatioREF_L = fitFunc_deltaT_L_TWCorr_EneRatioREF[index2]->Eval(anEvent->energyL/energyMeanREF);
	  float TW_corr_EneRatioREF_R = fitFunc_deltaT_R_TWCorr_EneRatioREF[index2]->Eval(anEvent->energyR/energyMeanREF);

          // -- EneRatioREF + TW corr
	  float EneRatioREF_TW_Corr_L = fitFunc_deltaT_EneRatioREF_TW_L[index2]->Eval(anEvent->energyL);
	  float EneRatioREF_TW_Corr_R = fitFunc_deltaT_EneRatioREF_TW_R[index2]->Eval(anEvent->energyR);

	  // -- TW + t1fine corr
	  //float TW_corr_t1fine_L = fitFunc_deltaT_TWCorrL_t1fineL[index2]->Eval(anEvent->t1fineL);
	  //float TW_corr_t1fine_R = fitFunc_deltaT_TWCorrR_t1fineR[index2]->Eval(anEvent->t1fineR);

	  float TW_corr_t1fine_L = GetProfileBinCorrection(p1_deltaT_TWCorrL_vs_t1fineL[index2],anEvent->t1fineL); //--> corr bin x bin
          float TW_corr_t1fine_R = GetProfileBinCorrection(p1_deltaT_TWCorrR_vs_t1fineR[index2],anEvent->t1fineR);

	  // -- EneRatio REF + t1fine corr
          //float EneR_corr_t1fine_L = fitFunc_deltaT_EneRCorrL_t1fineL[index2]->Eval(anEvent->t1fineL);
          //float EneR_corr_t1fine_R = fitFunc_deltaT_EneRCorrR_t1fineR[index2]->Eval(anEvent->t1fineR);

	  float EneR_corr_t1fine_L = GetProfileBinCorrection(p1_deltaT_EneRatioREFCorr_L_vs_t1fineL[index2],anEvent->t1fineL); //--> corr bin x bin
	  float EneR_corr_t1fine_R = GetProfileBinCorrection(p1_deltaT_EneRatioREFCorr_R_vs_t1fineR[index2],anEvent->t1fineR);


	  // -- Mean t1fine REF
	  float Mean_t1fineREF = 0.5*(anEvent->t1fineL_ext + anEvent->t1fineR_ext);

          std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));

          // -- delta T histograms
          if( h1_deltaT_TWCorr_EneRatioCorr_L[index2] == NULL )
	  {

	    // -- TW + EneRatio Corr	  
            h1_deltaT_TWCorr_EneRatioCorr_L[index2] = new TH1F(Form("h1_TWCorr_EneRatioCorr_L_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	    h1_deltaT_TWCorr_EneRatioCorr_R[index2] = new TH1F(Form("h1_TWCorr_EneRatioCorr_R_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);

            h2_deltaT_TWCorr_EneRatioCorr_vs_energy_L[index2] = new TH2F(Form("h2_deltaT_TWCorr_EneRatioCorr_vs_energy_L_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800., 2000, -12000., 12000.);
	    h2_deltaT_TWCorr_EneRatioCorr_vs_energy_R[index2] = new TH2F(Form("h2_deltaT_TWCorr_EneRatioCorr_vs_energy_R_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800., 2000, -12000., 12000.);
            p1_deltaT_TWCorr_EneRatioCorr_vs_energy_L[index2] = new TProfile(Form("p1_deltaT_TWCorr_EneRatioCorr_vs_energy_L_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800.);
	    p1_deltaT_TWCorr_EneRatioCorr_vs_energy_R[index2] = new TProfile(Form("p1_deltaT_TWCorr_EneRatioCorr_vs_energy_R_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800.);

            p1_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L[index2] = new TProfile(Form("p1_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L_%s",labelLR_energyBin.c_str()),"",100,0.,5.);
            p1_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R[index2] = new TProfile(Form("p1_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R_%s",labelLR_energyBin.c_str()),"",100,0.,5.);
            h2_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L[index2] = new TH2F(Form("h2_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);
	    h2_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R[index2] = new TH2F(Form("h2_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);
	  
	  

	  // -- EneRatio + TW Corr  
	  h1_deltaT_EneRatioREF_TW_Corr_L[index2] = new TH1F(Form("h1_deltaT_EneRatioREF_TW_Corr_L_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	  h1_deltaT_EneRatioREF_TW_Corr_R[index2] = new TH1F(Form("h1_deltaT_EneRatioREF_TW_Corr_R_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	  h2_deltaT_EneRatioREF_TW_Corr_vs_energy_L[index2] = new TH2F(Form("h2_deltaT_EneRatioREF_TW_vs_energy_L_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800., 2000, -12000., 12000.);
	  h2_deltaT_EneRatioREF_TW_Corr_vs_energy_R[index2] = new TH2F(Form("h2_deltaT_EneRatioREF_TW_vs_energy_R_%s",labelLR_energyBin.c_str()),"",50,ranges["R"][index1]->at(0),800., 2000, -12000., 12000.);
	  p1_deltaT_EneRatioREF_TW_Corr_vs_energy_L[index2] = new TProfile(Form("p1_deltaT_EneRatioREF_TW_vs_energy_L_%s",labelLR_energyBin.c_str()),"",50,ranges["L"][index1]->at(0),800.);
	  p1_deltaT_EneRatioREF_TW_Corr_vs_energy_R[index2] = new TProfile(Form("p1_deltaT_EneRatioREF_TW_vs_energy_R_%s",labelLR_energyBin.c_str()),"",50,ranges["R"][index1]->at(0),800.);

	  h2_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_L[index2] = new TH2F(Form("h2_deltaT_EneRatioREF_TW_vs_energyRatioREF_L_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);
	  h2_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_R[index2] = new TH2F(Form("h2_deltaT_EneRatioREF_TW_vs_energyRatioREF_R_%s",labelLR_energyBin.c_str()),"",100,0.,5., 2000, -12000., 12000.);
	  
	  p1_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_L[index2] = new TProfile(Form("p1_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_L_%s",labelLR_energyBin.c_str()),"",100,0.,5.);
	  p1_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_R[index2] = new TProfile(Form("p1_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_R_%s",labelLR_energyBin.c_str()),"",100,0.,5.);
	  
	  
	  // -- TW + t1fine corr
	  h1_deltaT_TW_t1fine_Corr_L[index2] = new TH1F(Form("h1_deltaT_TW_t1fine_Corr_L_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
          h1_deltaT_TW_t1fine_Corr_R[index2] = new TH1F(Form("h1_deltaT_TW_t1fine_Corr_R_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	  
	  // -- EneRatio REF + t1fine corr
	  h1_deltaT_EneRatioREF_t1fine_Corr_L[index2] = new TH1F(Form("h1_deltaT_EneRatioREF_t1fine_Corr_L_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
          h1_deltaT_EneRatioREF_t1fine_Corr_R[index2] = new TH1F(Form("h1_deltaT_EneRatioREF_t1fine_Corr_R_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	  

	  // -- Mean t1fine REF
          p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2] = new TProfile(Form("p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
          h2_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2] = new TH2F(Form("h2_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0.,1000., 2000, -12000., 12000.);
	  p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2] = new TProfile(Form("p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
          h2_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2] = new TH2F(Form("h2_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0.,1000., 2000, -12000., 12000.);
	  
	  p1_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2] = new TProfile(Form("p1_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
          h2_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2] = new TH2F(Form("h2_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0.,1000.,2000,-12000.,12000.);
          p1_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2] = new TProfile(Form("p1_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0,1000.);
          h2_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2] = new TH2F(Form("h2_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),"",50,0.,1000.,2000,-12000.,12000.);
	  }


	  if((deltaT_L_raw > deltaT_L_raw_low ) && (deltaT_L_raw < deltaT_L_raw_hig)){
          h1_deltaT_TWCorr_EneRatioCorr_L[index2] -> Fill(deltaT_L_raw - TW_corr_L - TW_corr_EneRatioREF_L);
          
	  h2_deltaT_TWCorr_EneRatioCorr_vs_energy_L[index2] -> Fill( anEvent->energyL, deltaT_L_raw - TW_corr_L - TW_corr_EneRatioREF_L);
          p1_deltaT_TWCorr_EneRatioCorr_vs_energy_L[index2] -> Fill( anEvent->energyL, deltaT_L_raw - TW_corr_L - TW_corr_EneRatioREF_L);
          
	  p1_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L[index2] -> Fill ( anEvent->energyL/energyMeanREF, deltaT_L_raw - TW_corr_L - TW_corr_EneRatioREF_L);
          h2_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L[index2] -> Fill( anEvent->energyL/energyMeanREF, deltaT_L_raw - TW_corr_L - TW_corr_EneRatioREF_L);
          

	  // -- EneRatio + TW Corr
	  h1_deltaT_EneRatioREF_TW_Corr_L[index2] -> Fill(deltaT_L_raw - EneRatio_corr_L - EneRatioREF_TW_Corr_L);
          h2_deltaT_EneRatioREF_TW_Corr_vs_energy_L[index2] -> Fill( anEvent->energyL, deltaT_L_raw - EneRatio_corr_L - EneRatioREF_TW_Corr_L);
	  p1_deltaT_EneRatioREF_TW_Corr_vs_energy_L[index2] -> Fill( anEvent->energyL, deltaT_L_raw - EneRatio_corr_L - EneRatioREF_TW_Corr_L);

	  h2_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_L[index2] -> Fill( anEvent->energyL/energyMeanREF, deltaT_L_raw - EneRatio_corr_L - EneRatioREF_TW_Corr_L);
		 
	  p1_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_L[index2] -> Fill( anEvent->energyL/energyMeanREF, deltaT_L_raw - EneRatio_corr_L - EneRatioREF_TW_Corr_L); 
	  
	  
	  // -- TW + t1fine Corr
	  h1_deltaT_TW_t1fine_Corr_L[index2] -> Fill(deltaT_L_raw - TW_corr_L - TW_corr_t1fine_L);
	  
	  // -- EneRatioREF + t1fine Corr
          h1_deltaT_EneRatioREF_t1fine_Corr_L[index2] -> Fill(deltaT_L_raw - EneRatio_corr_L - EneR_corr_t1fine_L);


	  // -- Mean t1fine REF
	  p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_L_raw - EneRatio_corr_L - EneR_corr_t1fine_L);
	  h2_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_L_raw - EneRatio_corr_L - EneR_corr_t1fine_L);
	  p1_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_L_raw - TW_corr_L - TW_corr_t1fine_L);
	  h2_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_L_raw - TW_corr_L - TW_corr_t1fine_L);
	  }


	  if((deltaT_R_raw > deltaT_R_raw_low ) && (deltaT_R_raw < deltaT_R_raw_hig)){
          h1_deltaT_TWCorr_EneRatioCorr_R[index2] -> Fill(deltaT_R_raw - TW_corr_R - TW_corr_EneRatioREF_R);
          
	  h2_deltaT_TWCorr_EneRatioCorr_vs_energy_R[index2] -> Fill( anEvent->energyR, deltaT_R_raw - TW_corr_R - TW_corr_EneRatioREF_R);
          p1_deltaT_TWCorr_EneRatioCorr_vs_energy_R[index2] -> Fill( anEvent->energyR, deltaT_R_raw - TW_corr_R - TW_corr_EneRatioREF_R);
          
	  p1_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R[index2] -> Fill ( anEvent->energyR/energyMeanREF, deltaT_R_raw - TW_corr_R - TW_corr_EneRatioREF_R);
          h2_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R[index2] -> Fill( anEvent->energyR/energyMeanREF, deltaT_R_raw - TW_corr_R - TW_corr_EneRatioREF_R);


	  // -- EneRatio + TW Corr
	  h1_deltaT_EneRatioREF_TW_Corr_R[index2] -> Fill(deltaT_R_raw - EneRatio_corr_R - EneRatioREF_TW_Corr_R);
          h2_deltaT_EneRatioREF_TW_Corr_vs_energy_R[index2] -> Fill( anEvent->energyR, deltaT_R_raw - EneRatio_corr_R - EneRatioREF_TW_Corr_R);
          p1_deltaT_EneRatioREF_TW_Corr_vs_energy_R[index2] -> Fill( anEvent->energyR, deltaT_R_raw - EneRatio_corr_R - EneRatioREF_TW_Corr_R);


          h2_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_R[index2] -> Fill( anEvent->energyR/energyMeanREF, deltaT_R_raw - EneRatio_corr_R - EneRatioREF_TW_Corr_R);

          p1_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_R[index2] -> Fill( anEvent->energyR/energyMeanREF, deltaT_R_raw - EneRatio_corr_R - EneRatioREF_TW_Corr_R);
	  
	  
          // -- TW + t1fine Corr
          h1_deltaT_TW_t1fine_Corr_R[index2] -> Fill(deltaT_R_raw - TW_corr_R -TW_corr_t1fine_R);

          // -- EneRatioREF + t1fine cCorr
          h1_deltaT_EneRatioREF_t1fine_Corr_R[index2] -> Fill(deltaT_R_raw - EneRatio_corr_R - EneR_corr_t1fine_R);

	  // -- Mean t1fine REF
          p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_R_raw - EneRatio_corr_R - EneR_corr_t1fine_R);
          h2_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_R_raw - EneRatio_corr_R - EneR_corr_t1fine_R);
          p1_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_R_raw - TW_corr_R - TW_corr_t1fine_R);
          h2_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2] -> Fill (Mean_t1fineREF ,deltaT_R_raw - TW_corr_R - TW_corr_t1fine_R);
	  }



        }
    }





    //------------------
    //--- draw plots loop 5

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

    std::map<double,float> CTRMeans_TW_EneRatio_L;
    std::map<double,float> CTRSigmas_TW_EneRatio_L;

    for(auto mapIt : h1_deltaT_TWCorr_EneRatioCorr_L)
    {
        double index = mapIt.first;

      FindSmallestInterval(vals,h1_deltaT_TWCorr_EneRatioCorr_L[index],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans_TW_EneRatio_L[index] = mean;
      CTRSigmas_TW_EneRatio_L[index] = effSigma;
    }


    std::map<double,float> CTRMeans_TW_EneRatio_R;
    std::map<double,float> CTRSigmas_TW_EneRatio_R;

    for(auto mapIt : h1_deltaT_TWCorr_EneRatioCorr_R)
    {
        double index = mapIt.first;

      FindSmallestInterval(vals,h1_deltaT_TWCorr_EneRatioCorr_R[index],0.68);
      float mean = vals[0];
      float min = vals[4];
      float max = vals[5];
      float delta = max-min;
      float sigma = 0.5*delta;
      float effSigma = sigma;
      CTRMeans_TW_EneRatio_R[index] = mean;
      CTRSigmas_TW_EneRatio_R[index] = effSigma;
    }




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
            
	   


			    if(!h1_deltaT_TWCorr_EneRatioCorr_L[index2]) continue;
			    std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));


			    // -- Delta T TW + EneRatio REF corr
	    
			    // Ch L:
			    c = new TCanvas(Form("c_h1_deltaT_TWCorr_EneRatioCorr_L_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_TWCorr_EneRatioCorr_L_%s",labelLR_energyBin.c_str()));
			    histo = h1_deltaT_TWCorr_EneRatioCorr_L[index2];
			    histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),histo->GetMean()+3.*histo->GetRMS());
			    histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
			    histo -> SetTitle(Form(";#Delta T_{L} TW+EneRatio corr [ps];entries"));
			    histo -> SetLineColor(kBrown);
			    histo -> SetLineWidth(2);
			    //c -> SetLogy();
			    histo -> Draw();
			    histo -> Write();

			    fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2] = new TF1(Form("fitFunc_deltaT_TWCorr_EneRatioCorr_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
			    histo -> Fit(fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2],"QNRS");
			    histo -> Fit(fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2],"QSR+","",fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2]->GetParameter(2),fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2]->GetParameter(2));
			    histo -> Fit(fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2],"QSR+","",fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2]->GetParameter(2),fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2]->GetParameter(2));

			    fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2] -> SetLineColor(kBlack);
			    fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2] -> SetLineWidth(2);
			    fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2] -> Draw("same");

			    latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2]->GetParameter(1), fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2]->GetParameter(2)));
			    latex -> SetNDC();
			    latex -> SetTextFont(42);
			    latex -> SetTextSize(0.04);
			    latex -> SetTextColor(kRed);
			    latex -> Draw("same");

			    c -> Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr/deltaT_TW_EneRatio_Corr_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
			    c -> Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr/deltaT_TW_EneRatio_Corr_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
			    delete c;
			    delete latex;


			    // Ch R:
                            c = new TCanvas(Form("c_h1_deltaT_TWCorr_EneRatioCorr_R_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_TWCorr_EneRatioCorr_R_%s",labelLR_energyBin.c_str()));
                            histo = h1_deltaT_TWCorr_EneRatioCorr_R[index2];
                            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),histo->GetMean()+3.*histo->GetRMS());
                            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
                            histo -> SetTitle(Form(";#Delta T_{R} TW+EneRatio corr [ps];entries"));
                            histo -> SetLineColor(kBrown);
                            histo -> SetLineWidth(2);
                            //c -> SetLogy();
                            histo -> Draw();
                            histo -> Write();

                            fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2] = new TF1(Form("fitFunc_deltaT_TWCorr_EneRatioCorr_R_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
                            histo -> Fit(fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2],"QNRS");
                            histo -> Fit(fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2],"QSR+","",fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2]->GetParameter(2),fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2]->GetParameter(2));
                            histo -> Fit(fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2],"QSR+","",fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2]->GetParameter(2),fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2]->GetParameter(2));

                            fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2] -> SetLineColor(kBlack);
                            fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2] -> Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2]->GetParameter(1), fitFunc_deltaT_TWCorr_EneRatioCorr_R[index2]->GetParameter(2)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");

                            c -> Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr/deltaT_TW_EneRatio_Corr_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr/deltaT_TW_EneRatio_Corr_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete c;
                            delete latex;




	    
			    // -- profile + TH2F DeltaT vs energy (with TW + Ene ratio corr)

            
			    //Ch L:
			    c = new TCanvas(Form("c_deltaT_TWCorr_EneRatioCorr_vs_energy_L_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_EneRatioCorr_vs_energy_L_%s",labelLR_energyBin.c_str()));
			    c->SetGridy();
            
                            h2 = h2_deltaT_TWCorr_EneRatioCorr_vs_energy_L[index2];
                            h2->SetTitle(";Energy (Ch L) [a.u.];#Delta T_{L} (TW + EneRatio Corr) [ps]");
                            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
			    //h2->GetYaxis()->SetRangeUser(CTRMeans_TW_EneRatio_L[index2] - 5*CTRSigmas_TW_EneRatio_L[index], CTRMeans_TW_EneRatio_L[index] + 5.*CTRSigmas_TW_EneRatio_L[index]);
                            h2->Draw("colz");

			    prof = p1_deltaT_TWCorr_EneRatioCorr_vs_energy_L[index2];
			    prof->SetLineColor(kRed);
			    prof->SetLineWidth(3);
			    prof->Draw("same");

			    latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
			    latex->SetNDC();
			    latex->SetTextFont(42);
			    latex->SetTextSize(0.04);
			    latex->SetTextColor(kRed);
			    latex->Draw("same");

			    c->Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr_vs_Energy/c_deltaT_TW_EneRatio_Corr_vs_energy_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
			    c->Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr_vs_Energy/c_deltaT_TW_EneRatio_Corr_vs_energy_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
		            delete latex;
		            delete c;


			    //Ch R:
                            c = new TCanvas(Form("c_deltaT_TWCorr_EneRatioCorr_vs_energy_R_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_EneRatioCorr_vs_energy_R_%s",labelLR_energyBin.c_str()));
                            c->SetGridy();

                            h2 = h2_deltaT_TWCorr_EneRatioCorr_vs_energy_R[index2];
                            h2->SetTitle(";Energy (Ch R) [a.u.];#Delta T_{R} (TW + EneRatio Corr) [ps]");
                            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
			    //h2->GetYaxis()->SetRangeUser(CTRMeans_TW_EneRatio_R[index2] - 5*CTRSigmas_TW_EneRatio_R[index], CTRMeans_TW_EneRatio_R[index] + 5.*CTRSigmas_TW_EneRatio_R[index]);
			    h2->Draw("colz");

                            prof = p1_deltaT_TWCorr_EneRatioCorr_vs_energy_R[index2];
                            prof->SetLineColor(kRed);
                            prof->SetLineWidth(3);
                            prof->Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex->SetNDC();
                            latex->SetTextFont(42);
                            latex->SetTextSize(0.04);
                            latex->SetTextColor(kRed);
                            latex->Draw("same");

                            c->Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr_vs_Energy/c_deltaT_TW_EneRatio_Corr_vs_energy_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c->Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr_vs_Energy/c_deltaT_TW_EneRatio_Corr_vs_energy_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;



			    // Profile + TH2F Delta T vs Energy Ratio (with TW + Ene ratio corr)

			    //Ch L:
			    c = new TCanvas(Form("c_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L_%s",labelLR_energyBin.c_str()));

			    prof = p1_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L[index2];
			    prof -> SetTitle(Form(";E_{L}/E_{REF};#Delta T_{L} TW + EneRatio Corr [ps]"));
			    prof -> GetYaxis() -> SetRangeUser(CTRMeans_TW_EneRatio_L[index2]-5.*CTRSigmas_TW_EneRatio_L[index2],CTRMeans_TW_EneRatio_L[index2]+5.*CTRSigmas_TW_EneRatio_L[index2]);
			    prof -> Draw("");
			    latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
			    latex -> SetNDC();
			    latex -> SetTextFont(42);
			    latex -> SetTextSize(0.04);
			    latex -> SetTextColor(kRed);
			    latex -> Draw("same");

			    /*
			    float fitXMin = ranges["L"][index1]->at(0);
			    float fitXMax = ranges["L"][index1]->at(1);

			    fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2] = new TF1(Form("fitFunc_deltaT_TWCorr_EneRatioCorr_L_%s",labelLR_energyBin.c_str()),"pol3",fitXMin,fitXMax);
			    prof -> Fit(fitFunc_deltaT_EneLCorr_L[index2],"QRS+");
			    fitFunc_deltaT_EneLCorr_L[index2] -> SetLineColor(kRed);
			    fitFunc_deltaT_EneLCorr_L[index2] -> SetLineWidth(2);
			    fitFunc_deltaT_EneLCorr_L[index2] -> Draw("same");

			    */
			    c -> Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr_vs_EneRatio/c_deltaT_TW_EneRatio_Corr_vs_EneRatio_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
			    c -> Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr_vs_EneRatio/c_deltaT_TW_EneRatio_Corr_vs_EneRatio_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
			    delete latex;
			    delete c;


			    //Ch R:
                            c = new TCanvas(Form("c_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R_%s",labelLR_energyBin.c_str()));

                            prof = p1_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R[index2];
                            prof -> SetTitle(Form(";E_{R}/E_{REF};#Delta T_{R} TW + EneRatio Corr [ps]"));
			    prof -> GetYaxis() -> SetRangeUser(CTRMeans_TW_EneRatio_R[index2]-5.*CTRSigmas_TW_EneRatio_R[index2],CTRMeans_TW_EneRatio_R[index2]+5.*CTRSigmas_TW_EneRatio_R[index2]);
                            prof -> Draw("");
                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");


                            /*
                            float fitXMin = ranges["L"][index1]->at(0);
                            float fitXMax = ranges["L"][index1]->at(1);

                            fitFunc_deltaT_TWCorr_EneRatioCorr_L[index2] = new TF1(Form("fitFunc_deltaT_TWCorr_EneRatioCorr_L_%s",labelLR_energyBin.c_str()),"pol3",fitXMin,fitXMax);
                            prof -> Fit(fitFunc_deltaT_EneLCorr_L[index2],"QRS+");
                            fitFunc_deltaT_EneLCorr_L[index2] -> SetLineColor(kRed);
                            fitFunc_deltaT_EneLCorr_L[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_EneLCorr_L[index2] -> Draw("same");

                            */
                            c -> Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr_vs_EneRatio/c_deltaT_TW_EneRatio_Corr_vs_EneRatio_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr_vs_EneRatio/c_deltaT_TW_EneRatio_Corr_vs_EneRatio_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;




			    // -- profile + TH2F DeltaT vs energy (with TW + Ene ratio corr)


                            //Ch L:
                            c = new TCanvas(Form("c_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L_%s",labelLR_energyBin.c_str()));
                            c->SetGridy();

                            h2 = h2_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L[index2];
                            h2->SetTitle(";E_{L}/E_{REF} [a.u.];#Delta T_{L} (TW + EneRatio Corr) [ps]");
                            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
                            h2->GetYaxis()->SetRangeUser(CTRMeans_TW_EneRatio_L[index2] - 5*CTRSigmas_TW_EneRatio_L[index2], CTRMeans_TW_EneRatio_L[index2] + 5.*CTRSigmas_TW_EneRatio_L[index2]);
                            h2->Draw("colz");

                            prof = p1_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L[index2];
                            prof->SetLineColor(kRed);
                            prof->SetLineWidth(3);
                            prof->Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex->SetNDC();
                            latex->SetTextFont(42);
                            latex->SetTextSize(0.04);
                            latex->SetTextColor(kRed);
                            latex->Draw("same");

                            c->Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr_vs_EneRatio_scatter/c_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c->Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr_vs_EneRatio_scatter/c_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;


			    //Ch R:
                            c = new TCanvas(Form("c_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R_%s",labelLR_energyBin.c_str()));
                            c->SetGridy();

                            h2 = h2_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R[index2];
                            h2->SetTitle(";E_{R}/E_{REF} [a.u.];#Delta T_{R} (TW + EneRatio Corr) [ps]");
                            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
                            h2->GetYaxis()->SetRangeUser(CTRMeans_TW_EneRatio_R[index2] - 5*CTRSigmas_TW_EneRatio_R[index2], CTRMeans_TW_EneRatio_R[index2] + 5.*CTRSigmas_TW_EneRatio_R[index2]);
                            h2->Draw("colz");

                            prof = p1_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R[index2];
                            prof->SetLineColor(kRed);
                            prof->SetLineWidth(3);
                            prof->Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex->SetNDC();
                            latex->SetTextFont(42);
                            latex->SetTextSize(0.04);
                            latex->SetTextColor(kRed);
                            latex->Draw("same");

                            c->Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr_vs_EneRatio_scatter/c_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c->Print(Form("%s/Loop5_deltaT_TW_EneRatio_Corr_vs_EneRatio_scatter/c_deltaT_TWCorr_EneRatioCorr_vs_energyRatioREF_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;





			    // -- Delta T EnergyRatioREF + TW corr vs energy

                            // Ch L:
                            c = new TCanvas(Form("c_h1_deltaT_EneRatioREF_TW_Corr_L_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_EneRatioREF_TW_Corr_L_%s",labelLR_energyBin.c_str()));
                            histo = h1_deltaT_EneRatioREF_TW_Corr_L[index2];
                            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),histo->GetMean()+3.*histo->GetRMS());
                            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
                            histo -> SetTitle(Form(";#Delta T_{L} EneRatio+TW corr [ps];entries"));
                            histo -> SetLineColor(kOrange);
                            histo -> SetLineWidth(2);
                            //c -> SetLogy();
                            histo -> Draw();
                            histo -> Write();

                            fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2] = new TF1(Form("fitFunc_deltaT_EneRatioREF_TW_Corr_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2],"QNRS");
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2],"QSR+","",fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2]->GetParameter(2));
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2],"QSR+","",fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2]->GetParameter(2));

                            fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2] -> SetLineColor(kBlack);
                            fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2] -> Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2]->GetParameter(1), fitFunc_deltaT_EneRatioREF_TW_Corr_L[index2]->GetParameter(2)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
			    latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");

                            c -> Print(Form("%s/Loop5_deltaT_EneRatio_TW_Corr/deltaT_EneRatio_TW_Corr_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop5_deltaT_EneRatio_TW_Corr/deltaT_EneRatio_TW_Corr_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete c;
                            delete latex;



                            // Ch R:
                            c = new TCanvas(Form("c_h1_deltaT_EneRatioREF_TW_Corr_R_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_EneRatioREF_TW_Corr_R_%s",labelLR_energyBin.c_str()));
                            histo = h1_deltaT_EneRatioREF_TW_Corr_R[index2];
                            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),histo->GetMean()+3.*histo->GetRMS());
                            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
                            histo -> SetTitle(Form(";#Delta T_{R} EneRatio+TW corr [ps];entries"));
                            histo -> SetLineColor(kOrange);
                            histo -> SetLineWidth(2);
                            //c -> SetLogy();
                            histo -> Draw();
                            histo -> Write();

                            fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2] = new TF1(Form("fitFunc_deltaT_EneRatioREF_TW_Corr_R_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2],"QNRS");
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2],"QSR+","",fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2]->GetParameter(2));
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2],"QSR+","",fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2]->GetParameter(2));

                            fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2] -> SetLineColor(kBlack);
                            fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2] -> Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2]->GetParameter(1), fitFunc_deltaT_EneRatioREF_TW_Corr_R[index2]->GetParameter(2)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");


                            c -> Print(Form("%s/Loop5_deltaT_EneRatio_TW_Corr/deltaT_EneRatio_TW_Corr_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop5_deltaT_EneRatio_TW_Corr/deltaT_EneRatio_TW_Corr_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete c;
                            delete latex;




                            // -- profile + TH2F DeltaT vs energy (with Ene ratio + TW corr)


                            //Ch L:
                            c = new TCanvas(Form("c_deltaT_EneRatioREF_TW_Corr_vs_energy_L_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREF_TW_Corr_vs_energy_L_%s",labelLR_energyBin.c_str()));
                            c->SetGridy();

                            h2 = h2_deltaT_EneRatioREF_TW_Corr_vs_energy_L[index2];
                            h2->SetTitle(";E_{L} [a.u.];#Delta T_{L} (EneRatio+TW Corr) [ps]");
                            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
                            //h2->GetYaxis()->SetRangeUser(CTRMeans_TW_EneRatio_L[index2] - 5*CTRSigmas_TW_EneRatio_L[index2], CTRMeans_TW_EneRatio_L[index2] + 5.*CTRSigmas_TW_EneRatio_L[index2]);
                            h2->Draw("colz");

                            prof = p1_deltaT_EneRatioREF_TW_Corr_vs_energy_L[index2];
                            prof->SetLineColor(kRed);
                            prof->SetLineWidth(3);
                            prof->Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex->SetNDC();
                            latex->SetTextFont(42);
                            latex->SetTextSize(0.04);
                            latex->SetTextColor(kRed);
                            latex->Draw("same");

                            c->Print(Form("%s/Loop5_deltaT_EneRatio_TW_Corr_vs_Energy_scatter/c_deltaT_EneRatioREF_TW_Corr_vs_energy_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c->Print(Form("%s/Loop5_deltaT_EneRatio_TW_Corr_vs_Energy_scatter/c_deltaT_EneRatioREF_TW_Corr_vs_energy_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;




                            //Ch R:
                            c = new TCanvas(Form("c_deltaT_EneRatioREF_TW_Corr_vs_energy_R_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREF_TW_Corr_vs_energy_R_%s",labelLR_energyBin.c_str()));
                            c->SetGridy();

                            h2 = h2_deltaT_EneRatioREF_TW_Corr_vs_energy_R[index2];
                            h2->SetTitle(";E_{R} [a.u.];#Delta T_{R} (EneRatio+TW Corr) [ps]");
                            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
                            //h2->GetYaxis()->SetRangeUser(CTRMeans_TW_EneRatio_L[index2] - 5*CTRSigmas_TW_EneRatio_L[index2], CTRMeans_TW_EneRatio_L[index2] + 5.*CTRSigmas_TW_EneRatio_L[index2]);
                            h2->Draw("colz");

                            prof = p1_deltaT_EneRatioREF_TW_Corr_vs_energy_R[index2];
                            prof->SetLineColor(kRed);
                            prof->SetLineWidth(3);
                            prof->Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex->SetNDC();
                            latex->SetTextFont(42);
                            latex->SetTextSize(0.04);
                            latex->SetTextColor(kRed);
                            latex->Draw("same");

                            c->Print(Form("%s/Loop5_deltaT_EneRatio_TW_Corr_vs_Energy_scatter/c_deltaT_EneRatioREF_TW_Corr_vs_energy_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c->Print(Form("%s/Loop5_deltaT_EneRatio_TW_Corr_vs_Energy_scatter/c_deltaT_EneRatioREF_TW_Corr_vs_energy_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;





                            // -- profile + TH2F DeltaT vs energyRatioREF (with Ene ratio + TW corr)


                            //Ch L:
                            c = new TCanvas(Form("c_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_L_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_L_%s",labelLR_energyBin.c_str()));
                            c->SetGridy();

                            h2 = h2_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_L[index2];
                            h2->SetTitle(";E_{L}/E_{REF} [a.u.];#Delta T_{L} (EneRatio+TW Corr) [ps]");
                            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
                            h2->GetYaxis()->SetRangeUser(CTRMeans_TW_EneRatio_L[index2] - 5*CTRSigmas_TW_EneRatio_L[index2], CTRMeans_TW_EneRatio_L[index2] + 5.*CTRSigmas_TW_EneRatio_L[index2]);
                            h2->Draw("colz");

                            prof = p1_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_L[index2];
                            prof->SetLineColor(kRed);
                            prof->SetLineWidth(3);
                            prof->Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex->SetNDC();
                            latex->SetTextFont(42);
                            latex->SetTextSize(0.04);
                            latex->SetTextColor(kRed);
                            latex->Draw("same");

                            c->Print(Form("%s/Loop5_deltaT_EneRatio_TW_Corr_vs_EneRatioREF_scatter/c_deltaT_EneRatioREF_TW_Corr_vs_EneRatioREF_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c->Print(Form("%s/Loop5_deltaT_EneRatio_TW_Corr_vs_EneRatioREF_scatter/c_deltaT_EneRatioREF_TW_Corr_vs_EneRatioREF_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;


                            //Ch R:
                            c = new TCanvas(Form("c_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_R_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_R_%s",labelLR_energyBin.c_str()));
                            c->SetGridy();

                            h2 = h2_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_R[index2];
                            h2->SetTitle(";E_{R}/E_{REF} [a.u.];#Delta T_{R} (EneRatio+TW Corr) [ps]");
                            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
                            h2->GetYaxis()->SetRangeUser(CTRMeans_TW_EneRatio_R[index2] - 5*CTRSigmas_TW_EneRatio_R[index2], CTRMeans_TW_EneRatio_R[index2] + 5.*CTRSigmas_TW_EneRatio_R[index2]);
                            h2->Draw("colz");

                            prof = p1_deltaT_EneRatioREF_TW_Corr_vs_energyRatioREF_R[index2];
                            prof->SetLineColor(kRed);
                            prof->SetLineWidth(3);
                            prof->Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex->SetNDC();
                            latex->SetTextFont(42);
                            latex->SetTextSize(0.04);
                            latex->SetTextColor(kRed);
                            latex->Draw("same");

                            c->Print(Form("%s/Loop5_deltaT_EneRatio_TW_Corr_vs_EneRatioREF_scatter/c_deltaT_EneRatioREF_TW_Corr_vs_EneRatioREF_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c->Print(Form("%s/Loop5_deltaT_EneRatio_TW_Corr_vs_EneRatioREF_scatter/c_deltaT_EneRatioREF_TW_Corr_vs_EneRatioREF_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;




                            // -- Delta T TW + t1fine corr

                            // Ch L:
                            c = new TCanvas(Form("c_h1_deltaT_TW_t1fine_Corr_L_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_TW_t1fine_Corr_L_%s",labelLR_energyBin.c_str()));
                            histo = h1_deltaT_TW_t1fine_Corr_L[index2];
                            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),histo->GetMean()+3.*histo->GetRMS());
                            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
                            histo -> SetTitle(Form(";#Delta T_{L} TW+t1fine corr [ps];entries"));
                            histo -> SetLineColor(kBlue);
                            histo -> SetLineWidth(2);
                            //c -> SetLogy();
                            histo -> Draw();
                            histo -> Write();

                            fitFunc_deltaT_TW_t1fineCorr_L[index2] = new TF1(Form("fitFunc_deltaT_TW_t1fineCorr_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
                            histo -> Fit(fitFunc_deltaT_TW_t1fineCorr_L[index2],"QNRS");
                            histo -> Fit(fitFunc_deltaT_TW_t1fineCorr_L[index2],"QSR+","",fitFunc_deltaT_TW_t1fineCorr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_TW_t1fineCorr_L[index2]->GetParameter(2),fitFunc_deltaT_TW_t1fineCorr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_TW_t1fineCorr_L[index2]->GetParameter(2));
                            histo -> Fit(fitFunc_deltaT_TW_t1fineCorr_L[index2],"QSR+","",fitFunc_deltaT_TW_t1fineCorr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_TW_t1fineCorr_L[index2]->GetParameter(2),fitFunc_deltaT_TW_t1fineCorr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_TW_t1fineCorr_L[index2]->GetParameter(2));

                            fitFunc_deltaT_TW_t1fineCorr_L[index2] -> SetLineColor(kBlack);
                            fitFunc_deltaT_TW_t1fineCorr_L[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_TW_t1fineCorr_L[index2] -> Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_TW_t1fineCorr_L[index2]->GetParameter(1), fitFunc_deltaT_TW_t1fineCorr_L[index2]->GetParameter(2)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");

                            c -> Print(Form("%s/Loop5_deltaT_TW_t1fine_Corr/deltaT_TW_t1fine_Corr_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop5_deltaT_TW_t1fine_Corr/deltaT_TW_t1fine_Corr_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete c;
                            delete latex;



                            // Ch R:
                            c = new TCanvas(Form("c_h1_deltaT_TW_t1fine_Corr_R_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_TW_t1fine_Corr_R_%s",labelLR_energyBin.c_str()));
                            histo = h1_deltaT_TW_t1fine_Corr_R[index2];
                            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),histo->GetMean()+3.*histo->GetRMS());
                            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
                            histo -> SetTitle(Form(";#Delta T_{R} TW+t1fine corr [ps];entries"));
                            histo -> SetLineColor(kBlue);
                            histo -> SetLineWidth(2);
                            //c -> SetLogy();
                            histo -> Draw();
                            histo -> Write();

                            fitFunc_deltaT_TW_t1fineCorr_R[index2] = new TF1(Form("fitFunc_deltaT_TW_t1fineCorr_R_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
                            histo -> Fit(fitFunc_deltaT_TW_t1fineCorr_R[index2],"QNRS");
                            histo -> Fit(fitFunc_deltaT_TW_t1fineCorr_R[index2],"QSR+","",fitFunc_deltaT_TW_t1fineCorr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_TW_t1fineCorr_R[index2]->GetParameter(2),fitFunc_deltaT_TW_t1fineCorr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_TW_t1fineCorr_R[index2]->GetParameter(2));
                            histo -> Fit(fitFunc_deltaT_TW_t1fineCorr_R[index2],"QSR+","",fitFunc_deltaT_TW_t1fineCorr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_TW_t1fineCorr_R[index2]->GetParameter(2),fitFunc_deltaT_TW_t1fineCorr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_TW_t1fineCorr_R[index2]->GetParameter(2));

                            fitFunc_deltaT_TW_t1fineCorr_R[index2] -> SetLineColor(kBlack);
                            fitFunc_deltaT_TW_t1fineCorr_R[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_TW_t1fineCorr_R[index2] -> Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_TW_t1fineCorr_R[index2]->GetParameter(1), fitFunc_deltaT_TW_t1fineCorr_R[index2]->GetParameter(2)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");

                            c -> Print(Form("%s/Loop5_deltaT_TW_t1fine_Corr/deltaT_TW_t1fine_Corr_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop5_deltaT_TW_t1fine_Corr/deltaT_TW_t1fine_Corr_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete c;
                            delete latex;





                            // -- Delta T EneRatio + t1fine corr

                            // Ch L:
                            c = new TCanvas(Form("c_h1_deltaT_EneRatioREF_t1fine_Corr_L_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_EneRatioREF_t1fine_Corr_L_%s",labelLR_energyBin.c_str()));
                            histo = h1_deltaT_EneRatioREF_t1fine_Corr_L[index2];
                            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),histo->GetMean()+3.*histo->GetRMS());
                            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
                            histo -> SetTitle(Form(";#Delta T_{L} EneRatioREF+t1fine corr [ps];entries"));
                            histo -> SetLineColor(kGreen);
                            histo -> SetLineWidth(2);
                            //c -> SetLogy();
                            histo -> Draw();
                            histo -> Write();

                            fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2] = new TF1(Form("fitFunc_deltaT_EneRatioREF_t1fineCorr_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2],"QNRS");
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2],"QSR+","",fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2]->GetParameter(2));
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2],"QSR+","",fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2]->GetParameter(2));

                            fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2] -> SetLineColor(kBlack);
                            fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2] -> Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2]->GetParameter(1), fitFunc_deltaT_EneRatioREF_t1fineCorr_L[index2]->GetParameter(2)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");

                            c -> Print(Form("%s/Loop5_deltaT_EneRatioREF_t1fine_Corr/deltaT_EneRatioREF_t1fine_Corr_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop5_deltaT_EneRatioREF_t1fine_Corr/deltaT_EneRatioREF_t1fine_Corr_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete c;
                            delete latex;



                            // Ch R:
                            c = new TCanvas(Form("c_h1_deltaT_EneRatioREF_t1fine_Corr_R_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_EneRatioREF_t1fine_Corr_R_%s",labelLR_energyBin.c_str()));
                            histo = h1_deltaT_EneRatioREF_t1fine_Corr_R[index2];
                            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),histo->GetMean()+3.*histo->GetRMS());
                            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
                            histo -> SetTitle(Form(";#Delta T_{R} EneRatioREF+t1fine corr [ps];entries"));
                            histo -> SetLineColor(kGreen);
                            histo -> SetLineWidth(2);
                            //c -> SetLogy();
                            histo -> Draw();
                            histo -> Write();

                            fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2] = new TF1(Form("fitFunc_deltaT_EneRatioREF_t1fineCorr_R_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2],"QNRS");
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2],"QSR+","",fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2]->GetParameter(2));
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2],"QSR+","",fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2]->GetParameter(2));

                            fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2] -> SetLineColor(kBlack);
                            fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2] -> Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2]->GetParameter(1), fitFunc_deltaT_EneRatioREF_t1fineCorr_R[index2]->GetParameter(2)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");

                            c -> Print(Form("%s/Loop5_deltaT_EneRatioREF_t1fine_Corr/deltaT_EneRatioREF_t1fine_Corr_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop5_deltaT_EneRatioREF_t1fine_Corr/deltaT_EneRatioREF_t1fine_Corr_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete c;
                            delete latex;




                            // -- profile + TH2F DeltaT (with EneRatio REF + t1fine corr) vs Mean t1fine REF

                            //Ch L:
                            c = new TCanvas(Form("c_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()));
                            c->SetGridy();

                            h2 = h2_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2];
                            h2->SetTitle(";Mean t1fine REF;#Delta T_{L} (EneRatio+t1fine Corr) [ps]");
                            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
                            h2->GetYaxis()->SetRangeUser(CTRMeans_TW_EneRatio_L[index2] - 5*CTRSigmas_TW_EneRatio_L[index2], CTRMeans_TW_EneRatio_L[index2] + 5.*CTRSigmas_TW_EneRatio_L[index2]);
                            h2->Draw("colz");

                            prof = p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2];
                            prof->SetLineColor(kRed);
                            prof->SetLineWidth(3);
                            prof->Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex->SetNDC();
                            latex->SetTextFont(42);
                            latex->SetTextSize(0.04);
                            latex->SetTextColor(kRed);
                            latex->Draw("same");

                            c->Print(Form("%s/Loop5_deltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF_scatter/c_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c->Print(Form("%s/Loop5_deltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF_scatter/c_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;


                            //Ch R:
                            c = new TCanvas(Form("c_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()));
                            c->SetGridy();

                            h2 = h2_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2];
                            h2->SetTitle(";Mean t1fine REF;#Delta T_{R} (EneRatio+t1fine Corr) [ps]");
                            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
                            h2->GetYaxis()->SetRangeUser(CTRMeans_TW_EneRatio_R[index2] - 5*CTRSigmas_TW_EneRatio_R[index2], CTRMeans_TW_EneRatio_R[index2] + 5.*CTRSigmas_TW_EneRatio_R[index2]);
                            h2->Draw("colz");

                            prof = p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2];
                            prof->SetLineColor(kRed);
                            prof->SetLineWidth(3);
                            prof->Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex->SetNDC();
                            latex->SetTextFont(42);
                            latex->SetTextSize(0.04);
                            latex->SetTextColor(kRed);
                            latex->Draw("same");

                            c->Print(Form("%s/Loop5_deltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF_scatter/c_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c->Print(Form("%s/Loop5_deltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF_scatter/c_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;




                            // Profile Delta T (w/ EnergyRatioREF + t1fine corr)  vs Mean t1fineREF

                            //Ch L:
                            c = new TCanvas(Form("c_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()));

                            prof = p1_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF[index2];
                            prof -> SetTitle(Form(";<t1fine_{REF}>;#Delta T_{L} EneRatioREF+t1fine Corr [ps]"));
                            prof -> GetYaxis() -> SetRangeUser(CTRMeans_TW_EneRatio_L[index2]-5.*CTRSigmas_TW_EneRatio_L[index2],CTRMeans_TW_EneRatio_L[index2]+5.*CTRSigmas_TW_EneRatio_L[index2]);
                            prof -> Draw("");
                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");

                            fitFunc_deltaT_EneRatio_t1fineCorr_L[index2] = new TF1(Form("fitFunc_deltaT_EneRatio_t1fineCorr_L_%s",labelLR_energyBin.c_str()),"pol3",300,850);
                            prof -> Fit(fitFunc_deltaT_EneRatio_t1fineCorr_L[index2],"QRS+");
                            fitFunc_deltaT_EneRatio_t1fineCorr_L[index2] -> SetLineColor(kRed);
                            fitFunc_deltaT_EneRatio_t1fineCorr_L[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_EneRatio_t1fineCorr_L[index2] -> Draw("same");

                            c -> Print(Form("%s/Loop5_deltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF/c_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop5_deltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF/c_deltaT_EneRatioREF_t1fine_L_Corr_vs_Meant1fineREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;


                            //Ch R:
                            c = new TCanvas(Form("c_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()));

                            prof = p1_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF[index2];
                            prof -> SetTitle(Form(";<t1fine_{REF}>;#Delta T_{R} EneRatioREF+t1fine Corr [ps]"));
                            prof -> GetYaxis() -> SetRangeUser(CTRMeans_TW_EneRatio_R[index2]- 600.,CTRMeans_TW_EneRatio_R[index2] + 600.);
                            prof -> Draw("");
                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");

                            fitFunc_deltaT_EneRatio_t1fineCorr_R[index2] = new TF1(Form("fitFunc_deltaT_EneRatio_t1fineCorr_R_%s",labelLR_energyBin.c_str()),"pol3",300,850);
                            prof -> Fit(fitFunc_deltaT_EneRatio_t1fineCorr_R[index2],"QRS+");
                            fitFunc_deltaT_EneRatio_t1fineCorr_R[index2] -> SetLineColor(kRed);
                            fitFunc_deltaT_EneRatio_t1fineCorr_R[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_EneRatio_t1fineCorr_R[index2] -> Draw("same");

                            c -> Print(Form("%s/Loop5_deltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF/c_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop5_deltaT_EneRatioREF_t1fine_Corr_vs_Meant1fineREF/c_deltaT_EneRatioREF_t1fine_R_Corr_vs_Meant1fineREF__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;





			    // -- profile + TH2F DeltaT + TW corr + t1fine DUT corr vs Mean t1fine REF 

                            //Ch L:
                            c = new TCanvas(Form("c_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()));
                            c->SetGridy();
                            h2 = h2_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2];
                            h2->SetTitle(";<t1fine_{REF}>;#Delta T_{L} TW + t1fineDUT Corr [ps]");
                            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
                            //h2->GetYaxis()->SetRangeUser(CTRMeans_TW_EneRatio_L[index2] - 5*CTRSigmas_TW_EneRatio_L[index], CTRMeans_TW_EneRatio_L[index] + 5.*CTRSigmas_TW_EneRatio_L[index]);
                            h2->Draw("colz");

                            prof = p1_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2];
                            prof->SetLineColor(kRed);
                            prof->SetLineWidth(3);
                            prof->Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex->SetNDC();
                            latex->SetTextFont(42);
                            latex->SetTextSize(0.04);
                            latex->SetTextColor(kRed);
                            latex->Draw("same");

                            c->Print(Form("%s/Loop5_deltaT_TW_t1fine_Corr_vs_Meant1fineREF/c_deltaT_TW_t1fine_Corr_vs_Meant1fineREF_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c->Print(Form("%s/Loop5_deltaT_TW_t1fine_Corr_vs_Meant1fineREF/c_deltaT_TW_t1fine_Corr_vs_Meant1fineREF_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;



                            //Ch R:
                            c = new TCanvas(Form("c_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()),Form("c_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF_%s",labelLR_energyBin.c_str()));
                            c->SetGridy();
                            h2 = h2_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2];
                            h2->SetTitle(";<t1fine_{REF}>;#Delta T_{R} TW + t1fineDUT Corr [ps]");
                            h2->GetYaxis()->SetRangeUser(h2->GetMean(2) - 600., h2->GetMean(2) + 600.);
                            //h2->GetYaxis()->SetRangeUser(CTRMeans_TW_EneRatio_L[index2] - 5*CTRSigmas_TW_EneRatio_L[index], CTRMeans_TW_EneRatio_L[index] + 5.*CTRSigmas_TW_EneRatio_L[index]);
                            h2->Draw("colz");

                            prof = p1_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2];
                            prof->SetLineColor(kRed);
                            prof->SetLineWidth(3);
                            prof->Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
                            latex->SetNDC();
                            latex->SetTextFont(42);
                            latex->SetTextSize(0.04);
                            latex->SetTextColor(kRed);
                            latex->Draw("same");

                            c->Print(Form("%s/Loop5_deltaT_TW_t1fine_Corr_vs_Meant1fineREF/c_deltaT_TW_t1fine_Corr_vs_Meant1fineREF_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c->Print(Form("%s/Loop5_deltaT_TW_t1fine_Corr_vs_Meant1fineREF/c_deltaT_TW_t1fine_Corr_vs_Meant1fineREF_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete latex;
                            delete c;


	  }
      }
    }









  //------------------------
  //--- 6nd loop over events
  for(auto mapIt : trees)
    {
      ModuleEventClass* anEvent = new ModuleEventClass();
      mapIt.second -> SetBranchAddress("event",&anEvent);

      int nEntries = mapIt.second->GetEntries();
      for(int entry = 0; entry < nEntries; ++entry)
        {
          if( entry%100000 == 0 ){
            std::cout << ">>> TIME-WALK loop 6: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
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
          long long deltaT_L_raw = anEvent->timeL - timeMean_ext;
          long long deltaT_R_raw = anEvent->timeR - timeMean_ext;

          float energyMeanREF = 0.5 *(anEvent->energyL_ext + anEvent->energyR_ext);

          float deltaT_L_raw_low = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);
          float deltaT_L_raw_hig = fitFunc_deltaT_Raw_L[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_L[index2]->GetParameter(2);

          float deltaT_R_raw_low = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) - 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);
          float deltaT_R_raw_hig = fitFunc_deltaT_Raw_R[index2]->GetParameter(1) + 5*fitFunc_deltaT_Raw_R[index2]->GetParameter(2);

          // -- TW corr only fitFunc(energy):
          float TW_corr_L = fitFunc_deltaT_EneLCorr_L[index2]->Eval(anEvent->energyL);
          float TW_corr_R = fitFunc_deltaT_EneRCorr_R[index2]->Eval(anEvent->energyR);

          // -- EneRatioREF correction (from delta T raw):
          float EneRatio_corr_L = fitFunc_deltaT_L_raw_ERatioREF[index2]->Eval(anEvent->energyL/energyMeanREF);
          float EneRatio_corr_R = fitFunc_deltaT_R_raw_ERatioREF[index2]->Eval(anEvent->energyR/energyMeanREF);

          // -- TW + t1fine corr
          float TW_corr_t1fine_L = fitFunc_deltaT_TWCorrL_t1fineL[index2]->Eval(anEvent->t1fineL);
          float TW_corr_t1fine_R = fitFunc_deltaT_TWCorrR_t1fineR[index2]->Eval(anEvent->t1fineR);

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

	  float TW_Mean_t1fineREF_corr_L = GetProfileBinCorrection(p1_deltaT_TW_t1fine_L_Corr_vs_Meant1fineREF[index2],Mean_t1fineREF); //--> corr bin x bin
          float TW_Mean_t1fineREF_corr_R = GetProfileBinCorrection(p1_deltaT_TW_t1fine_R_Corr_vs_Meant1fineREF[index2],Mean_t1fineREF);


          std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));


          // -- delta T histograms
          if( h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2] == NULL )
          {

            // -- EneRatioREF corr + t1fine DUT corr + Mean t1fine REF corr
            h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2] = new TH1F(Form("h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
            h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2] = new TH1F(Form("h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
          
	    // -- TW corr + t1fine DUT corr + Mean t1fine REF corr
	    h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2] = new TH1F(Form("h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
            h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2] = new TH1F(Form("h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R_%s",labelLR_energyBin.c_str()),"",2000,-12000,12000.);
	  }


          if((deltaT_L_raw > deltaT_L_raw_low ) && (deltaT_L_raw < deltaT_L_raw_hig)){
          h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2] -> Fill (deltaT_L_raw - EneRatio_corr_L - EneR_corr_t1fine_L - EneRatioREF_Mean_t1fineREF_corr_L);
	  h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2] -> Fill (deltaT_L_raw - TW_corr_L - TW_corr_t1fine_L - TW_Mean_t1fineREF_corr_L);
	  }


	  if((deltaT_R_raw > deltaT_R_raw_low ) && (deltaT_R_raw < deltaT_R_raw_hig)){
          h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2] -> Fill (deltaT_R_raw - EneRatio_corr_R - EneR_corr_t1fine_R - EneRatioREF_Mean_t1fineREF_corr_R );
	  h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2] -> Fill (deltaT_R_raw - TW_corr_R - TW_corr_t1fine_R - TW_Mean_t1fineREF_corr_R);
	  }

	}
    }




    //------------------
    //--- draw plots loop 5

    std::map<double,TF1*> fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L;
    std::map<double,TF1*> fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R;
    std::map<double,TF1*> fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L;
    std::map<double,TF1*> fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R;

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





                            // -- Delta T EneRatioREF corr + t1fine DUT corr + Mean t1fine REF

                            // Ch L:
                            c = new TCanvas(Form("c_h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L_%s",labelLR_energyBin.c_str()));
                            histo = h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2];
                            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),histo->GetMean()+3.*histo->GetRMS());
                            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
                            histo -> SetTitle(Form(";#Delta T_{L} EneRatio+t1fineDUT+t1fineREF [ps];entries"));
                            histo -> SetLineColor(kGreen);
                            histo -> SetLineWidth(2);
                            //c -> SetLogy();
                            histo -> Draw();
                            histo -> Write();

                            fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2] = new TF1(Form("fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2],"QNRS");
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2],"QSR+","",fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(2));
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2],"QSR+","",fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(2));

                            fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2] -> SetLineColor(kBlack);
                            fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2] -> Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(1), fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(2)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");

                            c -> Print(Form("%s/Loop6_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr/deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop6_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr/deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete c;
                            delete latex;



                            // Ch R:
                            c = new TCanvas(Form("c_h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R_%s",labelLR_energyBin.c_str()));
                            histo = h1_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2];
                            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),histo->GetMean()+3.*histo->GetRMS());
                            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
                            histo -> SetTitle(Form(";#Delta T_{R} EneRatio+t1fineDUT+t1fineREF [ps];entries"));
                            histo -> SetLineColor(kGreen);
                            histo -> SetLineWidth(2);
                            //c -> SetLogy();
                            histo -> Draw();
                            histo -> Write();

                            fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2] = new TF1(Form("fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2],"QNRS");
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2],"QSR+","",fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(2));
                            histo -> Fit(fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2],"QSR+","",fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(2),fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(2));

                            fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2] -> SetLineColor(kBlack);
                            fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2] -> Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(1), fitFunc_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(2)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");

                            c -> Print(Form("%s/Loop6_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr/deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop6_deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr/deltaT_EneRatioREF_t1fine_Meant1fineREF_Corr_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete c;
                            delete latex;








			     // -- Delta T TW corr + t1fine DUT corr + Mean t1fine REF


                            // Ch L:
                            c = new TCanvas(Form("c_h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L_%s",labelLR_energyBin.c_str()));
                            histo = h1_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2];
                            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),histo->GetMean()+3.*histo->GetRMS());
                            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
                            histo -> SetTitle(Form(";#Delta T_{L} TW+t1fineDUT+t1fineREF [ps];entries"));
                            histo -> SetLineColor(kBlue);
                            histo -> SetLineWidth(2);
                            //c -> SetLogy();
                            histo -> Draw();
                            histo -> Write();

                            fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2] = new TF1(Form("fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
                            histo -> Fit(fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2],"QNRS");
                            histo -> Fit(fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2],"QSR+","",fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(2),fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(2));
                            histo -> Fit(fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2],"QSR+","",fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(1)-2.*fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(2),fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(1)+2.*fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(2));

                            fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2] -> SetLineColor(kBlack);
                            fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2] -> Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d L}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(1), fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_L[index2]->GetParameter(2)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");

                            c -> Print(Form("%s/Loop6_deltaT_TW_t1fine_Meant1fineREF_Corr/deltaT_TW_t1fine_Meant1fineREF_Corr_L__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop6_deltaT_TW_t1fine_Meant1fineREF_Corr/deltaT_TW_t1fine_Meant1fineREF_Corr_L__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete c;
                            delete latex;


                            // Ch R:
                            c = new TCanvas(Form("c_h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R_%s",labelLR_energyBin.c_str()),Form("c_h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R_%s",labelLR_energyBin.c_str()));
                            histo = h1_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2];
                            histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-3.*histo->GetRMS(),histo->GetMean()+3.*histo->GetRMS());
                            histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
                            histo -> SetTitle(Form(";#Delta T_{R} TW+t1fineDUT+t1fineREF [ps];entries"));
                            histo -> SetLineColor(kBlue);
                            histo -> SetLineWidth(2);
                            //c -> SetLogy();
                            histo -> Draw();
                            histo -> Write();

                            fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2] = new TF1(Form("fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R_%s",labelLR_energyBin.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());
                            histo -> Fit(fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2],"QNRS");
                            histo -> Fit(fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2],"QSR+","",fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(2),fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(2));
                            histo -> Fit(fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2],"QSR+","",fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(1)-2.*fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(2),fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(1)+2.*fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(2));

                            fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2] -> SetLineColor(kBlack);
                            fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2] -> SetLineWidth(2);
                            fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2] -> Draw("same");

                            latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d R}{#splitline{V_{OV} = %.2f V, th. = %d DAC}""{Mean = %.3f, stdDev = %.3f}}",iBar, Vov, int(vth1), fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(1), fitFunc_deltaT_TW_t1fine_Meant1fineREF_Corr_R[index2]->GetParameter(2)));
                            latex -> SetNDC();
                            latex -> SetTextFont(42);
                            latex -> SetTextSize(0.04);
                            latex -> SetTextColor(kRed);
                            latex -> Draw("same");

                            c -> Print(Form("%s/Loop6_deltaT_TW_t1fine_Meant1fineREF_Corr/deltaT_TW_t1fine_Meant1fineREF_Corr_R__%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
                            c -> Print(Form("%s/Loop6_deltaT_TW_t1fine_Meant1fineREF_Corr/deltaT_TW_t1fine_Meant1fineREF_Corr_R__%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
                            delete c;
                            delete latex;








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
