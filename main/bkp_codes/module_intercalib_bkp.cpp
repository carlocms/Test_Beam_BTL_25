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
#include "macros/CMS_lumi.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <set>
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


// Compute the RMS and the mean from a TGraph plot
void computeMeanRMS(TGraph* g, double& mean, double& rms)
{
  int n = g->GetN();
  if(n == 0)
  {
    mean = 0.0;
    rms  = 0.0;
    return;
  }

  double sum = 0.0;
  double sum2 = 0.0;

  for(int i = 0; i < n; ++i)
  {
    double x, y;
    g->GetPoint(i, x, y);
    sum  += y;
    sum2 += y*y;
  }

  mean = sum / n;
  rms  = std::sqrt(sum2/n - mean*mean);
}



// --- For Landau Fit of the energy histograms 

struct FitResult {
    TF1* func;
    double mpv;
    double chi2_ndf;
};

FitResult DoLandauFit(TH1* histo,
                      double minE,
                      const std::string& name)
{
    FitResult result;
    result.func = nullptr;
    result.mpv = -1;
    result.chi2_ndf = -1;

    if (!histo) return result;

    // --- initial estimates
    float max = histo->GetBinCenter(histo->GetMaximumBin());
    float xmin = std::max((float)minE, max * 0.65f);
    float xmax = std::min(max * 2.5f, 940.f);

    // --- define function
    TF1* f = new TF1(name.c_str(),
                     "[0]*TMath::Landau(x,[1],[2])",
                     0., 1000.);

    // --- initial parameters
    f->SetParameters(
        histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10.,
        max,
        0.1 * max
    );

    f->SetParLimits(1, 0, 9999);
    f->SetParLimits(2, 0, 9999);

    // --- first fit
    histo->Fit(f, "QRS");

    // --- refine fit range
    if (f->GetParameter(1) > 0)
    {
        xmin = f->GetParameter(1) - 2 * std::abs(f->GetParameter(2));
        if (xmin < minE) xmin = minE;

        xmax = std::min(f->GetParameter(1) * 2.5, 940.);

        f->SetRange(xmin, xmax);

        f->SetParameters(
            histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10.,
            f->GetParameter(1),
            0.1 * f->GetParameter(1)
        );
    }

    // --- second fit
    histo->Fit(f, "QRS");

    // --- style
    f->SetLineColor(kBlack);
    f->SetLineWidth(2);

    // --- extract results
    if (f->GetNDF() > 0)
    {
        result.mpv = f->GetParameter(1);
        result.chi2_ndf = f->GetChisquare() / f->GetNDF();
    }

    result.func = f;
    return result;
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


  //--- Directory parameters
  std::string LO_mode = opts.GetOpt<std::string>("Intercalib.mode");
  std::string outputDirTxt = opts.GetOpt<std::string>("Intercalib.OutputDirTxt");
  std::string plotDir_base = opts.GetOpt<std::string>("Output.plotDir_intercalib");
  std::string plotDir = plotDir_base + "/mode_" + LO_mode;

  // --- From Loop 1:
  system(Form("mkdir -p %s",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop1_energy_RAW/",plotDir.c_str()));
 
  // --- From Loop 2:
  system(Form("mkdir -p %s/Loop2_energyLO_calib/Left/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_energyLO_calib/Right/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop2_scatter_energyLO_calib_L_vs_R/",plotDir.c_str()));

  // --- Chi2 check plots:
  system(Form("mkdir -p %s/Chi2_plots/", plotDir.c_str()));

  // --- From calibration factors plots:
  system(Form("mkdir -p %s/LO_calibration_factors/",plotDir.c_str()));
  system(Form("mkdir -p %s/TOFHIR_LO_calibration_factors/",plotDir.c_str())); 
  system(Form("mkdir -p %s/TOFHIR_calibration_factors/",plotDir.c_str()));
  system(Form("mkdir -p %s/Global_calibration_factors/",plotDir.c_str()));

  // --- From Loop 3:
  system(Form("mkdir -p %s/Loop3_energy_TOFHIR_LO_calib/Left/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_energy_TOFHIR_LO_calib/Right/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_energy_TOFHIR_calib/Left/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_energy_TOFHIR_calib/Right/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_scatter_energy_TOFHIR_LO_calib_L_vs_R/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_scatter_energy_TOFHIR_calib_L_vs_R/",plotDir.c_str()));

  // --- From MPV post calibration:
  system(Form("mkdir -p %s/MPV_RAW/",plotDir.c_str()));
  system(Form("mkdir -p %s/MPV_LO_calib/",plotDir.c_str()));
  system(Form("mkdir -p %s/MPV_TOFHIR_LO_calib/",plotDir.c_str()));
  system(Form("mkdir -p %s/MPV_TOFHIR_calib/",plotDir.c_str()));
  system(Form("mkdir -p %s/MPV_summary/", plotDir.c_str()));

  std::vector<std::string> LRLabels;
  LRLabels.push_back("L");
  LRLabels.push_back("R");
  LRLabels.push_back("L-R");

  std::vector<std::string> SideLabels;
  SideLabels.push_back("L");
  SideLabels.push_back("R");

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
          std::cout << "minEnergies:   bar " << bar <<  "   Vov " << ov << "   E " << minE[std::make_pair(bar,ov)] <<std::endl;
        }
    }
  else
    {
      for(unsigned int iBar = 0; iBar < 16; ++iBar)
        for(unsigned int ii = 0; ii < Vov.size(); ++ii)
          minE[std::make_pair(iBar, Vov[ii])] = map_energyMins[Vov[ii]];
    }


//---------------------------------
//--- open file LO intercalibration
//---------------------------------

std::string LOintercalibBaseName = opts.GetOpt<std::string>("Intercalib.LOintercalibFileName");

if( LOintercalibBaseName == "" )
{
    std::cerr << "[ERROR] LOintercalibFileName not provided in cfg" << std::endl;
    exit(1);
}

if( LO_mode != "global" && LO_mode != "separated_side" )
{
    std::cerr << "[ERROR] Invalid LO mode: " << LO_mode << std::endl;
    exit(1);
}

std::string LOintercalibFileName = LOintercalibBaseName + "_" + LO_mode + ".txt";
std::ifstream testFile(LOintercalibFileName);
std::map< std::pair<int, std::string>, double > LO_factors;

std::cout << "Reading LO intercalibration file: " << LOintercalibFileName << std::endl;

if( LOintercalibFileName != "" )
{
    std::ifstream LOintercalibFile(LOintercalibFileName);
    // --- check file open
    if( !LOintercalibFile.is_open() )
    {
        std::cerr << "[ERROR] Cannot open LO calibration file: "<< LOintercalibFileName << std::endl;
        exit(1);
    }

    std::string line;
    int nLines = 0;

    while( std::getline(LOintercalibFile, line) )
    {
        // --- skip empty or comment lines
        if( line.empty() ) continue;
        if( line[0] == '#' ) continue;

        std::istringstream ss(line);

        int bar;
        std::string side;
        float value;

        // --- check parsing success
        if( !(ss >> bar >> side >> value) )
        {
            std::cerr << "[WARNING] Malformed line skipped: " << line << std::endl;
            continue;
        }

        // --- normalize side string (safety)
        std::transform(side.begin(), side.end(), side.begin(), ::toupper);

        // --- check side validity
        if( side != "L" && side != "R" )
        {
            std::cerr << "[WARNING] Invalid side label (must be L or R): "<< side << " -> skipping line: " << line << std::endl;
            continue;
        }

        // --- store value
        auto key = std::make_pair(bar, side);

        // --- check duplicate
        if( LO_factors.find(key) != LO_factors.end() )
        {
            std::cerr << "[WARNING] Duplicate entry for bar " << bar << " side " << side << " -> overwriting previous value"<< std::endl;
        }

        LO_factors[key] = value;

        std::cout << "LO factor loaded: bar " << bar << " side " << side << " value " << value << std::endl;
        nLines++;
    }

    LOintercalibFile.close();

    std::cout << "[INFO] Loaded " << nLines << " LO calibration entries" << std::endl;

    // --- sanity check: expect 16 bars × 2 sides = 32 entries
    if( LO_factors.size() < 32 )
    {
        std::cerr << "[WARNING] Expected ~32 entries, found " << LO_factors.size() << std::endl;
    }
}
else
{
    std::cout << "[INFO] No LO calibration file provided -> using default = 1" << std::endl;

    for(unsigned int iBar = 0; iBar < 16; ++iBar)
    {
        for(auto sideLabel : SideLabels)
        {
            LO_factors[std::make_pair(iBar, sideLabel)] = 1.0;
        }
    }
}
  //----------------------------------------------------------------------------

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
  std::string outFileName = opts.GetOpt<std::string>("Output.outFileNameModuleIntercalib");
  TFile* outFile = TFile::Open(outFileName.c_str(),"RECREATE");
  outFile->cd();

  //Loop 2:
  std::map<double,TH1F*> h1_energy_LOcalib_L;
  std::map<double,TH1F*> h1_energy_LOcalib_R;
  std::map<double,TH2F*> h2_energy_LOcalib_L_vs_R;

  //Calibration factors:
  TGraph* g_LO_L = new TGraph();
  TGraph* g_LO_R = new TGraph();
  std::map< std::pair<float,int>, TGraph* > g_TOFHIR_LO_L;
  std::map< std::pair<float,int>, TGraph* > g_TOFHIR_LO_R;
  std::map< std::pair<float,int>, TGraph* > g_TOFHIR_L;
  std::map< std::pair<float,int>, TGraph* > g_TOFHIR_R;
  std::map< std::pair<float,int>, TGraph* > g_Global_calib_factors_L;
  std::map< std::pair<float,int>, TGraph* > g_Global_calib_factors_R;
  std::map< std::pair<float,int>, int > counter_TOFHIR_L;
  std::map< std::pair<float,int>, int > counter_TOFHIR_R;
  std::map< std::pair<float,int>, int > counter_L;
  std::map< std::pair<float,int>, int > counter_R;
  std::map< std::pair<float,int>, int > counter_Global_L;
  std::map< std::pair<float,int>, int > counter_Global_R;

  //Loop3:
  std::map<double,TH1F*> h1_energy_TOFHIR_LO_calib_L;
  std::map<double,TH1F*> h1_energy_TOFHIR_LO_calib_R;
  std::map<double,TH1F*> h1_energy_TOFHIR_calib_L;
  std::map<double,TH1F*> h1_energy_TOFHIR_calib_R;
  std::map<double,TH2F*> h2_energy_TOFHIR_LO_calib_L_vs_R;
  std::map<double,TH2F*> h2_energy_TOFHIR_calib_L_vs_R;

  // MPV post calibration:
  std::map< std::pair<float,int>, TGraph* > g_MPV_TOFHIR_LO_L;
  std::map< std::pair<float,int>, TGraph* > g_MPV_TOFHIR_LO_R;
  std::map< std::pair<float,int>, TGraph* > g_MPV_TOFHIR_L;
  std::map< std::pair<float,int>, TGraph* > g_MPV_TOFHIR_R;
  std::map< std::pair<float,int>, TGraph* > g_MPV_LO_L;
  std::map< std::pair<float,int>, TGraph* > g_MPV_LO_R;
  std::map< std::pair<float,int>, TGraph* > g_MPV_RAW_L;
  std::map< std::pair<float,int>, TGraph* > g_MPV_RAW_R;
  std::map< std::pair<float,int>, int > counter_MPV_TOFHIR_LO_L;
  std::map< std::pair<float,int>, int > counter_MPV_TOFHIR_LO_R;
  std::map< std::pair<float,int>, int > counter_MPV_TOFHIR_L;
  std::map< std::pair<float,int>, int > counter_MPV_TOFHIR_R;
  std::map< std::pair<float,int>, int > counter_MPV_LO_L;
  std::map< std::pair<float,int>, int > counter_MPV_LO_R;
  std::map< std::pair<float,int>, int > counter_MPV_RAW_L;
  std::map< std::pair<float,int>, int > counter_MPV_RAW_R;
 


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



  ///////////////////////
  //  draw 1st plots   //
  ///////////////////////
  //Retrieve the energy histograms from the root file from step 1 and apply a Landau fit. 
  //Extract energy ranges from the fit.
  //-----------------------------------------------------------------
  std::string source = opts.GetOpt<std::string>("Input.sourceName");
  std::string Na22 = "Na22";
  std::string Na22SingleBar = "Na22SingleBar";
  std::string Co60 = "Co60";
  std::string Co60SumPeak = "Co60SumPeak";
  std::string Laser = "Laser";
  std::string TB = "TB";
  std::string keepAll = "keepAll";
  std::vector<int> barList = opts.GetOpt<std::vector<int> >("Plots.barList");// list of bars to be analyzed read from cfg

  std::map< std::tuple<int,float,int,std::string>, double > MPV_RAW_map;


  for(auto stepLabel : stepLabels)
    {
      float Vov = map_Vovs[stepLabel];
      float vth1 = map_ths[stepLabel];
      std::string VovLabel(Form("Vov%.2f",Vov));
      std::string thLabel(Form("th%02.0f",vth1));

      // --- external bar energy:
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
	  c -> Print(Form("%s/Loop1_energy_RAW/c_energy_external__Vov%.2f_th%02.0f.png",plotDir.c_str(), Vov, vth1));
	  c -> Print(Form("%s/Loop1_energy_RAW/c_energy_external__Vov%.2f_th%02.0f.pdf",plotDir.c_str(), Vov, vth1));
	  delete c;
	  delete ftemp;
	}

      
      // ---DUT loop bars: 
      for(int iBar = 0; iBar < 16; ++iBar) {

	bool barFound = std::find(barList.begin(), barList.end(), iBar) != barList.end() ;
	if (!barFound) continue;

	int index( (10000*int(Vov*100.)) + (100*vth1) + iBar );

	// -- loop over L, R, LR
	for(auto LRLabel : LRLabels ) {

	  //label histo
	  std::string label(Form("bar%02d%s_%s",iBar,LRLabel.c_str(),stepLabel.c_str()));

	  latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d%s RAW}{V_{OV} = %.2f V, th. = %d DAC}",iBar,LRLabel.c_str(),Vov,int(vth1)));
	  if (LRLabel == "L-R") {
	    latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	  }
	  latex -> SetNDC();
	  latex -> SetTextFont(42);
	  latex -> SetTextSize(0.04);
	  latex -> SetTextColor(kRed);

	  
	  // --- draw energy histograms:
	  c = new TCanvas(Form("c_energy_%s",label.c_str()),Form("c_energy_%s",label.c_str()));
	  gPad -> SetLogy();

	  histo = (TH1F*)( inFile->Get(Form("h1_energy_%s",label.c_str())) );
	  if( !histo ) continue;
	  histo -> SetTitle(";energy [a.u.] (raw);entries");
	  histo -> SetLineColor(kRed);
	  histo -> SetLineWidth(2);
	  histo -> Draw();


	  // --- look for peaks and define energy ranges
	  if( source.compare(Na22) && source.compare(Na22SingleBar) && source.compare(Co60) && source.compare(Co60SumPeak) && source.compare(Laser) && source.compare(TB) && source.compare(keepAll) )
	    {
	      std::cout << " Source not found !!! " << std::endl;
	      return(0);
	    }

	  
	  ranges[LRLabel][index] = new std::vector<float>; //energy ranges

	  // -- if MIP peak, we don't use the spectrum analyzers - just langaus fit to the energy peak.
	  if(!source.compare(TB)){
	    
            float max = histo->GetBinCenter(histo->GetMaximumBin());
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

	    if (f_landau[index]->GetNDF() > 0 && LRLabel !="L-R")
	    {
		    double mpv_raw = f_landau[index]->GetParameter(1);
		    MPV_RAW_map[std::make_tuple(iBar,Vov,vth1, LRLabel)] = mpv_raw;
              }

	    f_landau[index] -> SetLineColor(kBlack);
	    f_landau[index] -> SetLineWidth(2);
	    f_landau[index] -> Draw("same");

	    if ( f_landau[index]->GetNDF() >0 && f_landau[index]->GetParameter(1) > minE[std::make_pair(iBar, Vov)] &&
		 (f_landau[index]->GetParameter(1) - 2.0 * std::abs(f_landau[index]->GetParameter(2))) >=  minE[std::make_pair(iBar, Vov)] &&
		 (f_landau[index]->GetParameter(1) - 2.0 * std::abs(f_landau[index]->GetParameter(2))) < 950) {
	      ranges[LRLabel][index] -> push_back( f_landau[index]->GetParameter(1) - 2.0 * std::abs(f_landau[index]->GetParameter(2)));
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
	  c -> Print(Form("%s/Loop1_energy_RAW/c_energy__%s.png",plotDir.c_str(),label.c_str()));
	  c -> Print(Form("%s/Loop1_energy_RAW/c_energy__%s.pdf",plotDir.c_str(),label.c_str()));
	  delete c;
	  delete latex;

	}// end loop over L, R, L-R labels

      }// -- end loop over bars

    } // -- end loop over stepLabels

  // ---  end 1st plots






  //////////////////////////////
  //   2nd loop over events   //
  //////////////////////////////
  //Energy histograms with application of the Light Output equalization factors per channel to individual events.
  //---------------------------------------------

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

      		accept[index1][entry] = false;

      		if(!ranges["L-R"][index1] ) continue;

      		int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

      		if( energyBinAverage < 1 ) continue;

      		accept[index1][entry] = true;

      		double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

      		//--------------------------------------------------------------------------------

      		auto keyL = std::make_pair(anEvent->barID, "L");
	    	auto keyR = std::make_pair(anEvent->barID, "R");
 
      		// safety check
	    	if( LO_factors.find(keyL) == LO_factors.end() || LO_factors.find(keyR) == LO_factors.end() )
		{
		      	std::cerr << "[ERROR] Missing LO factor for bar "<< anEvent->barID << std::endl;
			continue;	 
	       	}

		float energy_LOcalib_L = anEvent->energyL * LO_factors[keyL]; // energy equalization w.r.t the LO factor 	  
		float energy_LOcalib_R = anEvent->energyR * LO_factors[keyR];

		if( h1_energy_LOcalib_L[index2] == NULL )
		{
		std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
	        h1_energy_LOcalib_L[index2] = new TH1F(Form("h1_energy_LOcalib_L_%s",labelLR_energyBin.c_str()),"",512, 0., 1024.);
                h1_energy_LOcalib_R[index2] = new TH1F(Form("h1_energy_LOcalib_R_%s",labelLR_energyBin.c_str()),"",512, 0., 1024.);
                h2_energy_LOcalib_L_vs_R[index2] = new TH2F(Form("h2_energy_LOcalib_L_vs_R_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024., 512, 0., 1024.);
		}

      		h1_energy_LOcalib_L[index2] -> Fill(energy_LOcalib_L);
      		h1_energy_LOcalib_R[index2] -> Fill(energy_LOcalib_R);
      		h2_energy_LOcalib_L_vs_R[index2] -> Fill( energy_LOcalib_L, energy_LOcalib_R);

	} // end loop over entries
}





  ////////////////////////
  //   draw 2nd plots   //
  ////////////////////////
  

// --- maps:

std::map<double,TF1*> fitFunc_energy_LOcalib_L;
std::map<double,TF1*> fitFunc_energy_LOcalib_R;
    
std::map< std::tuple<int,float,int,std::string>, double > MPV_map;           // map [bar, vov, vth, side] -> MPV values
std::map< std::pair<float,int>, double > MPV_mean_map;                       // map [vov, vth] -> <MPV> all the channels 
std::map< std::tuple<float,int,std::string>, double > MPV_mean_map_side;
std::map< std::tuple<float,int,std::string>, std::vector<double> > MPV_values_map_side; //caso separated side
std::map< std::pair<float,int>, std::vector<double> > MPV_values_map;        // map [vov, vth] -> vector(MPVs)  to achieve <MPV> all the channels

std::map< std::tuple<int,float,int,std::string>, double > TOFHIR_LO_factors; // map [bar, vov, vth, side] -> <MPV> / MPV  
std::map< std::tuple<int,float,int,std::string>, double > TOFHIR_factors;      // map [bar, vov, vth, side] -> TOFHIR_LO_factor / TOFHIR_factor

std::map< std::tuple<int,float,int,std::string>, double > Global_factors;
std::map< std::tuple<int,float,int,std::string>, double > Chi2_values_LO_map;
 

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

	      if (!h1_energy_LOcalib_L[index2]) continue;
	      std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));

	      //----------------------------------------------------------------------------------------------------------

	      // --- energy distribution Side L + LO equalization:
	      c = new TCanvas(Form("c_energyLO_calib_L_%s",labelLR_energyBin.c_str()),Form("c_energyLO_calib_L_%s",labelLR_energyBin.c_str()));
	      histo = h1_energy_LOcalib_L[index2];
	      //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
	      histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
	      histo -> SetTitle(Form(";E_{L} + LO calib [a.u.];entries"));
	      histo -> SetLineColor(kRed);
	      histo -> SetLineWidth(2);
	      histo -> Draw();
	      writeExtraText = true;
              extraText = "Test Beam";
              CMS_lumi(c, 0, 0);
	      histo -> Write();

/*
	      // --- Landau fit
	      float max_L = histo->GetBinCenter(histo->GetMaximumBin());
	      double eneMin_L = minE[std::make_pair(iBar, Vov)];
	      float xmin_L = std::max(eneMin_L, max_L * 0.65);
	      float xmax_L = std::min(max_L * 2.5, 940.);

	      fitFunc_energy_LOcalib_L[index2] = new TF1(Form("f_landau_LO_bar%02d_%s", iBar, labelLR_energyBin.c_str()),"[0]*TMath::Landau(x,[1],[2])",0., 1000.);//0->xmin_L, 1000->xmax_L

	      fitFunc_energy_LOcalib_L[index2]->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, max_L,0.1 * max_L);

	      fitFunc_energy_LOcalib_L[index2]->SetParLimits(1, 0, 9999);
	      fitFunc_energy_LOcalib_L[index2]->SetParLimits(2, 0, 9999);

	      histo->Fit(fitFunc_energy_LOcalib_L[index2], "QRS"); //first fit

	      if (fitFunc_energy_LOcalib_L[index2]->GetParameter(1) > 0) // refine range
	      {  
		      xmin_L = fitFunc_energy_LOcalib_L[index2]->GetParameter(1)- 2 * std::abs(fitFunc_energy_LOcalib_L[index2]->GetParameter(2));

		      if (xmin_L < minE[std::make_pair(iBar, Vov)]) xmin_L = minE[std::make_pair(iBar, Vov)];
   
		      xmax_L = std::min(fitFunc_energy_LOcalib_L[index2]->GetParameter(1) * 2.5,940.);

		      fitFunc_energy_LOcalib_L[index2]->SetRange(xmin_L, xmax_L);
		      fitFunc_energy_LOcalib_L[index2]->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10,fitFunc_energy_LOcalib_L[index2]->GetParameter(1),0.1 * fitFunc_energy_LOcalib_L[index2]->GetParameter(1));
	      }

	      histo->Fit(fitFunc_energy_LOcalib_L[index2], "QRS"); //second fit

	      fitFunc_energy_LOcalib_L[index2]->SetLineColor(kBlack);
	      fitFunc_energy_LOcalib_L[index2]->SetLineWidth(2);
	      fitFunc_energy_LOcalib_L[index2]->Draw("same");

	      // --- maps fill:
	      if (fitFunc_energy_LOcalib_L[index2]->GetNDF() > 0)
	      {
		      double mpv_L = fitFunc_energy_LOcalib_L[index2]->GetParameter(1);
		      double Chi2_NDF_L = (fitFunc_energy_LOcalib_L[index2]->GetChisquare())/(fitFunc_energy_LOcalib_L[index2]->GetNDF());
    		      MPV_map[std::make_tuple(iBar,Vov,vth1, "L")] = mpv_L;
		      MPV_values_map[std::make_pair(Vov,int(vth1))].push_back(mpv_L);
		      MPV_values_map_side[std::make_tuple(Vov,vth1,"L")].push_back(mpv_L);
		      Chi2_values_LO_map[std::make_tuple(iBar,Vov,vth1, "L")] = Chi2_NDF_L;
		      //std::cout << "[MPV] bar " << iBar << " thr " << vth1 << " side L -> " << mpv_L << std::endl;
	      }
*/



	      auto fitRes = DoLandauFit(histo,minE[std::make_pair(iBar, Vov)],Form("f_landau_LO_bar%02d_%s", iBar, labelLR_energyBin.c_str()));


	      // draw
	      if (fitRes.func)fitRes.func->Draw("same");


	      // fill maps

	      if (fitRes.mpv > 0)

	      {
		      MPV_map[std::make_tuple(iBar,Vov,vth1,"L")] = fitRes.mpv;
		      MPV_values_map[std::make_pair(Vov,int(vth1))].push_back(fitRes.mpv);
		      MPV_values_map_side[std::make_tuple(Vov,vth1,"L")].push_back(fitRes.mpv);
		      Chi2_values_LO_map[std::make_tuple(iBar,Vov,vth1,"L")] = fitRes.chi2_ndf;
	      }




	      latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d L - with LO calibration}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
	      latex -> SetNDC();
	      latex -> SetTextFont(42);
	      latex -> SetTextSize(0.04);
	      latex -> SetTextColor(kBlack);
	      latex -> Draw("same");

	      c -> Print(Form("%s/Loop2_energyLO_calib/Left/c_energyLO_calib_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	      c -> Print(Form("%s/Loop2_energyLO_calib/Left/c_energyLO_calib_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
	      delete latex;
	      delete c;



              // --- energy distribution Side R + LO equalization:
              c = new TCanvas(Form("c_energyLO_calib_R_%s",labelLR_energyBin.c_str()),Form("c_energyLO_calib_R_%s",labelLR_energyBin.c_str()));
              histo = h1_energy_LOcalib_R[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";E_{R} + LO calib [a.u.];entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              // --- Landau fit:
              float max_R = histo->GetBinCenter(histo->GetMaximumBin());
              double eneMin_R = minE[std::make_pair(iBar, Vov)];
              float xmin_R = std::max(eneMin_R, max_R * 0.65);
              float xmax_R = std::min(max_R * 2.5, 940.);

              fitFunc_energy_LOcalib_R[index2] = new TF1(Form("f_landau_LO_bar%02d_%s", iBar, labelLR_energyBin.c_str()),"[0]*TMath::Landau(x,[1],[2])",0., 1000.);//0->xmin_R, 1000->xmax_R

              fitFunc_energy_LOcalib_R[index2]->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, max_R,0.1 * max_R);

              fitFunc_energy_LOcalib_R[index2]->SetParLimits(1, 0, 9999);
              fitFunc_energy_LOcalib_R[index2]->SetParLimits(2, 0, 9999);

              histo->Fit(fitFunc_energy_LOcalib_R[index2], "QRS");

              if (fitFunc_energy_LOcalib_R[index2]->GetParameter(1) > 0)
              {
                      xmin_R = fitFunc_energy_LOcalib_R[index2]->GetParameter(1)- 2 * std::abs(fitFunc_energy_LOcalib_R[index2]->GetParameter(2));

                      if (xmin_R < minE[std::make_pair(iBar, Vov)]) xmin_R = minE[std::make_pair(iBar, Vov)];

                      xmax_R = std::min(fitFunc_energy_LOcalib_R[index2]->GetParameter(1) * 2.5,940.);

                      fitFunc_energy_LOcalib_R[index2]->SetRange(xmin_R, xmax_R);
                      fitFunc_energy_LOcalib_R[index2]->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10,fitFunc_energy_LOcalib_R[index2]->GetParameter(1),0.1 * fitFunc_energy_LOcalib_R[index2]->GetParameter(1));
              }

              histo->Fit(fitFunc_energy_LOcalib_R[index2], "QRS");

              fitFunc_energy_LOcalib_R[index2]->SetLineColor(kBlack);
              fitFunc_energy_LOcalib_R[index2]->SetLineWidth(2);
              fitFunc_energy_LOcalib_R[index2]->Draw("same");

	      // --- maps fill:
	      if (fitFunc_energy_LOcalib_R[index2]->GetNDF() > 0)
              {
                      double mpv_R = fitFunc_energy_LOcalib_R[index2]->GetParameter(1);
		      double Chi2_NDF_R = (fitFunc_energy_LOcalib_R[index2]->GetChisquare())/(fitFunc_energy_LOcalib_R[index2]->GetNDF());
                      MPV_map[std::make_tuple(iBar,Vov,vth1, "R")] = mpv_R;
		      MPV_values_map[std::make_pair(Vov,int(vth1))].push_back(mpv_R);
		      MPV_values_map_side[std::make_tuple(Vov,vth1,"R")].push_back(mpv_R);
                      Chi2_values_LO_map[std::make_tuple(iBar,Vov,vth1, "R")] = Chi2_NDF_R;
                      //std::cout << "[MPV] bar " << iBar << " thr " << vth1 << " side R -> " << mpv_R << std::endl;
              }

              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d R - with LO calibration}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop2_energyLO_calib/Right/c_energyLO_calib_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop2_energyLO_calib/Right/c_energyLO_calib_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;



	      // --- scatter plot enrgy side L vs side R with LO equalization:
	      if(!h2_energy_LOcalib_L_vs_R[index2]) continue;

               c = new TCanvas(Form("c_h2_energyLO_calib_L_vs_R_%s",labelLR_energyBin.c_str()),Form("c_h2_energyLO_calib_L_vs_R_%s",labelLR_energyBin.c_str()));
               c -> SetGridy();

               h2 = h2_energy_LOcalib_L_vs_R[index2];
               //h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) -600., h2 -> GetMean(2)+600);
               h2 -> SetTitle(Form(";E_{L} + LO calib [a.u.];E_{R} + LO calib [a.u.]"));
               h2 -> Draw("colz");

	       latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
               latex -> SetNDC();
               latex -> SetTextFont(42);
               latex -> SetTextSize(0.04);
               latex -> SetTextColor(kBlack);
               latex -> Draw("same");

	       c -> Print(Form("%s/Loop2_scatter_energyLO_calib_L_vs_R/c_energyLO_calib_L_vs_R_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
	       c -> Print(Form("%s/Loop2_scatter_energyLO_calib_L_vs_R/c_energyLO_calib_L_vs_R_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
               delete c;
               delete latex;


            } // --- end loop over energy bins

        } // --- end loop ober bars

    } // --- end loop over stepLabels





///////////////////////////////////////////
//   Calibration factors determination   //
///////////////////////////////////////////

//--------------------------------------------------------------------
// Averaging MPV values <-> from energy distribution + LO equalization
//--------------------------------------------------------------------

if (LO_mode == "global")
{
    for (auto &it : MPV_values_map)
    {
        float vov = it.first.first;
        int thr   = it.first.second;
        auto &vec = it.second;

        if (vec.empty()) continue;

        double sum = 0.;
        int n = 0;

        for (auto val : vec)
        {
            if (val < 50 || val > 900) continue;
            sum += val;
            n++;
        }

        if (n == 0) continue;

        double mean = sum / n;
        MPV_mean_map[it.first] = mean;
    }
}
else if (LO_mode == "separated_side")
{
    for (auto &it : MPV_values_map_side)
    {
        float vov;
        int thr;
        std::string side;

        std::tie(vov, thr, side) = it.first;

        auto &vec = it.second;
        if (vec.empty()) continue;

        double sum = 0.;
        int n = 0;

        for (auto val : vec)
        {
            if (val < 50 || val > 900) continue;
            sum += val;
            n++;
        }

        if (n == 0) continue;

        double mean = sum / n;
        MPV_mean_map_side[it.first] = mean;
    }
}




//---------------------------------------
// Determination of the TOFHIR_LO_factors
//---------------------------------------

for (auto &it : MPV_map)
{
    int bar;
    float vov;
    int thr;
    std::string side;

    std::tie(bar, vov, thr, side) = it.first;
    double mpv = it.second;

    if (mpv <= 0) continue;

    double mean_mpv = -1;

    if (LO_mode == "global")
    {
        auto key = std::make_pair(vov, thr);
        if (!MPV_mean_map.count(key)) continue;
        mean_mpv = MPV_mean_map[key];
    }
    else if (LO_mode == "separated_side")
    {
        auto key = std::make_tuple(vov, thr, side);
        if (!MPV_mean_map_side.count(key)) continue;
        mean_mpv = MPV_mean_map_side[key];
    }

    if (mean_mpv <= 0) continue;

    double calib = mean_mpv / mpv;

    TOFHIR_LO_factors[it.first] = calib;

    std::cout << "[LO+TOFHIR CALIB] bar " << bar  << " thr " << thr << " side " << side << " --> " << calib << std::endl;
}




//------------------------------------
// Determination of the TOFHIR_factors
//------------------------------------

for (auto &it : TOFHIR_LO_factors)
{
    int bar;
    float vov;
    int thr;
    std::string side;

    std::tie(bar, vov, thr, side) = it.first;
    double TOFHIR_LO_calib = it.second;

    if (TOFHIR_LO_calib <= 0) continue;

    auto key_LO = std::make_pair(bar, side);
    if (!LO_factors.count(key_LO)) continue;

    double LO_factor_value = LO_factors[key_LO];

    double TOFHIR_calib = TOFHIR_LO_calib / LO_factor_value;

    TOFHIR_factors[it.first] = TOFHIR_calib;

    std::cout << "[TOFHIR CALIB] bar " << bar << " thr " << thr << " side " << side << " --> " << TOFHIR_calib << std::endl;
}



//------------------------------------------------
// Determination of the Global calibration factors
//------------------------------------------------

for (auto &it : TOFHIR_LO_factors)
{
    int bar;
    float vov;
    int thr;
    std::string side;

    std::tie(bar, vov, thr, side) = it.first;
    double TOFHIR_LO_calib = it.second;

    if (TOFHIR_LO_calib <= 0) continue;

    auto key_LO = std::make_pair(bar, side);
    if (!LO_factors.count(key_LO)) continue;

    double LO_factor_value = LO_factors[key_LO];

    double Global_calib = TOFHIR_LO_calib * LO_factor_value;

    Global_factors[it.first] = Global_calib;

    if (thr == 12) std::cout << "[Global CALIB] bar " << bar << " thr " << thr << " side " << side << " --> " << Global_calib << std::endl;
}





//------------------------------------
//   Txt output files calib factors
//------------------------------------

std::set< std::pair<float,int> > vov_thr_set;

for (auto &it : TOFHIR_LO_factors)
{
    float vov = std::get<1>(it.first);
    int   thr = std::get<2>(it.first);

    vov_thr_set.insert(std::make_pair(vov,thr));
}

// --- LOOP SU (Vov, thr)
for (auto &vt : vov_thr_set)
{
    float vov = vt.first;
    int   thr = vt.second;

    // --- FILE TOFHIR+LO:
    std::string outFileName_LO = Form("%s/TOFHIR_LO_calibration_factors_Vov%.2f_th%02d_%s.txt",outputDirTxt.c_str(), vov, thr, LO_mode.c_str());

    std::ofstream outFile_LO(outFileName_LO);

    if (!outFile_LO.is_open())
    {
        std::cerr << "ERROR: cannot open output file " << outFileName_LO << std::endl;
    }
    else
    {
        outFile_LO << "# mode: " << LO_mode << "\n";
        outFile_LO << "# bar side calib\n";

        // --- ordina per bar
        std::map<int, std::map<std::string,double>> ordered_LO;

        for (auto &it : TOFHIR_LO_factors)
        {
            int bar;
            float vov_i;
            int thr_i;
            std::string side;

            std::tie(bar, vov_i, thr_i, side) = it.first;

            if (vov_i != vov || thr_i != thr) continue;

            ordered_LO[bar][side] = it.second;
        }

        for (auto &b : ordered_LO)
        {
            for (auto &s : b.second)
            {
                outFile_LO << b.first << " " << s.first << " " << s.second << "\n";
            }
        }

        outFile_LO.close();
        std::cout << "[WRITE] TOFHIR+LO file saved: " << outFileName_LO << std::endl;
    }

    // --- FILE TOFHIR:
    std::string outFileName = Form(
        "%s/TOFHIR_calibration_factors_Vov%.2f_th%02d_%s.txt",
        outputDirTxt.c_str(), vov, thr, LO_mode.c_str()
    );

    std::ofstream outFile(outFileName);

    if (!outFile.is_open())
    {
        std::cerr << "ERROR: cannot open output file " << outFileName << std::endl;
    }
    else
    {
        outFile << "# mode: " << LO_mode << "\n";
        outFile << "# bar side calib\n";

        // --- ordina per bar
        std::map<int, std::map<std::string,double>> ordered;

        for (auto &it : TOFHIR_factors)
        {
            int bar;
            float vov_i;
            int thr_i;
            std::string side;

            std::tie(bar, vov_i, thr_i, side) = it.first;

            if (vov_i != vov || thr_i != thr) continue;

            ordered[bar][side] = it.second;
        }

        for (auto &b : ordered)
        {
            for (auto &s : b.second)
            {
                outFile << b.first << " " << s.first << " " << s.second << "\n";
            }
        }

        outFile.close();
        std::cout << "[WRITE] TOFHIR file saved: " << outFileName << std::endl;
    }
//}


    // --- FILE GLOBAL:
    std::string outFileName_Global = Form("%s/GLOBAL_calibration_factors_Vov%.2f_th%02d_%s.txt",outputDirTxt.c_str(), vov, thr, LO_mode.c_str());

    std::ofstream outFile_Global(outFileName_Global);

    if (!outFile_Global.is_open())
    {
        std::cerr << "ERROR: cannot open output file " << outFileName_Global << std::endl;
    }
    else
    {
        outFile << "# mode: " << LO_mode << "\n";
        outFile << "# bar side calib\n";

        // --- ordina per bar
        std::map<int, std::map<std::string,double>> ordered_Global;

        for (auto &it : Global_factors)
        {
            int bar;
            float vov_i;
            int thr_i;
            std::string side;

            std::tie(bar, vov_i, thr_i, side) = it.first;

            if (vov_i != vov || thr_i != thr) continue;


            ordered_Global[bar][side] = it.second;
        }

        for (auto &b : ordered_Global)
        {
            for (auto &s : b.second)
            {
                outFile_Global << b.first << " " << s.first << " " << s.second << "\n";
            }
        }

        outFile_Global.close();
        std::cout << "[WRITE] GLOBAL file saved: " << outFileName_Global << std::endl;
    }
}


//////////////////////////////////////////
//   Plot Callibration factors vs bar   //
//////////////////////////////////////////

//-------------------------------------------------------------------
//  MPV values from energy histograms withou any calibration (raw):
//------------------------------------------------------------------- 

for(auto& it : MPV_RAW_map)
{
  int bar;
  float vov;
  int thr;
  std::string side;
  double val;

  std::tie(bar, vov, thr, side) = it.first;
  val = it.second;

  std::pair<float,int> key = std::make_pair(vov, thr);

  // crea grafici se non esistono
  if( g_MPV_RAW_L.find(key) == g_MPV_RAW_L.end() )
  {
    g_MPV_RAW_L[key] = new TGraph();
    g_MPV_RAW_R[key] = new TGraph();
    counter_MPV_RAW_L[key] = 0;
    counter_MPV_RAW_R[key] = 0;
  }

  if(side == "L")
  {
    g_MPV_RAW_L[key]->SetPoint(counter_MPV_RAW_L[key], bar, val);
    counter_MPV_RAW_L[key]++;
  }
  else if(side == "R")
  {
    g_MPV_RAW_R[key]->SetPoint(counter_MPV_RAW_R[key], bar, val);
    counter_MPV_RAW_R[key]++;
  }
}


for(auto& it : g_MPV_RAW_L)
{
  float vov = it.first.first;
  int thr   = it.first.second;

  TGraph* gL = g_MPV_RAW_L[it.first];
  TGraph* gR = g_MPV_RAW_R[it.first];

  TCanvas* c = new TCanvas(Form("c_MPV_RAW_vov%.2f_th%d",vov,thr),Form("MPV w/o calib Vov=%.2f th=%d",vov,thr),800,600);

  gL->SetMarkerStyle(20);
  gL->SetMarkerColor(kRed);
  gL->SetLineColor(kRed);
  gL->SetLineWidth(2);
  gL->GetYaxis()->SetRangeUser(230.,440.);

  gR->SetMarkerStyle(21);
  gR->SetMarkerColor(kBlue);
  gR->SetLineColor(kBlue);
  gR->SetLineWidth(2);

  gL->SetTitle(Form(";Bar index; MPV raw (w/o calib) [a.u.]"));

  gL->Draw("APL");
  gR->Draw("PL SAME");


  // Mean and RMS values
  double mean_L, rms_L;
  double mean_R, rms_R;

  computeMeanRMS(gL, mean_L, rms_L);
  computeMeanRMS(gR, mean_R, rms_R);

  TLegend* leg = new TLegend(0.27,0.70,0.65,0.9);
  leg->SetHeader(Form("V_{OV} = %.2f V, vth = %i DAC", vov, thr), "C");
  leg->AddEntry(gL, Form("Side L (Mean=%.4f, RMS=%.4f)", mean_L, rms_L), "lp");
  leg->AddEntry(gR, Form("Side R (Mean=%.4f, RMS=%.4f)", mean_R, rms_R), "lp");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->Draw();

  c->SetGridy();
  c->SetTickx();
  c->SetTicky();

  c->Print(Form("%s/MPV_RAW/MPV_RAW_vov%.2f_th%d.pdf",plotDir.c_str(),vov,thr));
  c->Print(Form("%s/MPV_RAW/MPV_RAW_vov%.2f_th%d.png",plotDir.c_str(),vov,thr));

  delete leg;
  delete c;
}




//-------------------------------------------------------
// MPV values histograms energy + LO equalization
//-------------------------------------------------------

for(auto& it : MPV_map)
{
  int bar;
  float vov;
  int thr;
  std::string side;
  double val;

  std::tie(bar, vov, thr, side) = it.first;
  val = it.second;

  std::pair<float,int> key = std::make_pair(vov, thr);

  if( g_MPV_LO_L.find(key) == g_MPV_LO_L.end() )
  {
    g_MPV_LO_L[key] = new TGraph();
    g_MPV_LO_R[key] = new TGraph();
    counter_MPV_LO_L[key] = 0;
    counter_MPV_LO_R[key] = 0;
  }

  if(side == "L")
  {
    g_MPV_LO_L[key]->SetPoint(counter_MPV_LO_L[key], bar, val);
    counter_MPV_LO_L[key]++;
  }
  else if(side == "R")
  {
    g_MPV_LO_R[key]->SetPoint(counter_MPV_LO_R[key], bar, val);
    counter_MPV_LO_R[key]++;
  }
}


for(auto& it : g_MPV_LO_L)
{
  float vov = it.first.first;
  int thr   = it.first.second;

  TGraph* gL = g_MPV_LO_L[it.first];
  TGraph* gR = g_MPV_LO_R[it.first];

  TCanvas* c = new TCanvas(Form("c_MPV_LO_vov%.2f_th%d",vov,thr),Form("MPV with LO calib Vov=%.2f th=%d",vov,thr),800,600);

  // stile L
  gL->SetMarkerStyle(20);
  gL->SetMarkerColor(kRed);
  gL->SetLineColor(kRed);
  gL->SetLineWidth(2);
  gL->GetYaxis()->SetRangeUser(230.,440.);

  gR->SetMarkerStyle(21);
  gR->SetMarkerColor(kBlue);
  gR->SetLineColor(kBlue);
  gR->SetLineWidth(2);

  gL->SetTitle(Form(";Bar index; MPV w/ LO calib [a.u.]"));
  gL->Draw("APL");
  gR->Draw("PL SAME");


  // Mean and RMS values
  double mean_L, rms_L;
  double mean_R, rms_R;

  computeMeanRMS(gL, mean_L, rms_L);
  computeMeanRMS(gR, mean_R, rms_R);

  TLegend* leg = new TLegend(0.27,0.70,0.65,0.9);
  leg->SetHeader(Form("V_{OV} = %.2f V, vth = %i DAC", vov, thr), "C");
  leg->AddEntry(gL, Form("Side L (Mean=%.4f, RMS=%.4f)", mean_L, rms_L), "lp");
  leg->AddEntry(gR, Form("Side R (Mean=%.4f, RMS=%.4f)", mean_R, rms_R), "lp");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->Draw();

  c->SetGridy();
  c->SetTickx();
  c->SetTicky();

  c->Print(Form("%s/MPV_LO_calib/MPV_LO_calib_vov%.2f_th%d.pdf",plotDir.c_str(),vov,thr));
  c->Print(Form("%s/MPV_LO_calib/MPV_LO_calib_vov%.2f_th%d.png",plotDir.c_str(),vov,thr));

  delete leg;
  delete c;
}




//------------------------------------
// Plot of the LO equalization factors
//------------------------------------

int nL = 0;
int nR = 0;

for(int iBar = 0; iBar < 16; ++iBar)
{
  // --- L side
  if( LO_factors.count(std::make_pair(iBar,"L")) )
  {
    double val = LO_factors[std::make_pair(iBar,"L")];
    g_LO_L->SetPoint(nL, iBar, val);
    nL++;
  }

  // --- R side
  if( LO_factors.count(std::make_pair(iBar,"R")) )
  {
    double val = LO_factors[std::make_pair(iBar,"R")];
    g_LO_R->SetPoint(nR, iBar, val);
    nR++;
  }
}

// --- canvas
TCanvas* c_LO = new TCanvas("c_LO_factors","LO factors vs bar",800,600);
c_LO->cd();

g_LO_L->SetMarkerStyle(20);
g_LO_L->SetMarkerColor(kRed);
g_LO_L->SetLineColor(kRed);
g_LO_L->SetLineWidth(2);
g_LO_L->GetYaxis()->SetRangeUser(0.93, 1.09);

g_LO_R->SetMarkerStyle(21);
g_LO_R->SetMarkerColor(kBlue);
g_LO_R->SetLineColor(kBlue);
g_LO_R->SetLineWidth(2);

g_LO_L->SetTitle(";Bar index;LO calibration factor");

g_LO_L->Draw("APL");   // A = axes, P = points, L = line
g_LO_R->Draw("PL SAME");

// Mean and RMS values
double mean_L, rms_L;
double mean_R, rms_R;

computeMeanRMS(g_LO_L, mean_L, rms_L);
computeMeanRMS(g_LO_R, mean_R, rms_R);

TLegend* leg = new TLegend(0.20,0.70,0.65,0.9);
leg->AddEntry(g_LO_L, Form("Side L (<LO>_{L}=%.4f, RMS=%.4f)", mean_L, rms_L), "lp");
leg->AddEntry(g_LO_R, Form("Side R (<LO>_{R}=%.4f, RMS=%.4f)", mean_R, rms_R), "lp");
leg->SetFillStyle(0);
leg->SetBorderSize(0);
leg->SetTextSize(0.03);
leg->SetTextFont(42);
leg->Draw();

c_LO->SetGridy();
c_LO->SetTickx();
c_LO->SetTicky();

c_LO->Print(Form("%s/LO_calibration_factors/LO_factors_vs_bar.pdf", plotDir.c_str()));
c_LO->Print(Form("%s/LO_calibration_factors/LO_factors_vs_bar.png", plotDir.c_str()));
delete leg;
delete c_LO;






//-------------------------------
//  Plot of the TOFHIR+LO factors
//-------------------------------

for(auto& it : TOFHIR_LO_factors)
{
  int bar;
  float vov;
  int thr;
  std::string side;
  double val;

  std::tie(bar, vov, thr, side) = it.first;
  val = it.second;

  std::pair<float,int> key = std::make_pair(vov, thr);

  // crea grafici se non esistono
  if( g_TOFHIR_LO_L.find(key) == g_TOFHIR_LO_L.end() )
  {
    g_TOFHIR_LO_L[key] = new TGraph();
    g_TOFHIR_LO_R[key] = new TGraph();
    counter_L[key] = 0;
    counter_R[key] = 0;
  }

  if(side == "L")
  {
    g_TOFHIR_LO_L[key]->SetPoint(counter_L[key], bar, val);
    counter_L[key]++;
  }
  else if(side == "R")
  {
    g_TOFHIR_LO_R[key]->SetPoint(counter_R[key], bar, val);
    counter_R[key]++;
  }
}



for(auto& it : g_TOFHIR_LO_L)
{
  float vov = it.first.first;
  int thr   = it.first.second;

  TGraph* gL = g_TOFHIR_LO_L[it.first];
  TGraph* gR = g_TOFHIR_LO_R[it.first];

  TCanvas* c = new TCanvas(Form("c_TOFHIR_LO_vov%.2f_th%d",vov,thr),Form("TOFHIR+LO factors Vov=%.2f th=%d",vov,thr),800,600);

  gL->SetMarkerStyle(20);
  gL->SetMarkerColor(kRed);
  gL->SetLineColor(kRed);
  gL->SetLineWidth(2);
  gL->GetYaxis()->SetRangeUser(0.5, 1.5);//0.75, 1.35

  gR->SetMarkerStyle(21);
  gR->SetMarkerColor(kBlue);
  gR->SetLineColor(kBlue);
  gR->SetLineWidth(2);

  gL->SetTitle(Form(";Bar index;TOFHIR calibration"));

  gL->Draw("APL");
  gR->Draw("PL SAME");

  // Mean and RMS values
  double mean_L, rms_L;
  double mean_R, rms_R;

  computeMeanRMS(gL, mean_L, rms_L);
  computeMeanRMS(gR, mean_R, rms_R);

  TLegend* leg = new TLegend(0.27,0.70,0.65,0.9);
  leg->SetHeader(Form("V_{OV} = %.2f V, vth = %i DAC", vov, thr), "C");
  leg->AddEntry(gL, Form("Side L (Mean=%.4f, RMS=%.4f)", mean_L, rms_L), "lp");
  leg->AddEntry(gR, Form("Side R (Mean=%.4f, RMS=%.4f)", mean_R, rms_R), "lp");
  leg->SetFillStyle(0);  
  leg->SetBorderSize(0);  
  leg->SetTextSize(0.03);  
  leg->SetTextFont(42);
  leg->Draw();

  c->SetGridy();
  c->SetTickx();
  c->SetTicky();

  c->Print(Form("%s/TOFHIR_LO_calibration_factors/TOFHIR_LO_factors_vov%.2f_th%d.pdf",plotDir.c_str(),vov,thr));
  c->Print(Form("%s/TOFHIR_LO_calibration_factors/TOFHIR_LO_factors_vov%.2f_th%d.png",plotDir.c_str(),vov,thr));

  delete leg;
  delete c;
}




//-------------------------------
//  Plot of the Global factors
//-------------------------------

for(auto& it : Global_factors)
{
  int bar;
  float vov;
  int thr;
  std::string side;
  double val;

  std::tie(bar, vov, thr, side) = it.first;
  val = it.second;

  std::pair<float,int> key = std::make_pair(vov, thr);

  // crea grafici se non esistono
  if( g_Global_calib_factors_L.find(key) == g_Global_calib_factors_L.end() )
  { 
    g_Global_calib_factors_L[key] = new TGraph();
    g_Global_calib_factors_R[key] = new TGraph();
    counter_Global_L[key] = 0;
    counter_Global_R[key] = 0;
  }

  if(side == "L")
  {
    g_Global_calib_factors_L[key]->SetPoint(counter_Global_L[key], bar, val);
    counter_Global_L[key]++;
  }
  else if(side == "R")
  {
    g_Global_calib_factors_R[key]->SetPoint(counter_Global_R[key], bar, val);
    counter_Global_R[key]++;
  }
}


for(auto& it : g_Global_calib_factors_L)
{
  float vov = it.first.first;
  int thr   = it.first.second;

  TGraph* gL = g_Global_calib_factors_L[it.first];
  TGraph* gR = g_Global_calib_factors_R[it.first];

  TCanvas* c = new TCanvas(Form("c_Global_vov%.2f_th%d",vov,thr),Form("Global factors Vov=%.2f th=%d",vov,thr),800,600);

  gL->SetMarkerStyle(20);
  gL->SetMarkerColor(kRed);
  gL->SetLineColor(kRed);
  gL->SetLineWidth(2);
  gL->GetYaxis()->SetRangeUser(0.5, 1.5);//0.75, 1.35

  gR->SetMarkerStyle(21);
  gR->SetMarkerColor(kBlue);
  gR->SetLineColor(kBlue);
  gR->SetLineWidth(2);

  gL->SetTitle(Form(";Bar index;Global calib factors"));

  gL->Draw("APL");
  gR->Draw("PL SAME");

  // Mean and RMS values
  double mean_L, rms_L;
  double mean_R, rms_R;

  computeMeanRMS(gL, mean_L, rms_L);
  computeMeanRMS(gR, mean_R, rms_R);

  TLegend* leg = new TLegend(0.27,0.70,0.65,0.9);
  leg->SetHeader(Form("V_{OV} = %.2f V, vth = %i DAC", vov, thr), "C");
  leg->AddEntry(gL, Form("Side L (Mean=%.4f, RMS=%.4f)", mean_L, rms_L), "lp");
  leg->AddEntry(gR, Form("Side R (Mean=%.4f, RMS=%.4f)", mean_R, rms_R), "lp");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->Draw();

  c->SetGridy();
  c->SetTickx();
  c->SetTicky();

  c->Print(Form("%s/Global_calibration_factors/Global_factors_vov%.2f_th%d.pdf",plotDir.c_str(),vov,thr));
  c->Print(Form("%s/Global_calibration_factors/Global_factors_vov%.2f_th%d.png",plotDir.c_str(),vov,thr));

  delete leg;
  delete c;
}




//---------------------------
// Plot of the TOFHIR factors
//---------------------------
 
for(auto& it : TOFHIR_factors)
{
  int bar;
  float vov;
  int thr;
  std::string side;
  double val;

  std::tie(bar, vov, thr, side) = it.first;
  val = it.second;

  std::pair<float,int> key = std::make_pair(vov, thr);

  // crea grafici se non esistono
  if( g_TOFHIR_L.find(key) == g_TOFHIR_L.end() )
  {
    g_TOFHIR_L[key] = new TGraph();
    g_TOFHIR_R[key] = new TGraph();
    counter_TOFHIR_L[key] = 0;
    counter_TOFHIR_R[key] = 0;
  }


  if(side == "L")
  {
    g_TOFHIR_L[key]->SetPoint(counter_TOFHIR_L[key], bar, val);
    counter_TOFHIR_L[key]++;
  }
  else if(side == "R")
  {
    g_TOFHIR_R[key]->SetPoint(counter_TOFHIR_R[key], bar, val);
    counter_TOFHIR_R[key]++;
  }
}


for(auto& it : g_TOFHIR_L)
{
  float vov = it.first.first;
  int thr   = it.first.second;

  TGraph* gL = g_TOFHIR_L[it.first];
  TGraph* gR = g_TOFHIR_R[it.first];

  TCanvas* c = new TCanvas(Form("c_TOFHIR_vov%.2f_th%d",vov,thr),Form("TOFHIR factors Vov=%.2f th=%d",vov,thr),800,600);

  gL->SetMarkerStyle(20);
  gL->SetMarkerColor(kRed);
  gL->SetLineColor(kRed);
  gL->SetLineWidth(2);
  gL->GetYaxis()->SetRangeUser(0.75, 1.35);

  gR->SetMarkerStyle(21);
  gR->SetMarkerColor(kBlue);
  gR->SetLineColor(kBlue);
  gR->SetLineWidth(2);
  
  gL->SetTitle(Form(";Bar index;TOFHIR calibration factor"));

  gL->Draw("APL");
  gR->Draw("PL SAME");

  double mean_L, rms_L;
  double mean_R, rms_R;

  computeMeanRMS(gL, mean_L, rms_L);
  computeMeanRMS(gR, mean_R, rms_R);

  TLegend* leg = new TLegend(0.27,0.70,0.65,0.9);
  leg->SetHeader(Form("V_{OV} = %.2f V, vth = %i DAC", vov, thr), "C");
  leg->AddEntry(gL, Form("Side L (Mean=%.4f, RMS=%.4f)", mean_L, rms_L), "lp");
  leg->AddEntry(gR, Form("Side R (Mean=%.4f, RMS=%.4f)", mean_R, rms_R), "lp");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->Draw();

  c->SetGridy();
  c->SetTickx();
  c->SetTicky();

  c->Print(Form("%s/TOFHIR_calibration_factors/TOFHIR_factors_vov%.2f_th%d.pdf",plotDir.c_str(),vov,thr));
  c->Print(Form("%s/TOFHIR_calibration_factors/TOFHIR_factors_vov%.2f_th%d.png",plotDir.c_str(),vov,thr));

  delete leg;
  delete c;
}









//////////////////////////////
//   3nd loop over events   //
//////////////////////////////

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


	  auto keyLO_L = std::make_pair(anEvent->barID, "L");
          auto keyLO_R = std::make_pair(anEvent->barID, "R");

          if( LO_factors.find(keyLO_L) == LO_factors.end() || LO_factors.find(keyLO_R) == LO_factors.end() )
          {
                  std::cerr << "[ERROR] Missing LO factor for bar "<< anEvent->barID << std::endl;
                  continue;
          }
	  
	  
	  auto key_TOFHIR_LO_L = std::make_tuple(anEvent->barID, anEvent->Vov, anEvent->vth1, "L");
          auto key_TOFHIR_LO_R = std::make_tuple(anEvent->barID, anEvent->Vov, anEvent->vth1, "R");
 
          if( TOFHIR_LO_factors.find(key_TOFHIR_LO_L) == TOFHIR_LO_factors.end() || TOFHIR_LO_factors.find(key_TOFHIR_LO_R) == TOFHIR_LO_factors.end() )
          {
                  std::cerr << "[ERROR] Missing TOFHIR LO factor for bar "<< anEvent->barID << std::endl;
                  continue;
          }


	  auto key_TOFHIR_L = std::make_tuple(anEvent->barID, anEvent->Vov, anEvent->vth1, "L");
          auto key_TOFHIR_R = std::make_tuple(anEvent->barID, anEvent->Vov, anEvent->vth1, "R");

          if( TOFHIR_factors.find(key_TOFHIR_L) == TOFHIR_factors.end() || TOFHIR_factors.find(key_TOFHIR_R) == TOFHIR_factors.end() )
          {
                  std::cerr << "[ERROR] Missing TOFHIR factor for bar "<< anEvent->barID << std::endl;
                  continue;
          }


	  float energy_TOFHIR_LO_calib_L = (anEvent->energyL * LO_factors[keyLO_L]) * TOFHIR_LO_factors[key_TOFHIR_LO_L];
          float energy_TOFHIR_LO_calib_R = (anEvent->energyR * LO_factors[keyLO_R]) * TOFHIR_LO_factors[key_TOFHIR_LO_R];
 
	  //float energy_TOFHIR_calib_L = (anEvent->energyL) * TOFHIR_factors[key_TOFHIR_L];
          //float energy_TOFHIR_calib_R = (anEvent->energyR) * TOFHIR_factors[key_TOFHIR_R];

	  float energy_TOFHIR_calib_L = (anEvent->energyL) * TOFHIR_LO_factors[key_TOFHIR_LO_L];
	  float energy_TOFHIR_calib_R = (anEvent->energyR) * TOFHIR_LO_factors[key_TOFHIR_LO_R];

          if( h1_energy_TOFHIR_LO_calib_L[index2] == NULL )
            {
              std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
	 
	      h1_energy_TOFHIR_LO_calib_L[index2] = new TH1F(Form("h1_energy_TOFHIR_LO_calib_L_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024.);
              h1_energy_TOFHIR_LO_calib_R[index2] = new TH1F(Form("h1_energy_TOFHIR_LO_calib_R_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024.);
	      h2_energy_TOFHIR_LO_calib_L_vs_R[index2] = new TH2F(Form("h2_energy_TOFHIR_LO_calib_L_vs_R_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024., 512, 0., 1024.);
 
	      h1_energy_TOFHIR_calib_L[index2] = new TH1F(Form("h1_energy_TOFHIR_calib_L_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024.);
              h1_energy_TOFHIR_calib_R[index2] = new TH1F(Form("h1_energy_TOFHIR_calib_R_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024.);
	      h2_energy_TOFHIR_calib_L_vs_R[index2] = new TH2F(Form("h2_energy_TOFHIR_calib_L_vs_R_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024., 512, 0., 1024.);
            }

          
	  h1_energy_TOFHIR_LO_calib_L[index2] -> Fill(energy_TOFHIR_LO_calib_L);
	  h1_energy_TOFHIR_LO_calib_R[index2] -> Fill(energy_TOFHIR_LO_calib_R);
	  h2_energy_TOFHIR_LO_calib_L_vs_R[index2] -> Fill(energy_TOFHIR_LO_calib_L, energy_TOFHIR_LO_calib_R);

	  h1_energy_TOFHIR_calib_L[index2] -> Fill(energy_TOFHIR_calib_L);
          h1_energy_TOFHIR_calib_R[index2] -> Fill(energy_TOFHIR_calib_R);
	  h2_energy_TOFHIR_calib_L_vs_R[index2] -> Fill(energy_TOFHIR_calib_L, energy_TOFHIR_calib_R);

          }
      std::cout << std::endl;
    }






  /////////////////////////
  //   draw 3nd plots   //
  ////////////////////////

  std::map<double,TF1*> fitFunc_energy_TOFHIR_LO_calib_L;
  std::map<double,TF1*> fitFunc_energy_TOFHIR_LO_calib_R;
  std::map<double,TF1*> fitFunc_energy_TOFHIR_calib_L;
  std::map<double,TF1*> fitFunc_energy_TOFHIR_calib_R;

  std::map< std::tuple<int,float,int,std::string>, double > MPV_TOFHIR_LO_map;
  std::map< std::tuple<int,float,int,std::string>, double > MPV_TOFHIR_map;

  std::map< std::tuple<int,float,int,std::string>, double > Chi2_values_TOFHIR_LO_map;
  std::map< std::tuple<int,float,int,std::string>, double > Chi2_values_TOFHIR_map;

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

              if (!h1_energy_LOcalib_L[index2]) continue;
              std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));


              // -- energy histogram side L: LO + TOFHIR calibration
              c = new TCanvas(Form("c_energy_TOFHIR_LO_calib_L_%s",labelLR_energyBin.c_str()),Form("c_energy_TOFHIR_LO_calib_L_%s",labelLR_energyBin.c_str()));
              histo = h1_energy_TOFHIR_LO_calib_L[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";E_{L}+LO+TOFHIR calib [a.u.];entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              float max_L = histo->GetBinCenter(histo->GetMaximumBin());

              double eneMin_L = minE[std::make_pair(iBar, Vov)];
              float xmin_L = std::max(eneMin_L, max_L * 0.65);
              float xmax_L = std::min(max_L * 2.5, 940.);

              fitFunc_energy_TOFHIR_LO_calib_L[index2] = new TF1(Form("f_landau_TOFHIR_LO_bar%02d_%s", iBar, labelLR_energyBin.c_str()),"[0]*TMath::Landau(x,[1],[2])",0., 1000.);

              fitFunc_energy_TOFHIR_LO_calib_L[index2]->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, max_L,0.1 * max_L);

              fitFunc_energy_TOFHIR_LO_calib_L[index2]->SetParLimits(1, 0, 9999);
              fitFunc_energy_TOFHIR_LO_calib_L[index2]->SetParLimits(2, 0, 9999);

              histo->Fit(fitFunc_energy_TOFHIR_LO_calib_L[index2], "QRS");

              if (fitFunc_energy_TOFHIR_LO_calib_L[index2]->GetParameter(1) > 0)
              {
                      xmin_L = fitFunc_energy_TOFHIR_LO_calib_L[index2]->GetParameter(1)- 2 * std::abs(fitFunc_energy_TOFHIR_LO_calib_L[index2]->GetParameter(2));

                      if (xmin_L < minE[std::make_pair(iBar, Vov)]) xmin_L = minE[std::make_pair(iBar, Vov)];

                      xmax_L = std::min(fitFunc_energy_TOFHIR_LO_calib_L[index2]->GetParameter(1) * 2.5,940.);

                      fitFunc_energy_TOFHIR_LO_calib_L[index2]->SetRange(xmin_L, xmax_L);
                      fitFunc_energy_TOFHIR_LO_calib_L[index2]->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10,fitFunc_energy_TOFHIR_LO_calib_L[index2]->GetParameter(1),0.1 * fitFunc_energy_TOFHIR_LO_calib_L[index2]->GetParameter(1));
              }

              histo->Fit(fitFunc_energy_TOFHIR_LO_calib_L[index2], "QRS");
              fitFunc_energy_TOFHIR_LO_calib_L[index2]->SetLineColor(kBlack);
              fitFunc_energy_TOFHIR_LO_calib_L[index2]->SetLineWidth(2);
              fitFunc_energy_TOFHIR_LO_calib_L[index2]->Draw("same");

              if (fitFunc_energy_TOFHIR_LO_calib_L[index2]->GetNDF() > 0)
              {
                      double mpv_L = fitFunc_energy_TOFHIR_LO_calib_L[index2]->GetParameter(1);
                      double Chi2_NDF_L = (fitFunc_energy_TOFHIR_LO_calib_L[index2]->GetChisquare())/(fitFunc_energy_TOFHIR_LO_calib_L[index2]->GetNDF());
		      MPV_TOFHIR_LO_map[std::make_tuple(iBar,Vov,vth1, "L")] = mpv_L;
                      Chi2_values_TOFHIR_LO_map[std::make_tuple(iBar,Vov,vth1, "L")] = Chi2_NDF_L;
              }

              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d L - TOFHIR+LO calibration}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_energy_TOFHIR_LO_calib/Left/c_energyTOFHIR_LO_calib_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop3_energy_TOFHIR_LO_calib/Left/c_energyTOFHIR_LO_calib_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;


              // -- energy histogram side R: LO + TOFHIR calibration
              c = new TCanvas(Form("c_energy_TOFHIR_LO_calib_R_%s",labelLR_energyBin.c_str()),Form("c_energy_TOFHIR_LO_calib_R_%s",labelLR_energyBin.c_str()));
              histo = h1_energy_TOFHIR_LO_calib_R[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";E_{R}+LO+TOFHIR calib [a.u.];entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();
 
              float max_R = histo->GetBinCenter(histo->GetMaximumBin());
              double eneMin_R = minE[std::make_pair(iBar, Vov)];
              float xmin_R = std::max(eneMin_R, max_R * 0.65);
              float xmax_R = std::min(max_R * 2.5, 940.);

              fitFunc_energy_TOFHIR_LO_calib_R[index2] = new TF1(Form("f_landau_TOFHIR_LO_bar%02d_%s", iBar, labelLR_energyBin.c_str()),"[0]*TMath::Landau(x,[1],[2])",0., 1000.);

              fitFunc_energy_TOFHIR_LO_calib_R[index2]->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, max_R,0.1 * max_R);

              fitFunc_energy_TOFHIR_LO_calib_R[index2]->SetParLimits(1, 0, 9999);
              fitFunc_energy_TOFHIR_LO_calib_R[index2]->SetParLimits(2, 0, 9999);

              histo->Fit(fitFunc_energy_TOFHIR_LO_calib_R[index2], "QRS");

              if (fitFunc_energy_TOFHIR_LO_calib_R[index2]->GetParameter(1) > 0)
              {
                      xmin_R = fitFunc_energy_TOFHIR_LO_calib_R[index2]->GetParameter(1)- 2 * std::abs(fitFunc_energy_TOFHIR_LO_calib_R[index2]->GetParameter(2));

                      if (xmin_R < minE[std::make_pair(iBar, Vov)]) xmin_R = minE[std::make_pair(iBar, Vov)];

                      xmax_R = std::min(fitFunc_energy_TOFHIR_LO_calib_R[index2]->GetParameter(1) * 2.5,940.);

                      fitFunc_energy_TOFHIR_LO_calib_R[index2]->SetRange(xmin_R, xmax_R);
                      fitFunc_energy_TOFHIR_LO_calib_R[index2]->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10,fitFunc_energy_TOFHIR_LO_calib_R[index2]->GetParameter(1),0.1 * fitFunc_energy_TOFHIR_LO_calib_R[index2]->GetParameter(1));
              }

              histo->Fit(fitFunc_energy_TOFHIR_LO_calib_R[index2], "QRS");

              fitFunc_energy_TOFHIR_LO_calib_R[index2]->SetLineColor(kBlack);
              fitFunc_energy_TOFHIR_LO_calib_R[index2]->SetLineWidth(2);
              fitFunc_energy_TOFHIR_LO_calib_R[index2]->Draw("same");
              
              if (fitFunc_energy_TOFHIR_LO_calib_R[index2]->GetNDF() > 0)
              {
                      double mpv_R = fitFunc_energy_TOFHIR_LO_calib_R[index2]->GetParameter(1);
                      double Chi2_NDF_R = (fitFunc_energy_TOFHIR_LO_calib_R[index2]->GetChisquare())/(fitFunc_energy_TOFHIR_LO_calib_R[index2]->GetNDF());
		      MPV_TOFHIR_LO_map[std::make_tuple(iBar,Vov,vth1, "R")] = mpv_R;
		      Chi2_values_TOFHIR_LO_map[std::make_tuple(iBar,Vov,vth1, "R")] = Chi2_NDF_R;
              }

              
              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d R - TOFHIR+LO calibration}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_energy_TOFHIR_LO_calib/Right/c_energyTOFHIR_LO_calib_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop3_energy_TOFHIR_LO_calib/Right/c_energyTOFHIR_LO_calib_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;



              // --- scatter plot energy LO+TOFHIR calib L vs R side
              if(!h2_energy_TOFHIR_LO_calib_L_vs_R[index2]) continue;
              c = new TCanvas(Form("c_h2_energy_TOFHIR_LO_calib_L_vs_R_%s",labelLR_energyBin.c_str()),Form("c_h2_energy_TOFHIR_LO_calib_L_vs_R_%s",labelLR_energyBin.c_str()));
              c -> SetGridy();

              h2 = h2_energy_TOFHIR_LO_calib_L_vs_R[index2];
              //h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) -600., h2 -> GetMean(2)+600);
              h2 -> SetTitle(Form(";E_{L} (LO+TOFHIR) [a.u.];E_{R} (LO+TOFHIR) [a.u.]"));
              h2 -> Draw("colz");

              latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d - TOFHIR+LO calib}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_scatter_energy_TOFHIR_LO_calib_L_vs_R/c_energy_TOFHIR_LO_calib_L_vs_R_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop3_scatter_energy_TOFHIR_LO_calib_L_vs_R/c_energy_TOFHIR_LO_calib_L_vs_R_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete c;
              delete latex;




             // -- energy histrogram side L: only TOFHIR calibration
              c = new TCanvas(Form("c_energy_TOFHIR_calib_L_%s",labelLR_energyBin.c_str()),Form("c_energy_TOFHIR_calib_L_%s",labelLR_energyBin.c_str()));
              histo = h1_energy_TOFHIR_calib_L[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";E_{L}+TOFHIR calib [a.u.];entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              float maxTOFHIR_L = histo->GetBinCenter(histo->GetMaximumBin());
              double eneMinTOFHIR_L = minE[std::make_pair(iBar, Vov)];
              float xminTOFHIR_L = std::max(eneMinTOFHIR_L, maxTOFHIR_L * 0.65);
              float xmaxTOFHIR_L = std::min(maxTOFHIR_L * 2.5, 940.);

              fitFunc_energy_TOFHIR_calib_L[index2] = new TF1(Form("f_landau_TOFHIR_bar%02d_%s", iBar, labelLR_energyBin.c_str()),"[0]*TMath::Landau(x,[1],[2])",0., 1000.);

              fitFunc_energy_TOFHIR_calib_L[index2]->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, maxTOFHIR_L,0.1 * maxTOFHIR_L);

              fitFunc_energy_TOFHIR_calib_L[index2]->SetParLimits(1, 0, 9999);
              fitFunc_energy_TOFHIR_calib_L[index2]->SetParLimits(2, 0, 9999);

              histo->Fit(fitFunc_energy_TOFHIR_calib_L[index2], "QRS");

              if (fitFunc_energy_TOFHIR_calib_L[index2]->GetParameter(1) > 0)
              {
                      xminTOFHIR_L = fitFunc_energy_TOFHIR_calib_L[index2]->GetParameter(1)- 2 * std::abs(fitFunc_energy_TOFHIR_calib_L[index2]->GetParameter(2));

                      if (xminTOFHIR_L < minE[std::make_pair(iBar, Vov)]) xminTOFHIR_L = minE[std::make_pair(iBar, Vov)];

                      xmaxTOFHIR_L = std::min(fitFunc_energy_TOFHIR_calib_L[index2]->GetParameter(1) * 2.5,940.);

                      fitFunc_energy_TOFHIR_calib_L[index2]->SetRange(xminTOFHIR_L, xmaxTOFHIR_L);
                      fitFunc_energy_TOFHIR_calib_L[index2]->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10,fitFunc_energy_TOFHIR_calib_L[index2]->GetParameter(1),0.1 * fitFunc_energy_TOFHIR_calib_L[index2]->GetParameter(1));
              }

              histo->Fit(fitFunc_energy_TOFHIR_calib_L[index2], "QRS");
              fitFunc_energy_TOFHIR_calib_L[index2]->SetLineColor(kBlack);
              fitFunc_energy_TOFHIR_calib_L[index2]->SetLineWidth(2);
              fitFunc_energy_TOFHIR_calib_L[index2]->Draw("same");
              
              if (fitFunc_energy_TOFHIR_calib_L[index2]->GetNDF() > 0)
              {
                      double mpv_TOFHIR_L = fitFunc_energy_TOFHIR_calib_L[index2]->GetParameter(1);
                      double Chi2_NDF_TOFHIR_L = (fitFunc_energy_TOFHIR_calib_L[index2]->GetChisquare())/(fitFunc_energy_TOFHIR_calib_L[index2]->GetNDF());
		      MPV_TOFHIR_map[std::make_tuple(iBar,Vov,vth1, "L")] = mpv_TOFHIR_L;
		      Chi2_values_TOFHIR_map[std::make_tuple(iBar,Vov,vth1, "L")] = Chi2_NDF_TOFHIR_L;
              }

              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d L - TOFHIR+calibration}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_energy_TOFHIR_calib/Left/c_energyTOFHIR_calib_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop3_energy_TOFHIR_calib/Left/c_energyTOFHIR_calib_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;




             // -- energy histogram side R: TOFHIR calibration
             c = new TCanvas(Form("c_energy_TOFHIR_calib_R_%s",labelLR_energyBin.c_str()),Form("c_energy_TOFHIR_calib_R_%s",labelLR_energyBin.c_str()));
             histo = h1_energy_TOFHIR_calib_R[index2];
             //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
             histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
             histo -> SetTitle(Form(";E_{R}+TOFHIR calib [a.u.];entries"));
             histo -> SetLineColor(kRed);
             histo -> SetLineWidth(2);
             histo -> Draw();
             histo -> Write();

             float maxTOFHIR_R = histo->GetBinCenter(histo->GetMaximumBin());

             double eneMinTOFHIR_R = minE[std::make_pair(iBar, Vov)];
             float xminTOFHIR_R = std::max(eneMinTOFHIR_R, maxTOFHIR_R * 0.65);
             float xmaxTOFHIR_R = std::min(maxTOFHIR_R * 2.5, 940.);

             fitFunc_energy_TOFHIR_calib_R[index2] = new TF1(Form("f_landau_TOFHIR_bar%02d_%s", iBar, labelLR_energyBin.c_str()),"[0]*TMath::Landau(x,[1],[2])",0., 1000.);

             fitFunc_energy_TOFHIR_calib_R[index2]->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10, maxTOFHIR_R,0.1 * maxTOFHIR_R);

             fitFunc_energy_TOFHIR_calib_R[index2]->SetParLimits(1, 0, 9999);
             fitFunc_energy_TOFHIR_calib_R[index2]->SetParLimits(2, 0, 9999);

             histo->Fit(fitFunc_energy_TOFHIR_calib_R[index2], "QRS");

             if (fitFunc_energy_TOFHIR_calib_R[index2]->GetParameter(1) > 0)
              
	     {
                      xminTOFHIR_R = fitFunc_energy_TOFHIR_calib_R[index2]->GetParameter(1)- 2 * std::abs(fitFunc_energy_TOFHIR_calib_R[index2]->GetParameter(2));

	     	      if (xminTOFHIR_R < minE[std::make_pair(iBar, Vov)]) xminTOFHIR_R = minE[std::make_pair(iBar, Vov)];

                      xmaxTOFHIR_R = std::min(fitFunc_energy_TOFHIR_calib_R[index2]->GetParameter(1) * 2.5,940.);

                      fitFunc_energy_TOFHIR_calib_R[index2]->SetRange(xminTOFHIR_R, xmaxTOFHIR_R);
                      fitFunc_energy_TOFHIR_calib_R[index2]->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10,fitFunc_energy_TOFHIR_calib_R[index2]->GetParameter(1),0.1 * fitFunc_energy_TOFHIR_calib_R[index2]->GetParameter(1));
              }

	     histo->Fit(fitFunc_energy_TOFHIR_calib_R[index2], "QRS");

             fitFunc_energy_TOFHIR_calib_R[index2]->SetLineColor(kBlack);
             fitFunc_energy_TOFHIR_calib_R[index2]->SetLineWidth(2);
             fitFunc_energy_TOFHIR_calib_R[index2]->Draw("same");

	     if (fitFunc_energy_TOFHIR_calib_R[index2]->GetNDF() > 0)
              {       
                      double mpv_TOFHIR_R = fitFunc_energy_TOFHIR_calib_R[index2]->GetParameter(1);
                      double Chi2_NDF_TOFHIR_R = (fitFunc_energy_TOFHIR_calib_R[index2]->GetChisquare())/(fitFunc_energy_TOFHIR_calib_R[index2]->GetNDF());
		      MPV_TOFHIR_map[std::make_tuple(iBar,Vov,vth1, "R")] = mpv_TOFHIR_R;
		      Chi2_values_TOFHIR_map[std::make_tuple(iBar,Vov,vth1, "R")] = Chi2_NDF_TOFHIR_R;
              }


	     latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d R - TOFHIR+calibration}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
             latex -> SetNDC();
             latex -> SetTextFont(42);
             latex -> SetTextSize(0.04);
             latex -> SetTextColor(kBlack);
             latex -> Draw("same");

             c -> Print(Form("%s/Loop3_energy_TOFHIR_calib/Right/c_energyTOFHIR_calib_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
             c -> Print(Form("%s/Loop3_energy_TOFHIR_calib/Right/c_energyTOFHIR_calib_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
             delete latex;
             delete c;



             // --- scatter plot energy only TOFHIR calib L vs R side
             if(!h2_energy_TOFHIR_calib_L_vs_R[index2]) continue;

             c = new TCanvas(Form("c_h2_energy_TOFHIR_calib_L_vs_R_%s",labelLR_energyBin.c_str()),Form("c_h2_energy_TOFHIR_calib_L_vs_R_%s",labelLR_energyBin.c_str()));
             c -> SetGridy();

             h2 = h2_energy_TOFHIR_calib_L_vs_R[index2];
             //h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) -600., h2 -> GetMean(2)+600);
             h2 -> SetTitle(Form(";E_{L} (TOFHIR) [a.u.];E_{R} (TOFHIR) [a.u.]"));
             h2 -> Draw("colz");

             latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d - TOFHIR calib}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
             latex -> SetNDC();
             latex -> SetTextFont(42);
             latex -> SetTextSize(0.04);
             latex -> SetTextColor(kBlack);
             latex -> Draw("same");

             c -> Print(Form("%s/Loop3_scatter_energy_TOFHIR_calib_L_vs_R/c_energy_TOFHIR_calib_L_vs_R_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
             c -> Print(Form("%s/Loop3_scatter_energy_TOFHIR_calib_L_vs_R/c_energy_TOFHIR_calib_L_vs_R_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
             delete c;
             delete latex;

            } // --- end loop over energy bins

        } // --- end loop ober bars

    } // --- end loop over stepLabels







//////////////////////////////
//     Plot Chi2 values     //
//////////////////////////////

auto FillChi2Histos = [&](const std::map<std::tuple<int,float,int,std::string>, double>& chi2_map,
                         const std::string& label)
{
    std::set<std::pair<float,int>> vov_thr_set;

    for (const auto& it : chi2_map)
    {
        int bar;
        float vov;
        int thr;
        std::string side;

        std::tie(bar, vov, thr, side) = it.first;
        vov_thr_set.insert(std::make_pair(vov, thr));
    }

    // --- loop on (vov,thr)
    for (const auto& vt : vov_thr_set)
    {
        float vov = vt.first;
        int thr   = vt.second;

        TH1F* h_chi2 = new TH1F(Form("h_chi2_%s_Vov%.2f_th%02d", label.c_str(), vov, thr),Form(";#chi^{2}/NDF;Entries"),50, 0., 10.);

        for (const auto& it : chi2_map)
        {
            int bar;
            float vov_i;
            int thr_i;
            std::string side;

            std::tie(bar, vov_i, thr_i, side) = it.first;

            if (vov_i != vov || thr_i != thr) continue;

            double chi2 = it.second;

            if (chi2 <= 0 || chi2 > 50) continue; // selection for outliers

            h_chi2->Fill(chi2);
        }

        TCanvas* c = new TCanvas(Form("c_chi2_%s_Vov%.2f_th%02d", label.c_str(), vov, thr),"",800, 600);

        h_chi2->SetLineWidth(2);
        h_chi2->SetLineColor(kBlack);
        h_chi2->Draw("hist");

        double mean = h_chi2->GetMean();
        double rms  = h_chi2->GetRMS();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.60,0.85,Form("Mean = %.2f", mean));
        latex.DrawLatex(0.60,0.80,Form("RMS  = %.2f", rms));

        latex.DrawLatex(0.20,0.85,Form("V_{OV} = %.2f V, thr = %d", vov, thr));

        c->SetGridy();

        c->Print(Form("%s/Chi2_plots/c_chi2_%s_Vov%.2f_th%02d.pdf",plotDir.c_str(), label.c_str(), vov, thr));
        c->Print(Form("%s/Chi2_plots/c_chi2_%s_Vov%.2f_th%02d.png",plotDir.c_str(), label.c_str(), vov, thr));
        delete c;
    }
};

FillChi2Histos(Chi2_values_LO_map,           "LO");
FillChi2Histos(Chi2_values_TOFHIR_LO_map, "LO_TOFHIR");
FillChi2Histos(Chi2_values_TOFHIR_map,    "TOFHIR");







////////////////////////////////////////////////////////
//   Plot energy MPV values vs bar POST CALIBRATION   //
////////////////////////////////////////////////////////


//------------------------------
// MPV TOFHIR + LO coliabration
//------------------------------ 
for(auto& it : MPV_TOFHIR_LO_map)
{
  int bar;
  float vov;
  int thr;
  std::string side;
  double val;

  std::tie(bar, vov, thr, side) = it.first;
  val = it.second;

  std::pair<float,int> key = std::make_pair(vov, thr);

  if( g_MPV_TOFHIR_LO_L.find(key) == g_MPV_TOFHIR_LO_L.end() )
  {
    g_MPV_TOFHIR_LO_L[key] = new TGraph();
    g_MPV_TOFHIR_LO_R[key] = new TGraph();
    counter_MPV_TOFHIR_LO_L[key] = 0;
    counter_MPV_TOFHIR_LO_R[key] = 0;
  }


  if(side == "L")
  {
    g_MPV_TOFHIR_LO_L[key]->SetPoint(counter_MPV_TOFHIR_LO_L[key], bar, val);
    counter_MPV_TOFHIR_LO_L[key]++;
  }
  else if(side == "R")
  {
    g_MPV_TOFHIR_LO_R[key]->SetPoint(counter_MPV_TOFHIR_LO_R[key], bar, val);
    counter_MPV_TOFHIR_LO_R[key]++;
  }
}


for(auto& it : g_MPV_TOFHIR_LO_L)
{
  float vov = it.first.first;
  int thr   = it.first.second;

  TGraph* gL = g_MPV_TOFHIR_LO_L[it.first];
  TGraph* gR = g_MPV_TOFHIR_LO_R[it.first];

  TCanvas* c = new TCanvas(Form("c_MPV_TOFHIR_LO_vov%.2f_th%d",vov,thr),Form("TOFHIR + LO factors Vov=%.2f th=%d",vov,thr),800,600);

  gL->SetMarkerStyle(20);
  gL->SetMarkerColor(kRed);
  gL->SetLineColor(kRed);
  gL->SetLineWidth(2);
  gL->GetYaxis()->SetRangeUser(325., 340.);

  gR->SetMarkerStyle(21);
  gR->SetMarkerColor(kBlue);
  gR->SetLineColor(kBlue);
  gR->SetLineWidth(2);

  gL->SetTitle(Form(";Bar index;MPV (LO+TOFHIR) [a.u.]"));

  gL->Draw("APL");
  gR->Draw("PL SAME");

  double mean_L, rms_L;
  double mean_R, rms_R;

  computeMeanRMS(gL, mean_L, rms_L);
  computeMeanRMS(gR, mean_R, rms_R);

  TLegend* leg = new TLegend(0.27,0.70,0.65,0.9);
  leg->SetHeader(Form("V_{OV} = %.2f V, vth = %i DAC", vov, thr), "C");
  leg->AddEntry(gL, Form("Side L (Mean=%.4f, RMS=%.4f)", mean_L, rms_L), "lp");
  leg->AddEntry(gR, Form("Side R (Mean=%.4f, RMS=%.4f)", mean_R, rms_R), "lp");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->Draw();

  c->SetGridy();
  c->SetTickx();
  c->SetTicky();

  c->Print(Form("%s/MPV_TOFHIR_LO_calib/c_MPV_TOFHIR_LO_vov%.2f_th%d.pdf",plotDir.c_str(),vov,thr));
  c->Print(Form("%s/MPV_TOFHIR_LO_calib/c_MPV_TOFHIR_LO_vov%.2f_th%d.png",plotDir.c_str(),vov,thr));

  delete leg;
  delete c;
}



//----------------------------
//   MPV TOFHIR colibration
//---------------------------- 
for(auto& it : MPV_TOFHIR_map)
{
  int bar;
  float vov;
  int thr;
  std::string side;
  double val;

  std::tie(bar, vov, thr, side) = it.first;
  val = it.second;

  std::pair<float,int> key = std::make_pair(vov, thr);

  // create TGraph if not exist
  if( g_MPV_TOFHIR_L.find(key) == g_MPV_TOFHIR_L.end() )
  {
    g_MPV_TOFHIR_L[key] = new TGraph();
    g_MPV_TOFHIR_R[key] = new TGraph();
    counter_MPV_TOFHIR_L[key] = 0;
    counter_MPV_TOFHIR_R[key] = 0;
  }


  if(side == "L")
  {
    g_MPV_TOFHIR_L[key]->SetPoint(counter_MPV_TOFHIR_L[key], bar, val);
    counter_MPV_TOFHIR_L[key]++;
  }
  else if(side == "R")
  {
    g_MPV_TOFHIR_R[key]->SetPoint(counter_MPV_TOFHIR_R[key], bar, val);
    counter_MPV_TOFHIR_R[key]++;
  }
}

for(auto& it : g_MPV_TOFHIR_L)
{
  float vov = it.first.first;
  int thr   = it.first.second;

  TGraph* gL = g_MPV_TOFHIR_L[it.first];
  TGraph* gR = g_MPV_TOFHIR_R[it.first];

  TCanvas* c = new TCanvas(Form("c_MPV_TOFHIR_vov%.2f_th%d",vov,thr),Form("TOFHIR factors Vov=%.2f th=%d",vov,thr),800,600);

  gL->SetMarkerStyle(20);
  gL->SetMarkerColor(kRed);
  gL->SetLineColor(kRed);
  gL->SetLineWidth(2);
  gL->GetYaxis()->SetRangeUser(285, 385);

  gR->SetMarkerStyle(21);
  gR->SetMarkerColor(kBlue);
  gR->SetLineColor(kBlue);
  gR->SetLineWidth(2);

  gL->SetTitle(Form(";Bar index;MPV (TOFHIR) [a.u.]"));

  gL->Draw("APL");
  gR->Draw("PL SAME");

  double mean_L, rms_L;
  double mean_R, rms_R;

  computeMeanRMS(gL, mean_L, rms_L);
  computeMeanRMS(gR, mean_R, rms_R);

  TLegend* leg = new TLegend(0.27,0.70,0.65,0.9);
  leg->SetHeader(Form("V_{OV} = %.2f V, vth = %i DAC", vov, thr), "C");
  leg->AddEntry(gL, Form("Side L (Mean=%.4f, RMS=%.4f)", mean_L, rms_L), "lp");
  leg->AddEntry(gR, Form("Side R (Mean=%.4f, RMS=%.4f)", mean_R, rms_R), "lp");
  leg->SetFillStyle(0);
  leg->SetBorderSize(0);
  leg->SetTextSize(0.03);
  leg->SetTextFont(42);
  leg->Draw();

  c->SetGridy();
  c->SetTickx();
  c->SetTicky();

  c->Print(Form("%s/MPV_TOFHIR_calib/c_MPV_TOFHIR_vov%.2f_th%d.pdf",plotDir.c_str(),vov,thr));
  c->Print(Form("%s/MPV_TOFHIR_calib/c_MPV_TOFHIR_vov%.2f_th%d.png",plotDir.c_str(),vov,thr));

  delete leg;
  delete c;
}



// -----------------------------------------------
//  SUMMARY MPV plots (RAW, LO, LO+TOFHIR, TOFHIR)
// -----------------------------------------------

for (auto &it : g_MPV_RAW_L) // loop su tutte le chiavi (vov, vth)
{
    float vov = it.first.first;
    int vth   = it.first.second;

    if (!g_MPV_LO_L.count(it.first)) continue;
    if (!g_MPV_TOFHIR_LO_L.count(it.first)) continue;
    if (!g_MPV_TOFHIR_L.count(it.first)) continue;

    
    // --- LEFT SIDE 
    TCanvas* c_L = new TCanvas(Form("c_MPV_summary_L_Vov%.2f_th%02d", vov, vth),Form("MPV summary L Vov%.2f th%02d", vov, vth),800, 600);

    auto g_raw  = g_MPV_RAW_L[it.first];
    auto g_lo   = g_MPV_LO_L[it.first];
    auto g_lo_t = g_MPV_TOFHIR_LO_L[it.first];
    auto g_t    = g_MPV_TOFHIR_L[it.first];

    g_raw->SetMarkerStyle(20);
    g_raw->SetMarkerColor(kBlue);
    g_raw->SetLineColor(kBlue);

    g_t->SetMarkerStyle(21);
    g_t->SetMarkerColor(kOrange);
    g_t->SetLineColor(kOrange);

    g_lo->SetMarkerStyle(22);
    g_lo->SetMarkerColor(kRed);
    g_lo->SetLineColor(kRed);

    g_lo_t->SetMarkerStyle(23);
    g_lo_t->SetMarkerColor(kGray+2);
    g_lo_t->SetLineColor(kGray+2);

    g_raw->SetTitle(Form(";Bar index;MPV"));
    g_raw->GetYaxis()->SetRangeUser(150, 550); 

    g_raw->Draw("APL");
    g_lo->Draw("PL SAME");
    g_lo_t->Draw("PL SAME");
    g_t->Draw("PL SAME");

    TLegend* leg_L = new TLegend(0.20,0.70,0.50,0.88);
    leg_L->SetFillStyle(0);
    leg_L->SetBorderSize(0);
    leg_L->SetTextSize(0.03);

    leg_L->SetHeader(Form("Side L  |  V_{OV}=%.2f V  |  vth=%d", vov, vth), "C");

    leg_L->AddEntry(g_raw,  "RAW", "lp");
    leg_L->AddEntry(g_t,    "TOFHIR calib", "lp");
    leg_L->AddEntry(g_lo,   "LO calib", "lp");
    leg_L->AddEntry(g_lo_t, "TOFHIR+LO calib", "lp");

    leg_L->Draw();

    c_L->SetGridy();
    c_L->SetTickx();
    c_L->SetTicky();

    c_L->Print(Form("%s/MPV_summary/MPV_summary_L_Vov%.2f_th%02d.pdf", plotDir.c_str(), vov, vth));
    c_L->Print(Form("%s/MPV_summary/MPV_summary_L_Vov%.2f_th%02d.png", plotDir.c_str(), vov, vth));

    delete c_L;

    
    
    // --- RIGHT SIDE
    if (!g_MPV_RAW_R.count(it.first)) continue;
    if (!g_MPV_LO_R.count(it.first)) continue;
    if (!g_MPV_TOFHIR_LO_R.count(it.first)) continue;
    if (!g_MPV_TOFHIR_R.count(it.first)) continue;

    TCanvas* c_R = new TCanvas(Form("c_MPV_summary_R_Vov%.2f_th%02d", vov, vth),Form("MPV summary R Vov%.2f th%02d", vov, vth),800, 600);

    auto g_raw_R  = g_MPV_RAW_R[it.first];
    auto g_lo_R   = g_MPV_LO_R[it.first];
    auto g_lo_t_R = g_MPV_TOFHIR_LO_R[it.first];
    auto g_t_R    = g_MPV_TOFHIR_R[it.first];

    g_raw_R->SetMarkerStyle(20);
    g_raw_R->SetMarkerColor(kBlue);
    g_raw_R->SetLineColor(kBlue);

    g_t_R->SetMarkerStyle(21);
    g_t_R->SetMarkerColor(kOrange);
    g_t_R->SetLineColor(kOrange);

    g_lo_R->SetMarkerStyle(22);
    g_lo_R->SetMarkerColor(kRed);
    g_lo_R->SetLineColor(kRed);

    g_lo_t_R->SetMarkerStyle(23);
    g_lo_t_R->SetMarkerColor(kGray+2);
    g_lo_t_R->SetLineColor(kGray+2);

    g_raw_R->SetTitle(Form(";Bar index;MPV"));
    g_raw_R->GetYaxis()->SetRangeUser(150, 550);

    g_raw_R->Draw("APL");
    g_lo_R->Draw("PL SAME");
    g_lo_t_R->Draw("PL SAME");
    g_t_R->Draw("PL SAME");

    TLegend* leg_R = new TLegend(0.20,0.70,0.50,0.88);
    leg_R->SetFillStyle(0);
    leg_R->SetBorderSize(0);
    leg_R->SetTextSize(0.03);

    leg_R->SetHeader(Form("Side R  |  V_{OV}=%.2f V  |  vth=%d", vov, vth), "C");

    leg_R->AddEntry(g_raw_R,  "RAW", "lp");
    leg_R->AddEntry(g_lo_R,   "LO", "lp");
    leg_R->AddEntry(g_lo_t_R, "LO+TOFHIR", "lp");
    leg_R->AddEntry(g_t_R,    "TOFHIR", "lp");

    leg_R->Draw();

    c_R->SetGridy();
    c_R->SetTickx();
    c_R->SetTicky();

    c_R->Print(Form("%s/MPV_summary/MPV_summary_R_Vov%.2f_th%02d.pdf", plotDir.c_str(), vov, vth));
    c_R->Print(Form("%s/MPV_summary/MPV_summary_R_Vov%.2f_th%02d.png", plotDir.c_str(), vov, vth));

    delete c_R;
}



  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
