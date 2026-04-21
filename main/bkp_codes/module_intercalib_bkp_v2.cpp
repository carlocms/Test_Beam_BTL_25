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




//------------------------------------------------------------------------
//---- find energy bins
//------------------------------------------------------------------------
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

//------------------------------------------------------------------------
//---- Find fit range profile
//------------------------------------------------------------------------
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

//------------------------------------------------------------------------
// Compute the RMS and the mean from a TGraph plot
//------------------------------------------------------------------------
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



//------------------------------------------------------------------------
// --- For Landau Fit of the energy histograms 
//------------------------------------------------------------------------
struct FitResult {
    TF1* func;
    double mpv;
    double chi2_ndf;
};

FitResult DoLandauFit(TH1* histo,double minE,const std::string& name){
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
    TF1* f = new TF1(name.c_str(),"[0]*TMath::Landau(x,[1],[2])",0., 1000.);

    // --- initial parameters
    f->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10.,max,0.1 * max);

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

        f->SetParameters(histo->Integral(histo->GetMaximumBin(), histo->GetNbinsX())/10.,f->GetParameter(1),0.1 * f->GetParameter(1));
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




//------------------------------------------------------------------------
// --- Find range Y axis TGraph:
//------------------------------------------------------------------------
void SetAutoYRange(TGraph* g, double frac = 0.05)
{
    if (!g || g->GetN() == 0) return;

    double ymin = 1e9;
    double ymax = -1e9;

    for (int i = 0; i < g->GetN(); ++i)
    {
        double x, y;
        g->GetPoint(i, x, y);

        if (y < ymin) ymin = y;
        if (y > ymax) ymax = y;
    }

    if (ymax == ymin)
    {
        double delta = (ymin != 0) ? 0.05 * std::abs(ymin) : 0.01;
        ymin -= delta;
        ymax += delta;
    }
    else
    {
        double margin = frac * (ymax - ymin);
        ymin -= margin;
        ymax += margin;
    }

    // --- apply range
    g->GetYaxis()->SetRangeUser(ymin, ymax);
}



//------------------------------------------------------------------------
// --- helper function TGraph style:
//------------------------------------------------------------------------
void StyleGraph(TGraph* g, int color, int marker)
{
    g->SetMarkerStyle(marker);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
}



//------------------------------------------------------------------------
//--- helper function for Chi2 and MPV histos:
//------------------------------------------------------------------------
void FillAndPlotHistos(
    const std::map<std::tuple<int,float,int,std::string>, double>& input_map,
    const std::string& label,
    const std::string& varName,
    int nBins,
    double xmin,
    double xmax,
    double minCut,
    double maxCut,
    const std::string& outDir
)
{
    std::set<std::pair<float,int>> vov_thr_set;

    // --- find all (vov,thr)
    for (const auto& it : input_map)
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

        TH1F* h = new TH1F(Form("h_%s_%s_Vov%.2f_th%02d", varName.c_str(), label.c_str(), vov, thr),Form(";%s;Entries", varName.c_str()),nBins, xmin, xmax);

        // --- fill
        for (const auto& it : input_map)
        {
            int bar;
            float vov_i;
            int thr_i;
            std::string side;

            std::tie(bar, vov_i, thr_i, side) = it.first;

            if (vov_i != vov || thr_i != thr) continue;

            double val = it.second;

            if (val < minCut || val > maxCut) continue;

            h->Fill(val);
        }

        // --- canvas
        TCanvas* c = new TCanvas(Form("c_%s_%s_Vov%.2f_th%02d", varName.c_str(), label.c_str(), vov, thr),"", 800, 600);

        h->SetLineWidth(2);
        h->SetLineColor(kBlue);
        h->Draw("hist");

        double mean = h->GetMean();
        double rms  = h->GetRMS();

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);

        latex.DrawLatex(0.60,0.85,Form("Mean = %.3f", mean));
        latex.DrawLatex(0.60,0.80,Form("RMS  = %.3f", rms));

        latex.DrawLatex(0.20,0.85,Form("V_{OV} = %.2f V, thr = %d", vov, thr));

	c->Print(Form("%s/c_%s_Vov%.2f_th%02d.pdf",outDir.c_str(), label.c_str(), vov, thr));
        c->Print(Form("%s/c_%s_Vov%.2f_th%02d.png",outDir.c_str(), label.c_str(), vov, thr));
        delete c;
    }
}



//-------------------------------------------------------------------------------------------
// --- helper function TGaph Left and right same:
//------------------------------------------------------------------------
void BuildGraphs(
    const std::map<std::tuple<int,float,int,std::string>, double>& input_map,
    std::map<std::pair<float,int>, TGraph*>& gL_map,
    std::map<std::pair<float,int>, TGraph*>& gR_map
)
{
    std::map<std::pair<float,int>, int> counterL, counterR;

    for(const auto& it : input_map)
    {
        int bar;
        float vov;
        int thr;
        std::string side;
        double val;

        std::tie(bar, vov, thr, side) = it.first;
        val = it.second;

        std::pair<float,int> key = std::make_pair(vov, thr);

        if(gL_map.find(key) == gL_map.end())
        {
            gL_map[key] = new TGraph();
            gR_map[key] = new TGraph();
            counterL[key] = 0;
            counterR[key] = 0;
        }

        if(side == "L")
            gL_map[key]->SetPoint(counterL[key]++, bar, val);
        else if(side == "R")
            gR_map[key]->SetPoint(counterR[key]++, bar, val);
    }
}

//------------------------------------------------------------------------
// --- helper function TGaph Left and right same:
//------------------------------------------------------------------------
void PlotGraphs(
    const std::map<std::pair<float,int>, TGraph*>& gL_map,
    const std::map<std::pair<float,int>, TGraph*>& gR_map,
    const std::string& label,
    const std::string& yTitle,
    const std::string& outDir,
    double ymin = -1,
    double ymax = -1
)
{
    for(const auto& it : gL_map)
    {
        float vov = it.first.first;
        int thr   = it.first.second;

        TGraph* gL = gL_map.at(it.first);
        TGraph* gR = gR_map.at(it.first);

        TCanvas* c = new TCanvas(
            Form("c_%s_vov%.2f_th%d", label.c_str(), vov, thr),
            "", 800, 600
        );

        StyleGraph(gL, kRed, 20);
        StyleGraph(gR, kBlue, 21);

        // range
        if (ymin < ymax)
            gL->GetYaxis()->SetRangeUser(ymin, ymax);
        else
            SetAutoYRange(gL);

        gL->SetTitle(Form(";Bar index; %s", yTitle.c_str()));

        gL->Draw("APL");
        gR->Draw("PL SAME");

        // stats
        double mean_L, rms_L;
        double mean_R, rms_R;

        computeMeanRMS(gL, mean_L, rms_L);
        computeMeanRMS(gR, mean_R, rms_R);

        TLegend* leg = new TLegend(0.27,0.70,0.65,0.9);
        leg->SetHeader(Form("V_{OV} = %.2f V, vth = %d DAC", vov, thr), "C");

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

        c->Print(Form("%s/%s_vov%.2f_th%d.pdf",
                      outDir.c_str(), label.c_str(), vov, thr));
        c->Print(Form("%s/%s_vov%.2f_th%d.png",
                      outDir.c_str(), label.c_str(), vov, thr));

        delete leg;
        delete c;
    }
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
    std::cout << ">>> module_intercalibration::usage:   " << argv[0] << " configFile.cfg" << std::endl;
    return -1;
  }

  //--- parse the config file
  CfgManager opts;
  opts.ParseConfigFile(argv[1]);

  int debugMode = 0;
  if( argc > 2 ) debugMode = atoi(argv[2]);


  bool degub_eneIntercalib = false;

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

  // --- Chi2 and MPV check historgrams:
  system(Form("mkdir -p %s/Chi2_plots/", plotDir.c_str()));
  system(Form("mkdir -p %s/MPV_plots/", plotDir.c_str()));

  // --- From calibration factors plots:
  system(Form("mkdir -p %s/LO_calibration_factors/",plotDir.c_str()));
  system(Form("mkdir -p %s/RESIDUAL_TOFHIR_calibration_factors/",plotDir.c_str())); 
  system(Form("mkdir -p %s/TOFHIR_calibration_factors/",plotDir.c_str()));
  system(Form("mkdir -p %s/COMBINED_calibration_factors/",plotDir.c_str()));

  // --- From Loop 3:
  system(Form("mkdir -p %s/Loop3_energy_RESIDUAL_TOFHIR_calib/Left/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_energy_RESIDUAL_TOFHIR_calib/Right/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_energy_COMBINED_calib/Left/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_energy_COMBINED_calib/Right/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_energy_TOFHIR_calib/Left/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_energy_TOFHIR_calib/Right/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_scatter_energy_RESIDUAL_TOFHIR_calib_L_vs_R/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_scatter_energy_COMBINED_calib_L_vs_R/",plotDir.c_str()));
  system(Form("mkdir -p %s/Loop3_scatter_energy_TOFHIR_calib_L_vs_R/",plotDir.c_str()));

  // --- From MPV calibration:
  system(Form("mkdir -p %s/MPV_RAW/",plotDir.c_str()));
  system(Form("mkdir -p %s/MPV_LO_calib/",plotDir.c_str()));
  system(Form("mkdir -p %s/MPV_COMBINED_calib/",plotDir.c_str()));
  system(Form("mkdir -p %s/MPV_RESIDUAL_TOFHIR_calib/",plotDir.c_str()));
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

	if (degub_eneIntercalib){
        std::cout << "LO factor loaded: bar " << bar << " side " << side << " value " << value << std::endl;
	}
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
  std::map< std::pair<float,int>, TGraph* > g_RESIDUAL_TOFHIR_L;
  std::map< std::pair<float,int>, TGraph* > g_RESIDUAL_TOFHIR_R;
  std::map< std::pair<float,int>, TGraph* > g_TOFHIR_L;
  std::map< std::pair<float,int>, TGraph* > g_TOFHIR_R;
  std::map< std::pair<float,int>, TGraph* > g_COMBINED_calib_factors_L;
  std::map< std::pair<float,int>, TGraph* > g_COMBINED_calib_factors_R;
  std::map< std::pair<float,int>, int > counter_TOFHIR_L;
  std::map< std::pair<float,int>, int > counter_TOFHIR_R;
  std::map< std::pair<float,int>, int > counter_RESIDUAL_TOFHIR_L;
  std::map< std::pair<float,int>, int > counter_RESIDUAL_TOFHIR_R;
  std::map< std::pair<float,int>, int > counter_COMBINED_L;
  std::map< std::pair<float,int>, int > counter_COMBINED_R;

  //Loop3:
  std::map<double,TH1F*> h1_energy_COMBINED_calib_L;
  std::map<double,TH1F*> h1_energy_COMBINED_calib_R;
  std::map<double,TH1F*> h1_energy_RESIDUAL_TOFHIR_calib_L;
  std::map<double,TH1F*> h1_energy_RESIDUAL_TOFHIR_calib_R;
  std::map<double,TH1F*> h1_energy_TOFHIR_calib_L;
  std::map<double,TH1F*> h1_energy_TOFHIR_calib_R;
  std::map<double,TH2F*> h2_energy_RESIDUAL_TOFHIR_calib_L_vs_R;
  std::map<double,TH2F*> h2_energy_TOFHIR_calib_L_vs_R;
  std::map<double,TH2F*> h2_energy_COMBINED_calib_L_vs_R;

  // MPV calibration:
  std::map< std::pair<float,int>, TGraph* > g_MPV_RESIDUAL_TOFHIR_L;
  std::map< std::pair<float,int>, TGraph* > g_MPV_RESIDUAL_TOFHIR_R;
  std::map< std::pair<float,int>, TGraph* > g_MPV_COMBINED_L;
  std::map< std::pair<float,int>, TGraph* > g_MPV_COMBINED_R;
  std::map< std::pair<float,int>, TGraph* > g_MPV_TOFHIR_L;
  std::map< std::pair<float,int>, TGraph* > g_MPV_TOFHIR_R;
  std::map< std::pair<float,int>, TGraph* > g_MPV_LO_L;
  std::map< std::pair<float,int>, TGraph* > g_MPV_LO_R;
  std::map< std::pair<float,int>, TGraph* > g_MPV_RAW_L;
  std::map< std::pair<float,int>, TGraph* > g_MPV_RAW_R;
  std::map< std::pair<float,int>, int > counter_MPV_RESIDUAL_TOFHIR_L;
  std::map< std::pair<float,int>, int > counter_MPV_RESIDUAL_TOFHIR_R;
  std::map< std::pair<float,int>, int > counter_MPV_COMBINED_L;
  std::map< std::pair<float,int>, int > counter_MPV_COMBINED_R;
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
    
std::map< std::tuple<int,float,int,std::string>, double > MPV_LO_map;           // map [bar, vov, vth, side] -> MPV values
std::map< std::pair<float,int>, std::vector<double> > MPV_values_map;        // map [vov, vth] -> vector(MPVs)  to achieve <MPV> all the channels
std::map< std::tuple<float,int,std::string>, std::vector<double> > MPV_values_map_side; //caso separated side
std::map< std::tuple<int,float,int,std::string>, double > Chi2_values_LO_map; // mapping of the chi2 of the energy historgrams with LO correction


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

	      // --- Landau Fit
	      auto fitRes_LOcalib_L = DoLandauFit(histo,minE[std::make_pair(iBar, Vov)],Form("f_landau_LO_bar%02d_%s", iBar, labelLR_energyBin.c_str()));

	      // draw
	      if (fitRes_LOcalib_L.func)fitRes_LOcalib_L.func->Draw("same");

	      // fill maps
	      if (fitRes_LOcalib_L.mpv > 0)
	      {
		      MPV_LO_map[std::make_tuple(iBar,Vov,vth1,"L")] = fitRes_LOcalib_L.mpv;
		      MPV_values_map[std::make_pair(Vov,int(vth1))].push_back(fitRes_LOcalib_L.mpv);
		      MPV_values_map_side[std::make_tuple(Vov,vth1,"L")].push_back(fitRes_LOcalib_L.mpv);
		      Chi2_values_LO_map[std::make_tuple(iBar,Vov,vth1,"L")] = fitRes_LOcalib_L.chi2_ndf;
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

	      // --- Landau Fit
              auto fitRes_LOcalib_R = DoLandauFit(histo,minE[std::make_pair(iBar, Vov)],Form("f_landau_LO_bar%02d_%s", iBar, labelLR_energyBin.c_str()));

              // draw
              if (fitRes_LOcalib_R.func)fitRes_LOcalib_R.func->Draw("same");

              // fill maps
              if (fitRes_LOcalib_R.mpv > 0)
              {
                      MPV_LO_map[std::make_tuple(iBar,Vov,vth1,"R")] = fitRes_LOcalib_R.mpv;
                      MPV_values_map[std::make_pair(Vov,int(vth1))].push_back(fitRes_LOcalib_R.mpv);
                      MPV_values_map_side[std::make_tuple(Vov,vth1,"R")].push_back(fitRes_LOcalib_R.mpv);
                      Chi2_values_LO_map[std::make_tuple(iBar,Vov,vth1,"R")] = fitRes_LOcalib_R.chi2_ndf;
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

std::map< std::pair<float,int>, double > MPV_mean_map;                       // map [vov, vth] -> <MPV> all the channels 
std::map< std::tuple<float,int,std::string>, double > MPV_mean_map_side;

std::map< std::tuple<int,float,int,std::string>, double > RESIDUAL_TOFHIR_factors; // map [bar, vov, vth, side] -> <MPV> / MPV  
std::map< std::tuple<int,float,int,std::string>, double > TOFHIR_factors;      // map [bar, vov, vth, side] -> TOFHIR_LO_factor / TOFHIR_factor
std::map< std::tuple<int,float,int,std::string>, double > COMBINED_factors;



//--------------------------------------------------------------------
// Averaging MPV values <-> from energy distribution + LO equalization
//--------------------------------------------------------------------

if (LO_mode == "global") // averaging over all the channels
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
else if (LO_mode == "separated_side") // averaging over the channels of the same side (L, R)
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




//---------------------------------------------
// Determination of the RESIDUAL_TOFHIR_factors
//
// Residual term in the energy intercalibraton once filter out the LO component
//---------------------------------------------

for (auto &it : MPV_LO_map)
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

    double calib_factor = mean_mpv / mpv;

    RESIDUAL_TOFHIR_factors[it.first] = calib_factor;

    if (degub_eneIntercalib){
    std::cout << "[RESIDUAL TOFHIR CALIB]: " << "| bar " << bar  << "| thr " << thr << "| side " << side << " --> " << calib_factor << std::endl;
    }
}




//------------------------------------------------
// Determination of the COMBINED calibration factors
//------------------------------------------------

for (auto &it : RESIDUAL_TOFHIR_factors)
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

    double combined_calib = TOFHIR_LO_calib * LO_factor_value;

    COMBINED_factors[it.first] = combined_calib;

    if (degub_eneIntercalib){
	    std::cout << "[COMBINED CALIB]: "<< "| bar " << bar << "| thr " << thr << "| side " << side << " --> " << combined_calib << std::endl;
    }
}



//------------------------------------
// Determination of the TOFHIR_factors
//------------------------------------

for (auto &it : RESIDUAL_TOFHIR_factors)
{
    int bar;
    float vov;
    int thr;
    std::string side;

    std::tie(bar, vov, thr, side) = it.first;
    double RESIDUAL_TOFHIR_calib = it.second;

    if (RESIDUAL_TOFHIR_calib <= 0) continue;

    auto key_LO = std::make_pair(bar, side);
    if (!LO_factors.count(key_LO)) continue;

    double LO_factor_value = LO_factors[key_LO];

    double TOFHIR_calib = RESIDUAL_TOFHIR_calib / LO_factor_value;

    TOFHIR_factors[it.first] = TOFHIR_calib;

    if (degub_eneIntercalib){
    std::cout << "[TOFHIR CALIB]: " <<"| bar " << bar << "| thr " << thr << "| side " << side << " --> " << TOFHIR_calib << std::endl;
    }
}

std::cout << std::endl;
std::cout << "Determination Calibration factors, done!" << std::endl;


//------------------------------------
//   Txt output files calib factors
//------------------------------------

std::set< std::pair<float,int> > vov_thr_set;

for (auto &it : RESIDUAL_TOFHIR_factors)
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

    // --- FILE RESIDUAL TOFHIR:
    std::string outFileName_RESIDUAL_TOFHIR = Form("%s/RESIDUAL_TOFHIR_calibration_factors_Vov%.2f_th%02d_%s.txt",outputDirTxt.c_str(), vov, thr, LO_mode.c_str());

    std::ofstream outFile_RESIDUAL_TOFHIR(outFileName_RESIDUAL_TOFHIR);

    if (!outFile_RESIDUAL_TOFHIR.is_open())
    {
        std::cerr << "ERROR: cannot open output file " << outFileName_RESIDUAL_TOFHIR << std::endl;
    }
    else
    {
        outFile_RESIDUAL_TOFHIR << "# mode: " << LO_mode << "\n";
        outFile_RESIDUAL_TOFHIR << "# bar side calib\n";

        std::map<int, std::map<std::string,double>> ordered_RESIDUAL_TOFHIR;

        for (auto &it : RESIDUAL_TOFHIR_factors)
        {
            int bar;
            float vov_i;
            int thr_i;
            std::string side;

            std::tie(bar, vov_i, thr_i, side) = it.first;

            if (vov_i != vov || thr_i != thr) continue;

            ordered_RESIDUAL_TOFHIR[bar][side] = it.second;
        }

        for (auto &b : ordered_RESIDUAL_TOFHIR)
        {
            for (auto &s : b.second)
            {
                outFile_RESIDUAL_TOFHIR << b.first << " " << s.first << " " << s.second << "\n";
            }
        }

        outFile_RESIDUAL_TOFHIR.close();
        //std::cout << "[WRITE] RESIDUAL TOFHIR file saved: " << outFileName_RESIDUAL_TOFHIR << std::endl;
    }

    // --- FILE TOFHIR:
    std::string outFileName_TOFHIR = Form("%s/TOFHIR_calibration_factors_Vov%.2f_th%02d_%s.txt",outputDirTxt.c_str(), vov, thr, LO_mode.c_str());
    std::ofstream outFile_TOFHIR(outFileName_TOFHIR);

    if (!outFile_TOFHIR.is_open())
    {
        std::cerr << "ERROR: cannot open output file " << outFileName_TOFHIR << std::endl;
    }
    else
    {
        outFile_TOFHIR << "# mode: " << LO_mode << "\n";
        outFile_TOFHIR << "# bar side calib\n";

        std::map<int, std::map<std::string,double>> ordered_TOFHIR;

        for (auto &it : TOFHIR_factors)
        {
            int bar;
            float vov_i;
            int thr_i;
            std::string side;

            std::tie(bar, vov_i, thr_i, side) = it.first;

            if (vov_i != vov || thr_i != thr) continue;

            ordered_TOFHIR[bar][side] = it.second;
        }

        for (auto &b : ordered_TOFHIR)
        {
            for (auto &s : b.second)
            {
                outFile_TOFHIR << b.first << " " << s.first << " " << s.second << "\n";
            }
        }

        outFile_TOFHIR.close();
        //std::cout << "[WRITE] TOFHIR file saved: " << outFileName_TOFHIR << std::endl;
    }


    // --- FILE COMBINED factors:
    std::string outFileName_COMBINED = Form("%s/COMBINED_calibration_factors_Vov%.2f_th%02d_%s.txt",outputDirTxt.c_str(), vov, thr, LO_mode.c_str());

    std::ofstream outFile_COMBINED(outFileName_COMBINED);

    if (!outFile_COMBINED.is_open())
    {
        std::cerr << "ERROR: cannot open output file " << outFileName_COMBINED << std::endl;
    }
    else
    {
        outFile_COMBINED << "# mode: " << LO_mode << "\n";
        outFile_COMBINED << "# bar side calib\n";

        std::map<int, std::map<std::string,double>> ordered_COMBINED;

        for (auto &it : COMBINED_factors)
        {
            int bar;
            float vov_i;
            int thr_i;
            std::string side;

            std::tie(bar, vov_i, thr_i, side) = it.first;

            if (vov_i != vov || thr_i != thr) continue;

            ordered_COMBINED[bar][side] = it.second;
        }

        for (auto &b : ordered_COMBINED)
        {
            for (auto &s : b.second)
            {
                outFile_COMBINED << b.first << " " << s.first << " " << s.second << "\n";
            }
        }

        outFile_COMBINED.close();
        //std::cout << "[WRITE] COMBINED file saved: " << outFileName_COMBINED << std::endl;
    }
}
std::cout << std::endl;
std::cout << "[WRITE] Calibration factors in txt Files, done!" << std::endl;



//////////////////////////////////////////
//   Plot Callibration factors vs bar   //
//////////////////////////////////////////


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

StyleGraph(g_LO_L, kRed, 20);
g_LO_L->GetYaxis()->SetRangeUser(0.93, 1.09);
StyleGraph(g_LO_R, kBlue, 21);

g_LO_L->Draw("APL"); 
g_LO_R->Draw("PL SAME");

// Mean and RMS values
double mean_L, rms_L;
double mean_R, rms_R;

computeMeanRMS(g_LO_L, mean_L, rms_L);
computeMeanRMS(g_LO_R, mean_R, rms_R);

TLegend* leg = new TLegend(0.20,0.70,0.65,0.9);
leg->AddEntry(g_LO_L, Form("Side L (Mean=%.4f, RMS=%.4f)", mean_L, rms_L), "lp");
leg->AddEntry(g_LO_R, Form("Side R (Mean=%.4f, RMS=%.4f)", mean_R, rms_R), "lp");
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
//  Plot of the RESIDUAL factors
//-------------------------------
/*
for(auto& it : RESIDUAL_TOFHIR_factors)
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
  if( g_RESIDUAL_TOFHIR_L.find(key) == g_RESIDUAL_TOFHIR_L.end() )
  {
    g_RESIDUAL_TOFHIR_L[key] = new TGraph();
    g_RESIDUAL_TOFHIR_R[key] = new TGraph();
    counter_RESIDUAL_TOFHIR_L[key] = 0;
    counter_RESIDUAL_TOFHIR_R[key] = 0;
  }

  if(side == "L")
  {
    g_RESIDUAL_TOFHIR_L[key]->SetPoint(counter_RESIDUAL_TOFHIR_L[key], bar, val);
    counter_RESIDUAL_TOFHIR_L[key]++;
  }
  else if(side == "R")
  {
    g_RESIDUAL_TOFHIR_R[key]->SetPoint(counter_RESIDUAL_TOFHIR_R[key], bar, val);
    counter_RESIDUAL_TOFHIR_R[key]++;
  }
}


for(auto& it : g_RESIDUAL_TOFHIR_L)
{
  float vov = it.first.first;
  int thr   = it.first.second;

  TGraph* gL = g_RESIDUAL_TOFHIR_L[it.first];
  TGraph* gR = g_RESIDUAL_TOFHIR_R[it.first];

  TCanvas* c = new TCanvas(Form("c_RESIDUAL_TOFHIR_vov%.2f_th%d",vov,thr),Form("RESIDUAL TOFHIR factors Vov=%.2f th=%d",vov,thr),800,600);

  StyleGraph(gL, kRed, 20);
  gL->GetYaxis()->SetRangeUser(0.5, 1.5);//0.75, 1.35
  StyleGraph(gR, kBlue, 21);
  gL->SetTitle(Form(";Bar index;RESIDUAL TOFHIR calibration"));

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

  c->Print(Form("%s/RESIDUAL_TOFHIR_calibration_factors/RESIDUAL_TOFHIR_factors_vov%.2f_th%d.pdf",plotDir.c_str(),vov,thr));
  c->Print(Form("%s/RESIDUAL_TOFHIR_calibration_factors/RESIDUAL_TOFHIR_factors_vov%.2f_th%d.png",plotDir.c_str(),vov,thr));

  delete leg;
  delete c;
}
*/


BuildGraphs(RESIDUAL_TOFHIR_factors,g_RESIDUAL_TOFHIR_L,g_RESIDUAL_TOFHIR_R);
PlotGraphs(g_RESIDUAL_TOFHIR_L,g_RESIDUAL_TOFHIR_R,"RESIDUAL_TOFHIR","Residual TOFHIR calibration",plotDir + "/RESIDUAL_TOFHIR_calibration_factors",0.5, 1.5);


//-------------------------------
//  Plot of the COMBINED factors
//-------------------------------
BuildGraphs(COMBINED_factors,g_COMBINED_calib_factors_L,g_COMBINED_calib_factors_R);
PlotGraphs(g_COMBINED_calib_factors_L,g_COMBINED_calib_factors_R,"COMBINED_factors","Combined calibration",plotDir + "/COMBINED_calibration_factors",0.5, 1.5);

/*
for(auto& it : COMBINED_factors)
{
  int bar;
  float vov;
  int thr;
  std::string side;
  double val;

  std::tie(bar, vov, thr, side) = it.first;
  val = it.second;

  std::pair<float,int> key = std::make_pair(vov, thr);

  if( g_COMBINED_calib_factors_L.find(key) == g_COMBINED_calib_factors_L.end() )
  { 
    g_COMBINED_calib_factors_L[key] = new TGraph();
    g_COMBINED_calib_factors_R[key] = new TGraph();
    counter_COMBINED_L[key] = 0;
    counter_COMBINED_R[key] = 0;
  }

  if(side == "L")
  {
    g_COMBINED_calib_factors_L[key]->SetPoint(counter_COMBINED_L[key], bar, val);
    counter_COMBINED_L[key]++;
  }
  else if(side == "R")
  {
    g_COMBINED_calib_factors_R[key]->SetPoint(counter_COMBINED_R[key], bar, val);
    counter_COMBINED_R[key]++;
  }
}


for(auto& it : g_COMBINED_calib_factors_L)
{
  float vov = it.first.first;
  int thr   = it.first.second;

  TGraph* gL = g_COMBINED_calib_factors_L[it.first];
  TGraph* gR = g_COMBINED_calib_factors_R[it.first];

  TCanvas* c = new TCanvas(Form("c_Combined_vov%.2f_th%d",vov,thr),Form("Combined calibration factors Vov=%.2f th=%d",vov,thr),800,600);

  StyleGraph(gL, kRed, 20);
  gL->GetYaxis()->SetRangeUser(0.5, 1.5);//0.75, 1.35
  StyleGraph(gR, kBlue, 21);
  gL->SetTitle(Form(";Bar index;Combined calibration factors"));

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

  c->Print(Form("%s/COMBINED_calibration_factors/COMBINED_factors_vov%.2f_th%d.pdf",plotDir.c_str(),vov,thr));
  c->Print(Form("%s/COMBINED_calibration_factors/COMBINED_factors_vov%.2f_th%d.png",plotDir.c_str(),vov,thr));

  delete leg;
  delete c;
}
*/



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

  StyleGraph(gL, kRed, 20);
  //gL->GetYaxis()->SetRangeUser(0.75, 1.35);
  SetAutoYRange(gL);
  StyleGraph(gR, kBlue, 21);
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

std::cout << std::endl;
std::cout << "Calibration factors plots done!" << std::endl;







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
            std::cout << ">>> ENERGY INTERCALIBRATION loop 3: " << mapIt.first << " reading entry " << entry << " / " << nEntries << " (" << 100.*entry/nEntries << "%)" << "\r" << std::flush;
            //TrackProcess(cpu, mem, vsz, rss);
          }
          mapIt.second -> GetEntry(entry);

          bool barFound = std::find(barList.begin(), barList.end(), anEvent->barID) != barList.end() ;
          if (!barFound) continue;

          int index1( (10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );

          if( !accept[index1][entry] ) continue;

          int energyBinAverage = FindBin(0.5*(anEvent->energyL+anEvent->energyR),ranges["L-R"][index1])+1;

          double index2( (10000000*energyBinAverage+10000*int(anEvent->Vov*100.)) + (100*anEvent->vth1) + anEvent->barID );


          // --- Key to access LO factors
	  auto keyLO_L = std::make_pair(anEvent->barID, "L");
          auto keyLO_R = std::make_pair(anEvent->barID, "R");

          if( LO_factors.find(keyLO_L) == LO_factors.end() || LO_factors.find(keyLO_R) == LO_factors.end() )
          {
                  std::cerr << "[ERROR] Missing LO factor for bar "<< anEvent->barID << std::endl;
                  continue;
          }
	  

	  // --- Key to access RESIDUAL_TOFHIR factors
	  auto key_RESIDUAL_TOFHIR_L = std::make_tuple(anEvent->barID, anEvent->Vov, anEvent->vth1, "L");
          auto key_RESIDUAL_TOFHIR_R = std::make_tuple(anEvent->barID, anEvent->Vov, anEvent->vth1, "R");
 
          if( RESIDUAL_TOFHIR_factors.find(key_RESIDUAL_TOFHIR_L) == RESIDUAL_TOFHIR_factors.end() || RESIDUAL_TOFHIR_factors.find(key_RESIDUAL_TOFHIR_R) == RESIDUAL_TOFHIR_factors.end() )
          {
                  std::cerr << "[ERROR] Missing RESIDUAL TOFHIR factor for bar "<< anEvent->barID << std::endl;
                  continue;
          }
	  
	  
	  // --- Key to access RESIDUAL_TOFHIR factors
	  auto key_TOFHIR_L = std::make_tuple(anEvent->barID, anEvent->Vov, anEvent->vth1, "L");
          auto key_TOFHIR_R = std::make_tuple(anEvent->barID, anEvent->Vov, anEvent->vth1, "R");

          if( TOFHIR_factors.find(key_TOFHIR_L) == TOFHIR_factors.end() || TOFHIR_factors.find(key_TOFHIR_R) == TOFHIR_factors.end() )
          {
                  std::cerr << "[ERROR] Missing TOFHIR factor for bar "<< anEvent->barID << std::endl;
                  continue;
          }


	  float energy_COMBINED_calib_L = (anEvent->energyL * LO_factors[keyLO_L]) * RESIDUAL_TOFHIR_factors[key_RESIDUAL_TOFHIR_L];
          float energy_COMBINED_calib_R = (anEvent->energyR * LO_factors[keyLO_R]) * RESIDUAL_TOFHIR_factors[key_RESIDUAL_TOFHIR_R];
  
	  float energy_RESIDUAL_TOFHIR_calib_L = (anEvent->energyL) * RESIDUAL_TOFHIR_factors[key_RESIDUAL_TOFHIR_L];
	  float energy_RESIDUAL_TOFHIR_calib_R = (anEvent->energyR) * RESIDUAL_TOFHIR_factors[key_RESIDUAL_TOFHIR_R];

	  float energy_TOFHIR_calib_L = (anEvent->energyL) * TOFHIR_factors[key_TOFHIR_L];
          float energy_TOFHIR_calib_R = (anEvent->energyR) * TOFHIR_factors[key_TOFHIR_R];

          if( h1_energy_RESIDUAL_TOFHIR_calib_L[index2] == NULL )
            {
              std::string labelLR_energyBin(Form("bar%02d_Vov%.2f_th%02d_energyBin%02d",anEvent->barID,anEvent->Vov,anEvent->vth1,energyBinAverage));
	 
	      h1_energy_RESIDUAL_TOFHIR_calib_L[index2] = new TH1F(Form("h1_energy_RESIDUAL_TOFHIR_calib_L_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024.);
              h1_energy_RESIDUAL_TOFHIR_calib_R[index2] = new TH1F(Form("h1_energy_RESIDUAL_TOFHIR_calib_R_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024.);
	      h2_energy_RESIDUAL_TOFHIR_calib_L_vs_R[index2] = new TH2F(Form("h2_energy_RESIDUAL_TOFHIR_calib_L_vs_R_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024., 512, 0., 1024.);
 
	                    
	      
	      h1_energy_COMBINED_calib_L[index2] = new TH1F(Form("h1_energy_COMBINED_calib_L_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024.);
              h1_energy_COMBINED_calib_R[index2] = new TH1F(Form("h1_energy_COMBINED_calib_R_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024.);
              h2_energy_COMBINED_calib_L_vs_R[index2] = new TH2F(Form("h2_energy_COMBINED_calib_L_vs_R_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024., 512, 0., 1024.);

	      h1_energy_TOFHIR_calib_L[index2] = new TH1F(Form("h1_energy_TOFHIR_calib_L_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024.);
              h1_energy_TOFHIR_calib_R[index2] = new TH1F(Form("h1_energy_TOFHIR_calib_R_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024.);
	      h2_energy_TOFHIR_calib_L_vs_R[index2] = new TH2F(Form("h2_energy_TOFHIR_calib_L_vs_R_%s",labelLR_energyBin.c_str()),"", 512, 0., 1024., 512, 0., 1024.);
            }

          
	  h1_energy_RESIDUAL_TOFHIR_calib_L[index2] -> Fill(energy_RESIDUAL_TOFHIR_calib_L);
	  h1_energy_RESIDUAL_TOFHIR_calib_R[index2] -> Fill(energy_RESIDUAL_TOFHIR_calib_R);
	  h2_energy_RESIDUAL_TOFHIR_calib_L_vs_R[index2] -> Fill(energy_RESIDUAL_TOFHIR_calib_L, energy_RESIDUAL_TOFHIR_calib_R);

	  h1_energy_COMBINED_calib_L[index2] -> Fill(energy_COMBINED_calib_L);
          h1_energy_COMBINED_calib_R[index2] -> Fill(energy_COMBINED_calib_R);
          h2_energy_COMBINED_calib_L_vs_R[index2] -> Fill(energy_COMBINED_calib_L, energy_COMBINED_calib_R);

	  h1_energy_TOFHIR_calib_L[index2] -> Fill(energy_TOFHIR_calib_L);
          h1_energy_TOFHIR_calib_R[index2] -> Fill(energy_TOFHIR_calib_R);
	  h2_energy_TOFHIR_calib_L_vs_R[index2] -> Fill(energy_TOFHIR_calib_L, energy_TOFHIR_calib_R);

          }
      std::cout << std::endl;
    }






  /////////////////////////
  //   draw 3nd plots   //
  ////////////////////////

  std::map<double,TF1*> fitFunc_energy_RESIDUAL_TOFHIR_calib_L;
  std::map<double,TF1*> fitFunc_energy_RESIDUAL_TOFHIR_calib_R;
  std::map<double,TF1*> fitFunc_energy_COMBINED_calib_L;
  std::map<double,TF1*> fitFunc_energy_COMBINED_calib_R;
  std::map<double,TF1*> fitFunc_energy_TOFHIR_calib_L;
  std::map<double,TF1*> fitFunc_energy_TOFHIR_calib_R;

  std::map< std::tuple<int,float,int,std::string>, double > MPV_COMBINED_map;
  std::map< std::tuple<int,float,int,std::string>, double > MPV_RESIDUAL_TOFHIR_map;
  std::map< std::tuple<int,float,int,std::string>, double > MPV_TOFHIR_map;
 
  std::map< std::tuple<int,float,int,std::string>, double > Chi2_values_COMBINED_map;
  std::map< std::tuple<int,float,int,std::string>, double > Chi2_values_RESIDUAL_TOFHIR_map;
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

              if (!h1_energy_RESIDUAL_TOFHIR_calib_L[index2]) continue;
              std::string labelLR_energyBin(Form("%s_energyBin%02d",labelLR.c_str(),iEnergyBin));




              // -- energy histogram side L: COMBINED calibration
              c = new TCanvas(Form("c_energy_COMBINED_calib_L_%s",labelLR_energyBin.c_str()),Form("c_energy_COMBINED_calib_L_%s",labelLR_energyBin.c_str()));
              histo = h1_energy_COMBINED_calib_L[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";E_{L} [a.u.];entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              // --- Landau Fit
              auto fitRes_COMBINEDcalib_L = DoLandauFit(histo,minE[std::make_pair(iBar, Vov)],Form("f_landau_COMBINED_bar%02d_%s", iBar, labelLR_energyBin.c_str()));

              // draw
              if (fitRes_COMBINEDcalib_L.func)fitRes_COMBINEDcalib_L.func->Draw("same");

              // fill maps
              if (fitRes_COMBINEDcalib_L.mpv > 0)
              {       
                      MPV_COMBINED_map[std::make_tuple(iBar,Vov,vth1,"L")] = fitRes_COMBINEDcalib_L.mpv;
                      Chi2_values_COMBINED_map[std::make_tuple(iBar,Vov,vth1,"L")] = fitRes_COMBINEDcalib_L.chi2_ndf;
              }

              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d L - COMBINED calibration}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_energy_COMBINED_calib/Left/c_energy_COMBINED_calib_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop3_energy_COMBINED_calib/Left/c_energy_COMBINED_calib_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;


	      // -- energy histogram side R: COMBINED calibration
              c = new TCanvas(Form("c_energy_COMBINED_calib_R_%s",labelLR_energyBin.c_str()),Form("c_energy_COMBINED_calib_R_%s",labelLR_energyBin.c_str()));
              histo = h1_energy_COMBINED_calib_R[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";E_{R} [a.u.];entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              // --- Landau Fit
              auto fitRes_COMBINEDcalib_R = DoLandauFit(histo,minE[std::make_pair(iBar, Vov)],Form("f_landau_COMBINED_bar%02d_%s", iBar, labelLR_energyBin.c_str()));

              // draw
              if (fitRes_COMBINEDcalib_R.func)fitRes_COMBINEDcalib_R.func->Draw("same");

              // fill maps
              if (fitRes_COMBINEDcalib_R.mpv > 0)
              {
                      MPV_COMBINED_map[std::make_tuple(iBar,Vov,vth1,"R")] = fitRes_COMBINEDcalib_R.mpv;
                      Chi2_values_COMBINED_map[std::make_tuple(iBar,Vov,vth1,"R")] = fitRes_COMBINEDcalib_R.chi2_ndf;
              }

              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d R - COMBINED calibration}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_energy_COMBINED_calib/Right/c_energy_COMBINED_calib_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop3_energy_COMBINED_calib/Right/c_energy_COMBINED_calib_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;


              // --- scatter plot energy COMBINED calib L vs R side
              if(!h2_energy_COMBINED_calib_L_vs_R[index2]) continue;
              c = new TCanvas(Form("c_h2_energy_COMBINED_calib_L_vs_R_%s",labelLR_energyBin.c_str()),Form("c_h2_energy_COMBINED_calib_L_vs_R_%s",labelLR_energyBin.c_str()));
              c -> SetGridy();

              h2 = h2_energy_COMBINED_calib_L_vs_R[index2];
              //h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) -600., h2 -> GetMean(2)+600);
              h2 -> SetTitle(Form(";E_{L} [a.u.];E_{R} [a.u.]"));
              h2 -> Draw("colz");

              latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d - COMBINED calib}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_scatter_energy_COMBINED_calib_L_vs_R/c_energy_COMBINED_calib_L_vs_R_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop3_scatter_energy_COMBINED_calib_L_vs_R/c_energy_COMBINED_calib_L_vs_R_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete c;
              delete latex;





	               
	      // -- energy histogram side L: RESIDUAL TOFHIR calibration
              c = new TCanvas(Form("c_energy_RESIDUAL_TOFHIR_calib_L_%s",labelLR_energyBin.c_str()),Form("c_energy_RESIDUAL_TOFHIR_calib_L_%s",labelLR_energyBin.c_str()));
              histo = h1_energy_RESIDUAL_TOFHIR_calib_L[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";E_{L} [a.u.];entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              // --- Landau Fit
              auto fitRes_RESIDUAL_calib_L = DoLandauFit(histo,minE[std::make_pair(iBar, Vov)],Form("f_landau_RESIDUAL__bar%02d_%s", iBar, labelLR_energyBin.c_str()));

              // draw
              if (fitRes_RESIDUAL_calib_L.func)fitRes_RESIDUAL_calib_L.func->Draw("same");

              // fill maps
              if (fitRes_RESIDUAL_calib_L.mpv > 0)
              {
                      MPV_RESIDUAL_TOFHIR_map[std::make_tuple(iBar,Vov,vth1,"L")] = fitRes_RESIDUAL_calib_L.mpv;
                      Chi2_values_RESIDUAL_TOFHIR_map[std::make_tuple(iBar,Vov,vth1,"L")] = fitRes_RESIDUAL_calib_L.chi2_ndf;
              }

              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d L - RESIDUAL TOFHIR calibration}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_energy_RESIDUAL_TOFHIR_calib/Left/c_energy_RESIDUAL_TOFHIR_calib_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop3_energy_RESIDUAL_TOFHIR_calib/Left/c_energy_RESIDUAL_TOFHIR_calib_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));

              delete latex;
              delete c;




              // -- energy histogram side R: RESIDUAL TOFHIR calibration
              c = new TCanvas(Form("c_energy_RESIDUAL_TOFHIR_calib_R_%s",labelLR_energyBin.c_str()),Form("c_energy_RESIDUAL_TOFHIR_calib_R_%s",labelLR_energyBin.c_str()));
              histo = h1_energy_RESIDUAL_TOFHIR_calib_R[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";E_{R} [a.u.];entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();
 
	      // --- Landau Fit
              auto fitRes_RESIDUAL_calib_R = DoLandauFit(histo,minE[std::make_pair(iBar, Vov)],Form("f_landau_RESIDUAL__bar%02d_%s", iBar, labelLR_energyBin.c_str()));

              // draw
              if (fitRes_RESIDUAL_calib_R.func)fitRes_RESIDUAL_calib_R.func->Draw("same");

              // fill maps
              if (fitRes_RESIDUAL_calib_R.mpv > 0)
              {
                      MPV_RESIDUAL_TOFHIR_map[std::make_tuple(iBar,Vov,vth1,"R")] = fitRes_RESIDUAL_calib_R.mpv;
                      Chi2_values_RESIDUAL_TOFHIR_map[std::make_tuple(iBar,Vov,vth1,"R")] = fitRes_RESIDUAL_calib_R.chi2_ndf;
              }

              
              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d R - RESIDUAL TOFHIR calibration}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_energy_RESIDUAL_TOFHIR_calib/Right/c_energy_RESIDUAL_TOFHIR_calib_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop3_energy_RESIDUAL_TOFHIR_calib/Right/c_energy_RESIDUAL_TOFHIR_calib_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete latex;
              delete c;



              // --- scatter plot energy RESIDUAL TOFHIR calib L vs R side
              if(!h2_energy_RESIDUAL_TOFHIR_calib_L_vs_R[index2]) continue;
              c = new TCanvas(Form("c_h2_energy_RESIDUAL_TOFHIR_calib_L_vs_R_%s",labelLR_energyBin.c_str()),Form("c_h2_energy_RESIDUAL_TOFHIR_calib_L_vs_R_%s",labelLR_energyBin.c_str()));
              c -> SetGridy();

              h2 = h2_energy_RESIDUAL_TOFHIR_calib_L_vs_R[index2];
              //h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) -600., h2 -> GetMean(2)+600);
              h2 -> SetTitle(Form(";E_{L} [a.u.];E_{R} [a.u.]"));
              h2 -> Draw("colz");

              latex = new TLatex(0.40,0.85,Form("#splitline{bar %02d - RESIDUAL TOFHIR calib}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
              latex -> SetNDC();
              latex -> SetTextFont(42);
              latex -> SetTextSize(0.04);
              latex -> SetTextColor(kBlack);
              latex -> Draw("same");

              c -> Print(Form("%s/Loop3_scatter_energy_RESIDUAL_TOFHIR_calib_L_vs_R/c_energy_RESIDUAL_TOFHIR_calib_L_vs_R_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
              c -> Print(Form("%s/Loop3_scatter_energy_RESIDUAL_TOFHIR_calib_L_vs_R/c_energy_RESIDUAL_TOFHIR_calib_L_vs_R_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
              delete c;
              delete latex;





             // -- energy histrogram side L: TOFHIR calibration
              c = new TCanvas(Form("c_energy_TOFHIR_calib_L_%s",labelLR_energyBin.c_str()),Form("c_energy_TOFHIR_calib_L_%s",labelLR_energyBin.c_str()));
              histo = h1_energy_TOFHIR_calib_L[index2];
              //histo -> GetXaxis() -> SetRangeUser(histo->GetMean()-5.*histo->GetRMS(),histo->GetMean()+5.*histo->GetRMS());
              histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));
              histo -> SetTitle(Form(";E_{L} [a.u.];entries"));
              histo -> SetLineColor(kRed);
              histo -> SetLineWidth(2);
              histo -> Draw();
              histo -> Write();

              // --- Landau Fit
              auto fitRes_TOFHIR_calib_L = DoLandauFit(histo,minE[std::make_pair(iBar, Vov)],Form("f_landau_TOFHIR_bar%02d_%s", iBar, labelLR_energyBin.c_str()));

              // draw
              if (fitRes_TOFHIR_calib_L.func)fitRes_TOFHIR_calib_L.func->Draw("same");

              // fill maps
              if (fitRes_TOFHIR_calib_L.mpv > 0)
              {
                      MPV_TOFHIR_map[std::make_tuple(iBar,Vov,vth1,"L")] = fitRes_TOFHIR_calib_L.mpv;
                      Chi2_values_TOFHIR_map[std::make_tuple(iBar,Vov,vth1,"L")] = fitRes_TOFHIR_calib_L.chi2_ndf;
              }

              latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d L - TOFHIR calibration}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
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
             histo -> SetTitle(Form(";E_{R} [a.u.];entries"));
             histo -> SetLineColor(kRed);
             histo -> SetLineWidth(2);
             histo -> Draw();
             histo -> Write();
	     
              // --- Landau Fit
              auto fitRes_TOFHIR_calib_R = DoLandauFit(histo,minE[std::make_pair(iBar, Vov)],Form("f_landau_TOFHIR_bar%02d_%s", iBar, labelLR_energyBin.c_str()));

              // draw
              if (fitRes_TOFHIR_calib_R.func)fitRes_TOFHIR_calib_R.func->Draw("same");

              // fill maps
              if (fitRes_TOFHIR_calib_R.mpv > 0)
              {
                      MPV_TOFHIR_map[std::make_tuple(iBar,Vov,vth1,"R")] = fitRes_TOFHIR_calib_R.mpv;
                      Chi2_values_TOFHIR_map[std::make_tuple(iBar,Vov,vth1,"R")] = fitRes_TOFHIR_calib_R.chi2_ndf;
              }
	     
	     latex = new TLatex(0.20,0.85,Form("#splitline{bar %02d R - TOFHIR calibration}{V_{OV} = %.2f V, th. = %d DAC}",iBar,Vov,int(vth1)));
             latex -> SetNDC();
             latex -> SetTextFont(42);
             latex -> SetTextSize(0.04);
             latex -> SetTextColor(kBlack);
             latex -> Draw("same");

             c -> Print(Form("%s/Loop3_energy_TOFHIR_calib/Right/c_energyTOFHIR_calib_%s.png",plotDir.c_str(),labelLR_energyBin.c_str()));
             c -> Print(Form("%s/Loop3_energy_TOFHIR_calib/Right/c_energyTOFHIR_calib_%s.pdf",plotDir.c_str(),labelLR_energyBin.c_str()));
             delete latex;
             delete c;



             // --- scatter plot energy TOFHIR calib L vs R side
             if(!h2_energy_TOFHIR_calib_L_vs_R[index2]) continue;

             c = new TCanvas(Form("c_h2_energy_TOFHIR_calib_L_vs_R_%s",labelLR_energyBin.c_str()),Form("c_h2_energy_TOFHIR_calib_L_vs_R_%s",labelLR_energyBin.c_str()));
             c -> SetGridy();

             h2 = h2_energy_TOFHIR_calib_L_vs_R[index2];
             //h2 -> GetYaxis()->SetRangeUser(h2 -> GetMean(2) -600., h2 -> GetMean(2)+600);
             h2 -> SetTitle(Form(";E_{L} [a.u.];E_{R} [a.u.]"));
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

std::cout << std::endl;
std::cout << "Loop 3 ended" << std::endl;




///////////////////////////////////////////
//      Chi2 and MPV values hisograms    //
///////////////////////////////////////////

// --- Chi2 histograms:
FillAndPlotHistos(Chi2_values_LO_map,"Chi2_LO","#chi^{2}/NDF",50, 0., 10.,0., 50.,plotDir + "/Chi2_plots");
FillAndPlotHistos(Chi2_values_COMBINED_map,"Chi2_COMBINED","#chi^{2}/NDF",50, 0., 10.,0., 50.,plotDir + "/Chi2_plots");
FillAndPlotHistos(Chi2_values_RESIDUAL_TOFHIR_map,"Chi2_RESIDUAL_TOFHIR","#chi^{2}/NDF",50, 0., 10.,0., 50.,plotDir + "/Chi2_plots");
FillAndPlotHistos(Chi2_values_TOFHIR_map,"Chi2_TOFHIR","#chi^{2}/NDF",50, 0., 10.,0., 50.,plotDir + "/Chi2_plots");

// --- MPV histograms:
FillAndPlotHistos(MPV_RAW_map, "MPV_RAW", "MPV [a.u.]", 100, 0., 1000., 50., 900., plotDir+"/MPV_plots");
FillAndPlotHistos(MPV_LO_map, "MPV_LO", "MPV [a.u.]", 100, 0., 1000., 50., 900., plotDir+"/MPV_plots");
FillAndPlotHistos(MPV_RESIDUAL_TOFHIR_map, "MPV_RESIDUALTOFHIR", "MPV [a.u.]", 100, 0., 1000., 50., 900., plotDir+"/MPV_plots");
FillAndPlotHistos(MPV_TOFHIR_map, "MPV_TOFHIR", "MPV [a.u.]", 100, 0., 1000., 50., 900., plotDir+"/MPV_plots");
FillAndPlotHistos(MPV_COMBINED_map, "MPV_COMBINED", "MPV [a.u.]", 100, 0., 1000., 50., 900., plotDir+"/MPV_plots");

std::cout << std::endl;
std::cout << "[MVP AND Chi2 HISTOGRAMS] --> done!" << std::endl;




///////////////////////////////////////
//   Plot energy MPV values vs bar   //
///////////////////////////////////////

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

  StyleGraph(gL, kRed, 20);
  gL->GetYaxis()->SetRangeUser(230.,440.);
  StyleGraph(gR, kBlue, 21);
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
// MPV values histograms energy with LO equalization
//-------------------------------------------------------

for(auto& it : MPV_LO_map)
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

  StyleGraph(gL, kRed, 20);
  gL->GetYaxis()->SetRangeUser(230.,440.);
  StyleGraph(gR, kBlue, 21);

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



//------------------------------
// MPV COMBINED coliabration
//------------------------------ 
for(auto& it : MPV_COMBINED_map)
{
  int bar;
  float vov;
  int thr;
  std::string side;
  double val;

  std::tie(bar, vov, thr, side) = it.first;
  val = it.second;

  std::pair<float,int> key = std::make_pair(vov, thr);

  if( g_MPV_COMBINED_L.find(key) == g_MPV_COMBINED_L.end() )
  {
    g_MPV_COMBINED_L[key] = new TGraph();
    g_MPV_COMBINED_R[key] = new TGraph();
    counter_MPV_COMBINED_L[key] = 0;
    counter_MPV_COMBINED_R[key] = 0;
  }


  if(side == "L")
  {
    g_MPV_COMBINED_L[key]->SetPoint(counter_MPV_COMBINED_L[key], bar, val);
    counter_MPV_COMBINED_L[key]++;
  }
  else if(side == "R")
  {
    g_MPV_COMBINED_R[key]->SetPoint(counter_MPV_COMBINED_R[key], bar, val);
    counter_MPV_COMBINED_R[key]++;
  }
}


for(auto& it : g_MPV_COMBINED_L)
{
  float vov = it.first.first;
  int thr   = it.first.second;

  TGraph* gL = g_MPV_COMBINED_L[it.first];
  TGraph* gR = g_MPV_COMBINED_R[it.first];

  TCanvas* c = new TCanvas(Form("c_MPV_COMBINED_vov%.2f_th%d",vov,thr),Form("COMBINED factors Vov=%.2f th=%d",vov,thr),800,600);

  StyleGraph(gL, kRed, 20);
  gL->GetYaxis()->SetRangeUser(325., 340.);
  StyleGraph(gR, kBlue, 21);
  gL->SetTitle(Form(";Bar index;MPV COMBINED calibration [a.u.]"));

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

  c->Print(Form("%s/MPV_COMBINED_calib/c_MPV_COMBINED_vov%.2f_th%d.pdf",plotDir.c_str(),vov,thr));
  c->Print(Form("%s/MPV_COMBINED_calib/c_MPV_COMBINED_vov%.2f_th%d.png",plotDir.c_str(),vov,thr));

  delete leg;
  delete c;
}


//----------------------------
//   MPV RESIDUAL TOFHIR colibration
//---------------------------- 
for(auto& it : MPV_RESIDUAL_TOFHIR_map)
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
  if( g_MPV_RESIDUAL_TOFHIR_L.find(key) == g_MPV_RESIDUAL_TOFHIR_L.end() )
  {
    g_MPV_RESIDUAL_TOFHIR_L[key] = new TGraph();
    g_MPV_RESIDUAL_TOFHIR_R[key] = new TGraph();
    counter_MPV_RESIDUAL_TOFHIR_L[key] = 0;
    counter_MPV_RESIDUAL_TOFHIR_R[key] = 0;
  }


  if(side == "L")
  {
    g_MPV_RESIDUAL_TOFHIR_L[key]->SetPoint(counter_MPV_RESIDUAL_TOFHIR_L[key], bar, val);
    counter_MPV_RESIDUAL_TOFHIR_L[key]++;
  }
  else if(side == "R")
  {
    g_MPV_RESIDUAL_TOFHIR_R[key]->SetPoint(counter_MPV_RESIDUAL_TOFHIR_R[key], bar, val);
    counter_MPV_RESIDUAL_TOFHIR_R[key]++;
  }
}

for(auto& it : g_MPV_RESIDUAL_TOFHIR_L)
{
  float vov = it.first.first;
  int thr   = it.first.second;

  TGraph* gL = g_MPV_RESIDUAL_TOFHIR_L[it.first];
  TGraph* gR = g_MPV_RESIDUAL_TOFHIR_R[it.first];

  TCanvas* c = new TCanvas(Form("c_MPV_RESIDUAL_TOFHIR_vov%.2f_th%d",vov,thr),Form("RESIDUAL TOFHIR factors Vov=%.2f th=%d",vov,thr),800,600);

  StyleGraph(gL, kRed, 20);
  gL->GetYaxis()->SetRangeUser(285, 385);
  StyleGraph(gR, kBlue, 21);
  gL->SetTitle(Form(";Bar index;MPV RESIDUAL TOFHIR [a.u.]"));

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

  c->Print(Form("%s/MPV_RESIDUAL_TOFHIR_calib/c_MPV_RESIDUAL_TOFHIR_vov%.2f_th%d.pdf",plotDir.c_str(),vov,thr));
  c->Print(Form("%s/MPV_RESIDUAL_TOFHIR_calib/c_MPV_RESIDUAL_TOFHIR_vov%.2f_th%d.png",plotDir.c_str(),vov,thr));

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

  StyleGraph(gL, kRed, 20);
  gL->GetYaxis()->SetRangeUser(285, 385);
  StyleGraph(gR, kBlue, 21);
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

std::cout << std::endl;
std::cout << "MPV plots done!" << std::endl;

// -----------------------------------------------
//  SUMMARY MPV plots (RAW, LO, LO+TOFHIR, TOFHIR)
// -----------------------------------------------

for (auto &it : g_MPV_RAW_L) // loop su tutte le chiavi (vov, vth)
{
    float vov = it.first.first;
    int vth   = it.first.second;

    if (!g_MPV_LO_L.count(it.first)) continue;
    if (!g_MPV_COMBINED_L.count(it.first)) continue;
    if (!g_MPV_RESIDUAL_TOFHIR_L.count(it.first)) continue;

    
    // --- LEFT SIDE 
    TCanvas* c_L = new TCanvas(Form("c_MPV_summary_L_Vov%.2f_th%02d", vov, vth),Form("MPV summary L Vov%.2f th%02d", vov, vth),800, 600);

    auto g_raw_L  = g_MPV_RAW_L[it.first];
    auto g_lo_L   = g_MPV_LO_L[it.first];
    auto g_comb_L = g_MPV_COMBINED_L[it.first];
    auto g_res_L    = g_MPV_RESIDUAL_TOFHIR_L[it.first];

    StyleGraph(g_raw_L, kBlue, 20);
    StyleGraph(g_res_L, kOrange, 21);
    StyleGraph(g_lo_L, kRed, 22);
    StyleGraph(g_comb_L, kGray, 23);

    g_raw_L->SetTitle(Form(";Bar index;MPV"));
    g_raw_L->GetYaxis()->SetRangeUser(150, 550); 

    g_raw_L->Draw("APL");
    g_lo_L->Draw("PL SAME");
    g_comb_L->Draw("PL SAME");
    g_res_L->Draw("PL SAME");

    TLegend* leg_L = new TLegend(0.20,0.70,0.50,0.88);
    leg_L->SetFillStyle(0);
    leg_L->SetBorderSize(0);
    leg_L->SetTextSize(0.03);

    leg_L->SetHeader(Form("Side L  |  V_{OV}=%.2f V  |  vth=%d", vov, vth), "C");

    leg_L->AddEntry(g_raw_L,  "RAW", "lp");
    leg_L->AddEntry(g_res_L,    "RESIDUAL TOFHIR calib", "lp");
    leg_L->AddEntry(g_lo_L,   "LO calib", "lp");
    leg_L->AddEntry(g_comb_L, "COMBINED calib", "lp");

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
    if (!g_MPV_COMBINED_R.count(it.first)) continue;
    if (!g_MPV_RESIDUAL_TOFHIR_R.count(it.first)) continue;

    TCanvas* c_R = new TCanvas(Form("c_MPV_summary_R_Vov%.2f_th%02d", vov, vth),Form("MPV summary R Vov%.2f th%02d", vov, vth),800, 600);

    auto g_raw_R  = g_MPV_RAW_R[it.first];
    auto g_lo_R   = g_MPV_LO_R[it.first];
    auto g_comb_R = g_MPV_COMBINED_R[it.first];
    auto g_res_R    = g_MPV_RESIDUAL_TOFHIR_R[it.first];

    StyleGraph(g_raw_R, kBlue, 20);
    StyleGraph(g_res_R, kOrange, 21);
    StyleGraph(g_lo_R, kRed, 22);
    StyleGraph(g_comb_R, kGray, 23);

    g_raw_R->SetTitle(Form(";Bar index;MPV"));
    g_raw_R->GetYaxis()->SetRangeUser(150, 550);

    g_raw_R->Draw("APL");
    g_lo_R->Draw("PL SAME");
    g_comb_R->Draw("PL SAME");
    g_res_R->Draw("PL SAME");

    TLegend* leg_R = new TLegend(0.20,0.70,0.50,0.88);
    leg_R->SetFillStyle(0);
    leg_R->SetBorderSize(0);
    leg_R->SetTextSize(0.03);

    leg_R->SetHeader(Form("Side R  |  V_{OV}=%.2f V  |  vth=%d", vov, vth), "C");

    leg_R->AddEntry(g_raw_R,  "RAW", "lp");
    leg_R->AddEntry(g_lo_R,   "LO calib", "lp");
    leg_R->AddEntry(g_comb_R, "COMBINED calib", "lp");
    leg_R->AddEntry(g_res_R,  "RESIDUAL TOFHIR calib", "lp");

    leg_R->Draw();

    c_R->SetGridy();
    c_R->SetTickx();
    c_R->SetTicky();

    c_R->Print(Form("%s/MPV_summary/MPV_summary_R_Vov%.2f_th%02d.pdf", plotDir.c_str(), vov, vth));
    c_R->Print(Form("%s/MPV_summary/MPV_summary_R_Vov%.2f_th%02d.png", plotDir.c_str(), vov, vth));

    delete c_R;
}
std::cout << std::endl;
std::cout << "Summary MPV plots done!" << std::endl;


  int bytes = outFile -> Write();
  std::cout << "============================================"  << std::endl;
  std::cout << "nr of  B written:  " << int(bytes)             << std::endl;
  std::cout << "nr of KB written:  " << int(bytes/1024.)       << std::endl;
  std::cout << "nr of MB written:  " << int(bytes/1024./1024.) << std::endl;
  std::cout << "============================================"  << std::endl;
}
