#include "interface/AnalysisUtils.h"
#include "interface/FitUtils.h"


ClassImp(ModuleEventWithRefClass)

void TrackProcess(float* cpu, float* mem, float* vsz, float* rss)
{
  std::string dummy1, dummy2, dummy3, dummy4, dummy5;
  std::string time;  
  //---get cpu/mem info      
  int pid = getpid();
  std::string ps_command = "ps up "+std::string(Form("%d",pid))+" >.proc.tmp";
  system(ps_command.c_str());
  std::ifstream proc_tmp(".proc.tmp", std::ios::in);
  getline(proc_tmp, dummy1);
  proc_tmp >> dummy1 >> dummy2 >> cpu[0] >> mem[0] >> vsz[0] >> rss[0]
           >> dummy3 >> dummy4 >> dummy5 >> time;
  vsz[0] = vsz[0]/1000;
  rss[0] = rss[0]/1000;
  proc_tmp.close();
  if(cpu[0]>cpu[1])
    cpu[1] = cpu[0];
  if(mem[0]>mem[1])
    mem[1] = mem[0];
  if(vsz[0]>vsz[1])
    vsz[1] = vsz[0];
  if(rss[0]>rss[1])
    rss[1] = rss[0];  
  //---print statistics
  std::cout << "\033[0m""-----Machine stats---current/max-----" << std::endl
            << "CPU(%): " << cpu[0] << "/" << cpu[1] << std::endl
            << "MEM(%): " << mem[0] << "/" << mem[1] << std::endl
            << "VSZ(M): " << vsz[0] << "/" << vsz[1] << std::endl
            << "RSS(M): " << rss[0] << "/" << rss[1] << std::endl
            << "time lasted: " << time << std::endl;
}

std::vector<std::string> GetTokens(const std::string& input, const char& sep)
{
  std::stringstream ss(input);
  std::string token;
  std::vector<std::string> tokens;  
  while( std::getline(ss,token,sep) ) {
    tokens.push_back(token);
  }  
  return tokens;
}

float DeltaEta(const float& eta1, const float& eta2)
{
  return fabs( eta1 - eta2 );
}

float DeltaPhi(const float& phi1, const float& phi2)
{
  float dphi = fabs( phi1 - phi2 );
  if( dphi > PI ) dphi = 2*PI - dphi;
  return dphi;
}

float DeltaR(const float& eta1, const float& phi1, const float& eta2, const float& phi2)
{
  return sqrt( DeltaEta(eta1,eta2)*DeltaEta(eta1,eta2) +
               DeltaPhi(phi1,phi2)*DeltaPhi(phi1,phi2) );
}

float FindXMaximum(TH1F* histo, const float& xMin, const float& xMax, const bool& checkDerivative)
{
  float max = -999999999.;
  int binMax = -1;
  for(int bin = 2; bin <= histo->GetNbinsX()-1; ++bin)
  {
    if( histo->GetBinCenter(bin) < xMin ) continue;
    if( histo->GetBinCenter(bin) > xMax ) continue;
    if( histo->GetBinContent(bin) > max )
    {
      float delta_min = ( histo->GetBinContent(bin-1) - histo->GetBinContent(bin) ) / histo->GetBinContent(bin);
      float delta_max = ( histo->GetBinContent(bin+1) - histo->GetBinContent(bin) ) / histo->GetBinContent(bin);
      if( (fabs(delta_min) < 0.2 && fabs(delta_max) < 0.2) || !checkDerivative )
	{
	  max = histo->GetBinContent(bin); binMax = bin;
	}
    };
  }
  return histo->GetBinCenter(binMax);
}

int FindBin(const float& val, const std::vector<float>* ranges)
{
  for(unsigned int ii = 0; ii < ranges->size()-1; ++ii)
  {
    if( val >= ranges->at(ii) && val < ranges->at(ii+1) ) return ii;
  }  
  return -1;
}

// -- read TOFHIR calibration csv and extract the calibMap to be used in step1 ---
std::map<std::tuple<int,std::string,int,int>, float> loadCalibrationMap(const std::string& filename)
{
  std::map<std::tuple<int,std::string,int,int>, float> calibMap;
  std::ifstream fin(filename);
  if (!fin.is_open()) {
    std::cerr << "[ERROR] Error in opening file: " << filename << std::endl;
    return calibMap;
  }
  std::string line;
  // skip header
  getline(fin, line);
  while (getline(fin, line)){
    std::stringstream ss(line);
    std::string token;
    int bar, thr;
    float vov;
    float calib;
    std::string side;
    getline(ss, token, ','); bar = atoi(token.c_str());
    getline(ss, token, ','); side = token;
    getline(ss, token, ','); vov = atof(token.c_str());
    getline(ss, token, ','); thr = atoi(token.c_str());
    getline(ss, token, ','); calib = atof(token.c_str());
    int vov_int = int(vov * 100 + 0.5);
    calibMap[std::make_tuple(bar, side, vov_int, thr)] = calib;
  }
  std::cout << ">>> Loaded calibration file: " << filename << std::endl;
  std::cout << "    Calibration entries loaded: " << calibMap.size() << std::endl;
  return calibMap;
}

// this is currently not used in the analysis (minE is set manually in the config, allowing to use a simple landau fit of the MIP peak)
Double_t langaufun(Double_t *x, Double_t *par)
{  
  // --------------------
  //  Numerical convolution of a Landau energy-loss distribution with a Gaussian resolution function.
  //  Used to fit energy spectra including detector smearing effects.
  //  Fit parameters:
  //  - par[0] = Width (scale) parameter of Landau density
  //  - par[1] = Most Probable (MP, location) parameter of Landau density
  //  - par[2] = Total area (integral -inf to inf, normalization constant)
  //  - par[3] = Width (sigma) of convoluted Gaussian function
  //  In the Landau distribution (represented by the CERNLIB approximation),
  //   the maximum is located at x=-0.22278298 with the location parameter=0.
  //  This shift is corrected within this function, so that the actual
  //   maximum is identical to the MP parameter.
  // --------------------
  
  // constants
  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location  
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

void GetEnergyBins(TH1F *h, std::vector<float> *r, std::map<int, float> & b){
  // --------------------
  // Divide the input energy histogram into sub-ranges within the values in vector r (energy ranges)
  //  and computes the mean energy in each sub-range, storing the result in map b.
  // --------------------
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

void drawDeltaT(TCanvas *& c, TH1F *histo, TF1 *& fitFunc, std::string xaxis_label, std::string latex_label, std::string drawSame ){
  // --------------------  
  // Plots the time difference histogram, refines a Gaussian fit near the maximum,
  //  and write effective sigma and fitted Gaussian resolution.
  // --------------------  
  c->cd();
  histo -> SetTitle(Form("; %s #Deltat [ps];entries", xaxis_label.c_str()));
  histo -> Draw(drawSame.c_str());
  float* vals = new float[6];
  FindSmallestInterval(vals,histo,0.68);
  float min = vals[4];
  float max = vals[5];
  float delta = max-min;
  float sigma = 0.5*delta;
  float effSigma = sigma;
  float fitXMin = histo->GetBinCenter(histo->GetMaximumBin()) - 200.;
  float fitXMax = histo->GetBinCenter(histo->GetMaximumBin()) + 200.;
  fitFunc -> SetParameters(1,histo->GetMean(),histo->GetRMS());
  fitFunc -> SetRange(fitXMin, fitXMax);
  histo -> Fit(fitFunc,"QNRSL");
  fitFunc -> SetRange(fitFunc->GetParameter(1)-1.0*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+1.0*fitFunc->GetParameter(2));
  histo -> Fit(fitFunc,"QNRSL");
  fitFunc -> SetRange(fitFunc->GetParameter(1)-2.5*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+2.5*fitFunc->GetParameter(2));
  histo -> Fit(fitFunc,"QRSL+");
  fitFunc -> SetLineColor( histo -> GetLineColor() + 1 );
  fitFunc -> SetLineWidth(3);
  fitFunc -> Draw("same");
  histo -> SetMaximum(histo->GetMaximum()+0.1*histo->GetMaximum());
  histo -> GetXaxis() -> SetRangeUser(fitFunc->GetParameter(1)-7.*fitFunc->GetParameter(2),fitFunc->GetParameter(1)+7.*fitFunc->GetParameter(2));
  TLatex* latex = new TLatex(0.55,0.85,Form("#splitline{#sigma_{%s}^{eff} = %.0f ps}{#sigma_{%s}^{gaus} = %.0f ps}",latex_label.c_str(),effSigma, latex_label.c_str(),fabs(fitFunc->GetParameter(2))));
  if (drawSame == "same")
    latex = new TLatex(0.20,0.85,Form("#splitline{#sigma_{%s}^{eff} = %.0f ps}{#sigma_{%s}^{gaus} = %.0f ps}",latex_label.c_str(),effSigma, latex_label.c_str(),fabs(fitFunc->GetParameter(2))));
  latex -> SetNDC();
  latex -> SetTextFont(42);
  latex -> SetTextSize(0.04);
  latex -> SetTextColor( histo -> GetLineColor() );
  latex -> Draw("same");
  delete vals;
  delete latex;
}

int findBarAndSide(int ch, int asic, std::vector<unsigned int>& channelMapping, std::string& side)
{
  // ----------------------------------
  // Given a channel, ASIC and channelMapping, this function returns the corresponding bar and side
  // ----------------------------------
  int localChannel = ch - 32*asic;
  for (size_t iBar = 0; iBar < channelMapping.size()/2; iBar++) {
    if (channelMapping[iBar*2] == localChannel) {
      side = "L";
      return iBar;
    }
    if (channelMapping[iBar*2+1] == localChannel) {
      side = "R";
      return iBar;
    }
  }
  side = "UNKNOWN";
  return -1;
}

void applyCalibration(float& energy, int bar, const std::string& side, float Vov, float vth, const CalibMap& calibMap)
{
  int vov_int = int(Vov * 100 + 0.5);
  int thr_int = int(vth + 0.5);
  CalibKey key = std::make_tuple(bar, side, vov_int, thr_int);
  auto it = calibMap.find(key);
  if (it != calibMap.end()) {
    energy *= it->second;
  }
}
