#ifndef ENE_INTERCALIB_UTILS_H
#define ENE_INTERCALIB_UTILS_H

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/FitUtils.h"
#include "interface/AnalysisUtils.h"
#include "interface/AmplitudeWalk_Utils.h"

#include <utility>
#include <cmath>
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
#include <set>

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



/*
//------------------------------------------------------------------------
//---- find energy bins
//------------------------------------------------------------------------
inline void GetEnergyBins(TH1F *h, std::vector<float> *r, std::map<int, float> & b){

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
inline std::pair<float,float> GetDynamicFitRange(TProfile* prof, float marginFraction = 0.05)
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

*/

//------------------------------------------------------------------------
// Compute the RMS and the mean from a TGraph plot
//------------------------------------------------------------------------
inline void computeMeanRMS(TGraph* g, double& mean, double& rms)
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

inline FitResult DoLandauFit(TH1* histo,double minE,const std::string& name){
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
inline void SetAutoYRange(TGraph* g, double frac = 0.05)
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

/*
//------------------------------------------------------------------------
// --- helper function TGraph style:
//------------------------------------------------------------------------
inline void StyleGraph(TGraph* g, int color, int marker)
{
    g->SetMarkerStyle(marker);
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetLineWidth(2);
}

*/





//------------------------------------------------------------------------
//--- helper function for Chi2 and MPV histos:
//------------------------------------------------------------------------
inline void FillAndPlotHistos(
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
inline void BuildGraphs(
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
inline void PlotGraphs(
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

        TCanvas* c = new TCanvas(Form("c_%s_vov%.2f_th%d", label.c_str(), vov, thr),"", 800, 600);

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

        c->Print(Form("%s/%s_vov%.2f_th%d.pdf",outDir.c_str(), label.c_str(), vov, thr));
        c->Print(Form("%s/%s_vov%.2f_th%d.png",outDir.c_str(), label.c_str(), vov, thr));

        delete leg;
        delete c;
    }
}




//------------------------------------------------------------------------
// --- helper function TGaph with MPV Ratio (Left/right):
//------------------------------------------------------------------------
inline void BuildMPVRatioGraphs(
    const std::map<std::pair<float,int>, TGraph*>& gL_map,
    const std::map<std::pair<float,int>, TGraph*>& gR_map,
    std::map<std::pair<float,int>, TGraph*>& gRatio_map
)
{
    for(const auto& it : gL_map)
    {
        auto key = it.first;

        TGraph* gL = gL_map.at(key);
        TGraph* gR = gR_map.at(key);

        TGraph* gRatio = new TGraph();

        int n = 0;

        for(int i = 0; i < gL->GetN(); ++i)
        {
            double xL, yL;
            double xR, yR;

            gL->GetPoint(i, xL, yL);
            gR->GetPoint(i, xR, yR);

            if (xL != xR) continue; // sicurezza

            if (yR == 0) continue;

            double ratio = yL / yR;

            gRatio->SetPoint(n, xL, ratio);
            n++;
        }

        gRatio_map[key] = gRatio;
    }
}




//------------------------------------------------------------------------
// --- helper function plot of the TGaph with MPV Ratio (Left/right):
//------------------------------------------------------------------------
inline void PlotMPVRatioGraphs(
    const std::map<std::pair<float,int>, TGraph*>& gRatio_map,
    const std::string& label,
    const std::string& outDir
)
{
    for(const auto& it : gRatio_map)
    {
        float vov = it.first.first;
        int thr   = it.first.second;

        TGraph* g = it.second;

        TCanvas* c = new TCanvas(Form("c_MPV_ratio_%s_vov%.2f_th%d", label.c_str(), vov, thr),"", 800, 600);

        StyleGraph(g, kBlue, 20);

        SetAutoYRange(g);

        g->SetTitle(";Bar index; MPV_{L} / MPV_{R}");

        g->Draw("APL");

        TLine* line = new TLine(0,1,16,1);
        line->SetLineStyle(2);
        line->Draw("same");

        double mean, rms;
        computeMeanRMS(g, mean, rms);

        TLatex latex;
        latex.SetNDC();
        latex.SetTextSize(0.04);
        latex.DrawLatex(0.20,0.85,Form("Mean = %.3f", mean));
        latex.DrawLatex(0.20,0.80,Form("RMS  = %.3f", rms));

        latex.DrawLatex(0.20,0.75,Form("V_{OV} = %.2f V, thr = %d", vov, thr));

        c->SetGridy();

        c->Print(Form("%s/MPV_ratio_%s_vov%.2f_th%d.pdf",outDir.c_str(), label.c_str(), vov, thr));
        c->Print(Form("%s/MPV_ratio_%s_vov%.2f_th%d.png",outDir.c_str(), label.c_str(), vov, thr));

        delete c;
    }
}






void SaveCalibrationFactorsToFile(
    const std::map<std::tuple<int,float,int,std::string>, double>& factors,
    const std::string& outputDir,
    const std::string& calibLabel,
    const std::string& LO_mode
)
{
    // ----------------------------------------
    // output file name
    // ----------------------------------------

    std::string outFileName =
        Form("%s/%s_calibration_factors_%s.csv",
             outputDir.c_str(),
             calibLabel.c_str(),
             LO_mode.c_str());

    std::ofstream outFile(outFileName);

    if( !outFile.is_open() )
    {
        std::cerr << "ERROR: cannot open output file "
                  << outFileName << std::endl;
        return;
    }

    // ----------------------------------------
    // CSV header
    // ----------------------------------------

    outFile << "bar,side,vov,th,calib\n";

    // ----------------------------------------
    // ordered map
    // ----------------------------------------

    std::map<
        std::tuple<float,int,int,std::string>,
        double
    > orderedMap;

    for(const auto& it : factors)
    {
        int bar;
        float vov;
        int thr;
        std::string side;

        std::tie(bar,vov,thr,side) = it.first;

        orderedMap[
            std::make_tuple(vov,thr,bar,side)
        ] = it.second;
    }

    // ----------------------------------------
    // write ordered entries
    // ----------------------------------------

    for(const auto& it : orderedMap)
    {
        float vov;
        int thr;
        int bar;
        std::string side;

        std::tie(vov,thr,bar,side) = it.first;

        double calib = it.second;

        outFile
            << bar   << ","
            << side  << ","
            << vov   << ","
            << thr   << ","
            << calib << "\n";
    }

    outFile.close();

    std::cout << "[WRITE] Saved file: "
              << outFileName << std::endl;
}




#endif
