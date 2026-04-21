#ifndef AMPLITUDEWALK_UTILS_H
#define AMPLITUDEWALK_UTILS_H

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/FitUtils.h"
#include "interface/AnalysisUtils.h"

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


// HELPER FUNCTIONS:


//---- find energy bins
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




inline float GetProfileBinCorrection(TProfile* prof, float x)
{
    if (!prof) return 0.0;

    int bin = prof->FindBin(x);

    if (bin < 1 || bin > prof->GetNbinsX())
        return 0.0;

    if (prof->GetBinEntries(bin) <= 0)
        return 0.0;

    return prof->GetBinContent(bin);
}





//-----------------------------------------------------
// --- Helper function for computing the CTR From Histo
// ----------------------------------------------------
inline void ComputeCTRFromHistoMap(
    std::map<double, TH1F*>& hMap,
    std::map<double, float>& CTRMeans,
    std::map<double, float>& CTRSigmas
)
{
    float vals[10];

    for (auto& mapIt : hMap)
    {
        double index = mapIt.first;
        TH1F* histo = mapIt.second;

        if (!histo) continue;

        FindSmallestInterval(vals, histo, 0.68);

        float mean  = vals[0];
        float min   = vals[4];
        float max   = vals[5];
        float sigma = 0.5 * (max - min);

        CTRMeans[index]  = mean;
        CTRSigmas[index] = sigma;
    }
}



//-----------------------------------------------------
// --- Helper function latex for histrograms
// ----------------------------------------------------
inline std::string BuildLatexSplit(
    int iBar,
    float Vov,
    int vth1,
    TH1F* histo,
    const std::string& side // "L", "R", oppure ""
    //float xmin,
    //float xmax
)
{
    if (!histo) return "";

    return Form(
        "#splitline{bar %02d %s}"
        "{#splitline{V_{OV} = %.2f V, th. = %d DAC}"
        "{Mean = %.3f, stdDev = %.3f}}",
        iBar,
        side.c_str(),
        Vov,
        vth1,
        histo->GetMean(),
        histo->GetRMS()
    );
}




//------------------------------------------------------------------------
// Compute the RMS and the mean from a TGraph plot
//------------------------------------------------------------------------
inline void computeMeanRMS(TGraphErrors* g, double& mean, double& rms)
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






// FIT FUNCTIONS:

//---- Find fit range profile
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





//-----------------------------------------------------
// --- Helper function plotting histograms (no fit)
// ----------------------------------------------------
inline void DrawHistogram(
    TH1F* histo,
    const std::string& canvasName,
    const std::string& xTitle,
    const std::string& plotDir,
    const std::string& subDir,
    const std::string& label,
    const std::string& latexText,
    float xMin = -999,
    float xMax = -999,
    bool useLogY = false
)
{
    if (!histo) return;

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str());

    float mean = histo->GetMean();
    float rms  = histo->GetRMS();


    // --- Range X
    if (xMin < xMax)
        histo->GetXaxis()->SetRangeUser(xMin, xMax);
    else
        histo->GetXaxis()->SetRangeUser(mean - 5.*rms, mean + 5.*rms);

    // --- Range Y
    histo->SetMaximum(1.25 * histo->GetBinContent(histo->FindBin(FindXMaximum(histo, mean - 2.*rms, mean + 2.*rms))));

    // --- Stile
    histo->SetTitle(Form(";%s;entries", xTitle.c_str()));
    histo->SetLineColor(kRed);
    histo->SetLineWidth(2);

    if (useLogY) c->SetLogy();

    histo->Draw();
    histo->Write();

    // --- Latex (LEFT aligned)
    TLatex* latex = new TLatex(0.20, 0.90, latexText.c_str());

    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.04);
    latex->SetTextColor(kBlack);
    latex->SetTextAlign(13);
    latex->Draw("same");

    // --- Save
    c->Print(Form("%s/%s/%s.png", plotDir.c_str(), subDir.c_str(), label.c_str()));
    c->Print(Form("%s/%s/%s.pdf", plotDir.c_str(), subDir.c_str(), label.c_str()));


    delete latex;
    delete c;
}






//-----------------------------------------------------
// --- Helper function plotting histograms (with fit)
// ----------------------------------------------------
inline TF1* DrawAndFitHistogram(TH1F* histo,
                         const std::string& canvasName,
                         const std::string& xTitle,
                         const std::string& plotDir,
                         const std::string& subDir,
                         int iBar,
                         float Vov,
                         int vth1,
                         const std::string& LRlabel = "",
			 // -- optional
			 int Fit_LineColour = kBlack,
			 float Xaxis_Range = 5.)
{
    if (!histo) return nullptr;

    // Canvas
    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str());

    histo -> GetXaxis() -> SetRangeUser(histo->GetMean() - Xaxis_Range*histo->GetRMS(),histo->GetMean() + Xaxis_Range*histo->GetRMS());
    histo -> SetMaximum(1.25*histo->GetBinContent(histo->FindBin(FindXMaximum(histo,histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS()))));

    histo->SetTitle(Form(";%s;entries", xTitle.c_str()));
    histo->SetLineColor(kRed);
    histo->SetLineWidth(2);
    histo->Draw();
    histo->Write();


    // --- Fit Gaussiano iterativo
    TF1* fitFunc = new TF1(Form("fit_%s", canvasName.c_str()),"gaus",histo->GetMean()-2.*histo->GetRMS(),histo->GetMean()+2.*histo->GetRMS());


    // 1° fit
    histo->Fit(fitFunc, "QNRS");

    // 2° fit (range centrato)
    histo->Fit(fitFunc, "QSR+","",fitFunc->GetParameter(1) - 2.*fitFunc->GetParameter(2),fitFunc->GetParameter(1) + 2.*fitFunc->GetParameter(2));

    // 3° fit (refinement)
    histo->Fit(fitFunc, "QSR+","",fitFunc->GetParameter(1) - 2.*fitFunc->GetParameter(2),fitFunc->GetParameter(1) + 2.*fitFunc->GetParameter(2));

    //fitFunc->SetLineColor(kBlack);
    fitFunc->SetLineColor(Fit_LineColour);
    fitFunc->SetLineWidth(2);
    fitFunc->Draw("same");

    TLatex* latex = new TLatex(0.20, 0.85,
        Form("#splitline{bar %02d %s}{#splitline{V_{OV} = %.2f V, th. = %d DAC}{Mean = %.3f, stdDev = %.3f}}",
             iBar,
             LRlabel.c_str(),
             Vov,
             vth1,
             histo->GetMean(),
             histo->GetRMS()
        )
    );


    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.04);
    latex->SetTextColor(kBlack);
    latex->Draw("same");

    c->Print(Form("%s/%s/%s.png", plotDir.c_str(), subDir.c_str(), canvasName.c_str()));
    c->Print(Form("%s/%s/%s.pdf", plotDir.c_str(), subDir.c_str(), canvasName.c_str()));

    delete latex;
    delete c;

    return fitFunc;
}








//-----------------------------------------------------
// --- Helper function plotting TProfile (with fit)
// ----------------------------------------------------
inline TF1* DrawAndFitProfile(TProfile* profile,
                         std::map<double,float>& CTRMeans_map,
                         std::map<double,float>& CTRSigmas_map,
                         std::map<std::string, std::map<int, std::vector<float>*> >& ranges_map,
                         bool DynamicFitRange,
                         double map_index1,
                         double map_index2,
                         const std::string& polFunc,
                         const std::string& canvasName,
                         const std::string& xTitle,
                         const std::string& yTitle,
                         const std::string& plotDir,
                         const std::string& subDir,
                         int iBar,
                         float Vov,
                         int vth1,
                         const std::string& LRlabel = "",

                         // --- optional settings
                         bool useCustomXRange = false,
                         float xMin = 0., float xMax = 0.,
                         bool useCustomYRange = false,
                         float yMin = 0., float yMax = 0.
                         )
{
        if (!profile) return nullptr;


        // Canvas
        TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str());

        float mean = CTRMeans_map.at(map_index2);
        float sigma = CTRSigmas_map.at(map_index2);

        profile ->SetTitle(Form(";%s;%s", xTitle.c_str(), yTitle.c_str()));
        profile -> GetYaxis() -> SetRangeUser(mean - 5.*sigma, mean + 5.*sigma);
        if (useCustomXRange) profile -> GetXaxis() -> SetRangeUser(xMin,xMax);
        if (useCustomYRange) profile -> GetYaxis() -> SetRangeUser(yMin,yMax);
        profile -> Draw("");

        TLatex* latex = new TLatex(0.20, 0.85,
        Form("#splitline{bar %02d %s}{V_{OV} = %.2f V, th. = %d DAC}",iBar,LRlabel.c_str(),Vov,vth1));
        latex->SetNDC();
        latex->SetTextFont(42);
        latex->SetTextSize(0.04);
        latex->SetTextColor(kBlack);
        latex->Draw("same");

        float fitXMin = ranges_map[LRlabel][map_index1]->at(0);
        float fitXMax = ranges_map[LRlabel][map_index1]->at(1);

        if (DynamicFitRange)
        {
                auto rangeFit = GetDynamicFitRange(profile, 0.05);
                fitXMin = rangeFit.first;
                fitXMax = rangeFit.second;
        }


        TF1* fitFunc = new TF1(Form("fit_%s", canvasName.c_str()),polFunc.c_str(),fitXMin,fitXMax);

        profile -> Fit(fitFunc,"QRS+");
        fitFunc -> SetLineColor(kRed);
        fitFunc -> SetLineWidth(2);
        fitFunc -> Draw("same");

        c->Print(Form("%s/%s/%s.png", plotDir.c_str(), subDir.c_str(), canvasName.c_str()));
        c->Print(Form("%s/%s/%s.pdf", plotDir.c_str(), subDir.c_str(), canvasName.c_str()));

        delete latex;
        delete c;

        return fitFunc;
}








//---------------------------------------------------------
// --- Helper function plotting TH2F and TProfile togather
// --------------------------------------------------------

inline TF1* DrawTH2WithProfile(
    TH2F* h2,
    TProfile* prof,
    const std::string& canvasName,
    const std::string& xTitle,
    const std::string& yTitle,
    const std::string& plotDir,
    const std::string& subDir,
    int iBar,
    float Vov,
    int vth1,
    const std::string& LRlabel = "",

    // --- optional settings
    bool useCustomXRange = false,
    float xMin = 0., float xMax = 0.,
    bool useCustomYRange = false,
    float yMin = 0., float yMax = 0.,
    bool useGridY = true
)
{
    if (!h2 && !prof) return nullptr;

    // Canvas
    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str());
    if (useGridY) c->SetGridy();


    // --- Titoli
    h2->SetTitle(Form(";%s;%s", xTitle.c_str(), yTitle.c_str()));

    // --- Range asse X (opzionale)
    if (useCustomXRange)
        h2->GetXaxis()->SetRangeUser(xMin, xMax);

    // --- Range asse Y
    if (useCustomYRange)
    {
        h2->GetYaxis()->SetRangeUser(yMin, yMax);
    }
    else
    {
        // default intelligente: centrato sulla media
        float meanY = h2->GetMean(2);
        h2->GetYaxis()->SetRangeUser(meanY - 600., meanY + 600.);
    }

    // --- Draw TH2
    h2->Draw("colz");

    // --- Overlay profile
    if(prof){
    prof->SetLineColor(kRed);
    prof->SetLineWidth(3);
    prof->Draw("same");
    }


    // --- Latex
    TLatex* latex = new TLatex(0.40, 0.85,
        Form("#splitline{bar %02d %s}{V_{OV} = %.2f V, th. = %d DAC}",
             iBar,
             LRlabel.c_str(),
             Vov,
             vth1)
    );

    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.04);
    latex->SetTextColor(kRed);
    latex->Draw("same");

    // --- Save
    c->Print(Form("%s/%s/%s.png", plotDir.c_str(), subDir.c_str(), canvasName.c_str()));
    c->Print(Form("%s/%s/%s.pdf", plotDir.c_str(), subDir.c_str(), canvasName.c_str()));

    delete latex;
    delete c;

    return nullptr;
}


//-------------------------------------------------------------------------------------------
// --- helper function TGaph Left and right same:
//------------------------------------------------------------------------
inline void BuildGraphs(
    const std::map<std::tuple<int,float,int,std::string>, std::pair<double,double>>& input_map,
    std::map<std::pair<float,int>, TGraphErrors*>& gL_map,
    std::map<std::pair<float,int>, TGraphErrors*>& gR_map
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
        double eval;

        std::tie(bar, vov, thr, side) = it.first;
        val = it.second.first;
        eval = it.second.second;

        std::pair<float,int> key = std::make_pair(vov, thr);

        if(gL_map.find(key) == gL_map.end())
        {
            gL_map[key] = new TGraphErrors();
            gR_map[key] = new TGraphErrors();
            counterL[key] = 0;
            counterR[key] = 0;
        }

            if(side == "L")
        {
            int i = counterL[key]++;
            gL_map[key]->SetPoint(i, bar+1, val);
            gL_map[key]->SetPointError(i, 0., eval);
        }
        else if(side == "R")
        {
            int i = counterR[key]++;
            gR_map[key]->SetPoint(i, bar+1, val);
            gR_map[key]->SetPointError(i, 0., eval);
        }

    }
}


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




//------------------------------------------------------------------------
// --- helper function TGaph Left and right same:
//------------------------------------------------------------------------
inline void PlotGraphs(
    const std::map<std::pair<float,int>, TGraphErrors*>& gL_map,
    const std::map<std::pair<float,int>, TGraphErrors*>& gR_map,
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

        TGraphErrors* gL = gL_map.at(it.first);
        TGraphErrors* gR = gR_map.at(it.first);

        TCanvas* c = new TCanvas(Form("c_%s_vov%.2f_th%d", label.c_str(), vov, thr),"", 800, 600);

        StyleGraph(gL, kRed, 20);
        StyleGraph(gR, kBlue, 21);


        gL->GetYaxis()->SetRangeUser(ymin, ymax);
        gL->GetXaxis()->SetRangeUser(0., 17.);

        gL->SetTitle(Form(";Bar index; %s", yTitle.c_str()));


        gL->Draw("APL");
        gR->Draw("PL SAME");

        // stats
        double mean_L, rms_L;
        double mean_R, rms_R;

        computeMeanRMS(gL, mean_L, rms_L);
        computeMeanRMS(gR, mean_R, rms_R);

        TLegend* leg = new TLegend(0.20,0.70,0.60,0.9);
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
// --- helper function to plot togather different amplitude walk function:
//------------------------------------------------------------------------

inline void DrawOverlayFitsFromIndices(
    const std::map<double, TF1*>& fitMap,
    const std::vector<double>& selectedIndices,
    const std::string& canvasName,
    const std::string& outDir,
    //const std::string& subdir,
    const std::string& header = ""
)
{
    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);

    TLegend* leg = new TLegend(0.50,0.65,0.8,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.03);

    int color = 1;
    bool first = true;

    for(double idx : selectedIndices)
    {
        if(fitMap.find(idx) == fitMap.end()) continue;

        TF1* func = fitMap.at(idx);
        if(!func) continue;


        // stile
        func->SetLineColor(color);
        func->SetLineWidth(2);

        float bar = idx - 13001200;

        // draw
        if(first)
        {
            func->SetTitle(";Energy [a.u.];#Delta T raw [ps]");
            func->Draw();
            first = false;
        }
        else
        {
            func->Draw("same");
        }

        leg->AddEntry(func, Form("bar = %.0f", bar), "l");

        color++;
        if(color == 10) color = 1;
    }

    if(header != "")
        leg->SetHeader(header.c_str(), "C");

    leg->Draw();

    c->SetGrid();


    c->Print(Form("%s/%s.pdf", outDir.c_str(), canvasName.c_str()));
    c->Print(Form("%s/%s.png", outDir.c_str(), canvasName.c_str()));

    delete leg;
    delete c;
}






inline TF1* DrawAndFitProfile_Iterative(TProfile* profile,
    std::map<double,float>& CTRMeans_map,
    std::map<double,float>& CTRSigmas_map,
    std::map<std::string, std::map<int, std::vector<float>*> >& ranges_map,
    bool DynamicFitRange,
    double map_index1,
    double map_index2,
    const std::string& polFunc,
    const std::string& canvasName,
    const std::string& xTitle,
    const std::string& yTitle,
    const std::string& plotDir,
    const std::string& subDir,
    int iBar,
    float Vov,
    int vth1,
    const std::string& LRlabel = "",

    // optional
    bool useCustomXRange = false,
    float xMin = 0., float xMax = 0.,
    bool useCustomYRange = false,
    float yMin = 0., float yMax = 0.,

    // new parameter
    float nSigmaCut = 2.5
)
{
    if (!profile) return nullptr;

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str());

    float mean  = CTRMeans_map.at(map_index2);
    float sigma = CTRSigmas_map.at(map_index2);

    profile->SetTitle(Form(";%s;%s", xTitle.c_str(), yTitle.c_str()));
    profile->GetYaxis()->SetRangeUser(mean - 5.*sigma, mean + 5.*sigma);

    if (useCustomXRange) profile->GetXaxis()->SetRangeUser(xMin,xMax);
    if (useCustomYRange) profile->GetYaxis()->SetRangeUser(yMin,yMax);

    profile->Draw();

    // --- Latex
    TLatex* latex = new TLatex(0.20, 0.85,
        Form("#splitline{bar %02d %s}{V_{OV} = %.2f V, th. = %d DAC}",
        iBar, LRlabel.c_str(), Vov, vth1));

    latex->SetNDC();
    latex->SetTextFont(42);
    latex->SetTextSize(0.04);
    latex->Draw("same");

    // --- Fit range
    float fitXMin = ranges_map[LRlabel][map_index1]->at(0);
    float fitXMax = ranges_map[LRlabel][map_index1]->at(1);

    if (DynamicFitRange)
    {
        auto rangeFit = GetDynamicFitRange(profile, 0.05);
        fitXMin = rangeFit.first;
        fitXMax = rangeFit.second;
    }

    // =========================================================
    // --- 1° FIT (grezzo)
    // =========================================================
    TF1* fitFunc = new TF1(Form("fit1_%s", canvasName.c_str()),
                          polFunc.c_str(), fitXMin, fitXMax);

    profile->Fit(fitFunc,"QRS");

    // =========================================================
    // --- COSTRUZIONE TGraph filtrato
    // =========================================================
    TGraphErrors* gClean = new TGraphErrors();

    int point = 0;
    double nMin = 30;


    for(int i = 1; i <= profile->GetNbinsX(); ++i)

    {
      double x   = profile->GetBinCenter(i);
      double y   = profile->GetBinContent(i);
      double err = profile->GetBinError(i);
      double n   = profile->GetBinEntries(i);

      if(n < nMin) continue;
      if(err <= 0) continue;
      if(x < fitXMin || x > fitXMax) continue;

      double yFit = fitFunc->Eval(x);

      double sigmaEff = std::max(err, 20.0);

      if(fabs(y - yFit) > nSigmaCut * sigmaEff) continue;

      gClean->SetPoint(point, x, y);
      gClean->SetPointError(point, 0., err);
      point++;

    }
    // =========================================================
    // --- 2° FIT (pulito)
    // =========================================================
    TF1* fitFunc2 = new TF1(Form("fit2_%s", canvasName.c_str()),
                           polFunc.c_str(), fitXMin, fitXMax);

    gClean->Fit(fitFunc2,"QRS");

    // =========================================================
    // --- DRAW
    // =========================================================
    fitFunc2->SetLineColor(kRed);
    fitFunc2->SetLineWidth(2);
    fitFunc2->Draw("same");

    // opzionale: mostra punti usati
    gClean->SetMarkerStyle(20);
    gClean->SetMarkerColor(kBlue);
    gClean->Draw("P same");

    // =========================================================
    // --- SAVE
    // =========================================================
    c->Print(Form("%s/%s/%s.png", plotDir.c_str(), subDir.c_str(), canvasName.c_str()));
    c->Print(Form("%s/%s/%s.pdf", plotDir.c_str(), subDir.c_str(), canvasName.c_str()));

    delete latex;
    delete c;
    delete fitFunc;   // tieni solo quello finale

    return fitFunc2;
}






#endif
