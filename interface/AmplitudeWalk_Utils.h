#ifndef AMPLITUDEWALK_UTILS_H
#define AMPLITUDEWALK_UTILS_H

#include "CfgManager/interface/CfgManager.h"
#include "CfgManager/interface/CfgManagerT.h"
#include "interface/FitUtils.h"
#include "interface/AnalysisUtils.h"

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
#include "TKey.h"

//////////////////////////////////////////////////////////////////
//--- FUNCTIONS USED IN THE STEP 1 FOR ENERGY INTERCALIBRATION  //
//////////////////////////////////////////////////////////////////
std::map<std::tuple<int,std::string,int>,double> calibMap_DUT;
std::map<std::tuple<int,std::string,int>,double> calibMap_REF;


// --- Loading energy intercalibration factors:
inline void LoadEnergyCalibrationFactors(
    const std::string& fileName,
    std::map<std::tuple<int,std::string,int>, double>& calibMap
)
{
    std::ifstream fin(fileName);

    if(!fin)
    {
        std::cout << "[ERROR] Cannot open file "<< fileName << std::endl;
        return;
    }

    int bar;
    int vth;
    std::string side;
    double coeff;

    while(fin >> bar >> side >> vth >> coeff)
    {
        calibMap[std::make_tuple(bar, side, vth)] = coeff;

        std::cout << "Loaded: "
                  << "bar=" << bar
                  << " side=" << side
                  << " vth=" << vth
                  << " coeff=" << coeff
                  << std::endl;
    }

    fin.close();
}



// --- Application of the energy interalibration:
inline float ApplyEnergyCalibration(
    float energyRaw,
    int bar,
    const std::string& side,
    int vth,
    const std::map<std::tuple<int,std::string,int>, double>& calibMap
)
{
    auto key = std::make_tuple(bar, side, vth);

    auto it = calibMap.find(key);

    if(it == calibMap.end())
    {
        std::cout << "[WARNING] Missing calibration for "
                  << "bar=" << bar
                  << " side=" << side
                  << " vth=" << vth
                  << std::endl;

        return energyRaw;
    }

    return energyRaw * it->second;
}









//////////////////////////////////////////////////////////////////
//--- FUNCTIONS USED IN THE AMPLITUDE WALK CORRECTION STUDIES:  //
//////////////////////////////////////////////////////////////////

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
// --- Function for computing the CTR From Histo
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
// --- Function latex for histrograms
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
// Compute the RMS and the mean from a TGraphErrors plot
//------------------------------------------------------------------------
inline void computeMeanRMSTGraphErrors(TGraphErrors* g, double& mean, double& rms)
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
// Compute the fit range for TProfile
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





//-----------------------------------------------------
// --- Function plotting histograms TH1F (no fit)
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
// --- Function plotting histograms TH1F (with fit)
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
    TF1* fitFunc = new TF1(Form("fit_%s", canvasName.c_str()),"gaus",histo->GetMean()-1.*histo->GetRMS(),histo->GetMean()+1.*histo->GetRMS());


    // 1° fit
    histo->Fit(fitFunc, "QNRS");

    // 2° fit (range centrato)
    histo->Fit(fitFunc, "QSR+","",fitFunc->GetParameter(1) - 1.*fitFunc->GetParameter(2),fitFunc->GetParameter(1) + 1.*fitFunc->GetParameter(2));

    // 3° fit (refinement)
    histo->Fit(fitFunc, "QSR+","",fitFunc->GetParameter(1) - 1.*fitFunc->GetParameter(2),fitFunc->GetParameter(1) + 1.*fitFunc->GetParameter(2));

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
// --- Function plotting TProfile (with fit)
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
        //float fitXMax = 550.;
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








//-----------------------------------------------------------
// --- Function plotting TH2F and TProfile in the same canvas
// ----------------------------------------------------------

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




//-------------------------------------------------------------------------
// --- Function build TGaphErrors Left and right side in the same canvas:
//-------------------------------------------------------------------------
inline void BuildGraphsErrors(
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
            gL_map[key]->SetPoint(i, bar, val);
            gL_map[key]->SetPointError(i, 0., eval);
        }
        else if(side == "R")
        {
            int i = counterR[key]++;
            gR_map[key]->SetPoint(i, bar, val);
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




//-------------------------------------------------------------------------------
// --- Function plot and save TGaphErrors Left and right side in the same canvas:
//-------------------------------------------------------------------------------
inline void PlotGraphsErrors(
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
        //gL->GetXaxis()->SetRangeUser(0., 17.);
        gL->GetXaxis()->SetLimits(-1, 16);
        gL->GetXaxis()->SetNdivisions(17, 0, 0, kTRUE);
        
	gL->SetTitle(Form(";Bar index; %s", yTitle.c_str()));

        gL->Draw("APL");
        gR->Draw("PL SAME");

        // stats
        double mean_L, rms_L;
        double mean_R, rms_R;

        computeMeanRMSTGraphErrors(gL, mean_L, rms_L);
        computeMeanRMSTGraphErrors(gR, mean_R, rms_R);

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


/*
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

        computeMeanRMSTGraphErrors(gL, mean_L, rms_L);
        computeMeanRMSTGraphErrors(gR, mean_R, rms_R);

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

*/

//-------------------------------------------------------------------------------------------
// --- helper function TGaph Left and right same:
//------------------------------------------------------------------------
inline void BuildGraphs_mono(
    const std::map<std::tuple<int,float,int>, std::pair<double,double>>& input_map,
    std::map<std::pair<float,int>, TGraphErrors*>& gL_map
)
{
    std::map<std::pair<float,int>, int> counterL;

    for(const auto& it : input_map)
    {
        int bar;
        float vov;
        int thr;
        double val;
        double eval;

        std::tie(bar, vov, thr) = it.first;
        val = it.second.first;
        eval = it.second.second;

        std::pair<float,int> key = std::make_pair(vov, thr);

        if(gL_map.find(key) == gL_map.end())
        {
            gL_map[key] = new TGraphErrors();
            counterL[key] = 0;
            
        }
            int i = counterL[key]++;
            gL_map[key]->SetPoint(i, bar+1, val);
            gL_map[key]->SetPointError(i, 0., eval);

    }
}



//------------------------------------------------------------------------
// --- helper function TGaph Left and right same:
//------------------------------------------------------------------------
inline void PlotGraphs_mono(
    const std::map<std::pair<float,int>, TGraphErrors*>& gL_map,
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

        TCanvas* c = new TCanvas(Form("c_%s_vov%.2f_th%d", label.c_str(), vov, thr),"", 800, 600);

        StyleGraph(gL, kRed, 20);

        gL->GetYaxis()->SetRangeUser(ymin, ymax);
        gL->GetXaxis()->SetRangeUser(0., 17.);

        gL->SetTitle(Form(";Bar index; %s", yTitle.c_str()));

        gL->Draw("APL");

        // stats
        double mean_L, rms_L;

        computeMeanRMSTGraphErrors(gL, mean_L, rms_L);

        TLegend* leg = new TLegend(0.20,0.70,0.60,0.9);
        leg->SetHeader(Form("V_{OV} = %.2f V, vth = %d DAC", vov, thr), "C");

        leg->AddEntry(gL, Form("Side L (Mean=%.4f, RMS=%.4f)", mean_L, rms_L), "lp");

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









//-----------------------------------------------------------------------------------------
// --- Function to plot togather different amplitude walk function from different channels:
//-----------------------------------------------------------------------------------------

inline void DrawOverlayFitsFromIndices(
    const std::map<double, TF1*>& fitMap,
    const std::vector<double>& selectedIndices,
    const std::string& canvasName,
    const std::string& outDir,
    //const std::string& subdir,
    const std::string& header = "",
    float yMin = 0.,
    float yMax = 3000.
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

        func->SetLineColor(color);
        func->SetLineWidth(2);

        //float bar = idx - 13001100;

	int idx_int = int(idx);
        int bar = idx_int % 100;

        //leg->AddEntry(func, Form("bar = %d", bar), "l");
        if(first)
        {
            func->SetTitle(";Energy [a.u.];#Delta T raw [ps]");
	    if(yMax > yMin)
                func->SetRange(func->GetXmin(), func->GetXmax());
	    
	    func->Draw();
	    if(yMax > yMin)
            {
                gPad->Update();
                func->GetHistogram()->GetYaxis()->SetRangeUser(yMin, yMax);
            }
            first = false;
        }
        else
        {
            func->Draw("same");
        }

        leg->AddEntry(func, Form("bar = %.d", bar), "l");

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



//-----------------------------------------------------------------------------------------
// --- Wrapper function ao plot togather different amplitude walk function from different channels:
//-----------------------------------------------------------------------------------------


inline void DrawOverlayFits_CorrFunction(
    const std::map<double, TF1*>& fitMap,
    const std::string& plotDir,
    const std::string& sideLabel // "L" o "R"
)
{
    // mappa: vth → lista indici
    std::map<int, std::vector<double>> vthMap;

    // ----------------------------------
    // --- Raggruppa per vth
    // ----------------------------------
    for(const auto& [idx, func] : fitMap)
    {
        if(!func) continue;

        int idx_int = int(idx);

        int bar = idx_int % 100;
        int vth = (idx_int / 100) % 100;

	if(bar < 4 || bar > 10) continue;

        vthMap[vth].push_back(idx);
    }

    // ----------------------------------
    // --- Loop sui vth
    // ----------------------------------
    for(const auto& [vth, indices] : vthMap)
    {
        if(indices.empty()) continue;

        std::string canvasName = Form("c_overlay_%s_vth%d", sideLabel.c_str(), vth);
        //std::string outDir     = plotDir + "/Loop3/TW_func_comparison/Side_" + sideLabel;
        std::string outDir     = plotDir + "_" + sideLabel;
        std::string header = Form("Side %s - vth = %d DAC", sideLabel.c_str(), vth);

        DrawOverlayFitsFromIndices(
            fitMap,
            indices,
            canvasName,
            outDir,
            header,
            -500., 500.
        );
    }
}



//-----------------------------------------------------------------------------------------
// --- Plot stessa barra, diversi VTH (una canvas per barra)
//-----------------------------------------------------------------------------------------
inline void DrawOverlayFits_SameBar_DifferentVth(
    const std::map<double, TF1*>& fitMap,
    const std::string& outDir,
    const std::string& sideLabel,   // "L" o "R"
    float yMin = 0.,
    float yMax = 0.
)
{
    // --- mappa: bar -> lista di idx
    std::map<int, std::vector<double>> barMap;

    // =========================================================
    // --- Raggruppa per barra
    // =========================================================
    for(const auto& [idx, func] : fitMap)
    {
        if(!func) continue;

        int idx_int = int(idx);

        int bar = idx_int % 100;
        int vth = (idx_int / 100) % 100;

	if (vth == 40) continue;
        barMap[bar].push_back(idx);
    }

    // =========================================================
    // --- Loop sulle barre
    // =========================================================
    for(const auto& [bar, indices] : barMap)
    {
        TCanvas* c = new TCanvas(
            Form("c_bar%02d_%s", bar, sideLabel.c_str()),
            Form("Bar %02d - %s", bar, sideLabel.c_str()),
            800,600
        );

        TLegend* leg = new TLegend(0.55,0.5,0.85,0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextSize(0.03);

        int color = 1;
        bool first = true;

        // --- ordina per vth (importante!)
        std::vector<std::pair<int,double>> sorted;
        for(double idx : indices)
        {
            int idx_int = int(idx);
            int vth = (idx_int / 100) % 100;
            sorted.push_back({vth, idx});
        }

        std::sort(sorted.begin(), sorted.end());

        // =========================================================
        // --- Draw
        // =========================================================
        for(const auto& [vth, idx] : sorted)
        {
            TF1* func = fitMap.at(idx);
            if(!func) continue;

            func->SetLineColor(color);
            func->SetLineWidth(2);

            if(first)
            {
                func->SetTitle(";Energy [a.u.];#Delta T raw [ps]");
                func->Draw();

                if(yMax > yMin)
                {
                    gPad->Update();
                    func->GetHistogram()->GetYaxis()->SetRangeUser(yMin, yMax);
                }

                first = false;
            }
            else
            {
                func->Draw("same");
            }

            leg->AddEntry(func, Form("vth = %d", vth), "l");

            color++;
            if(color == 10) color = 1;
        }

        // header
        leg->SetHeader(Form("Bar %02d - Side %s", bar, sideLabel.c_str()), "C");

        leg->Draw();
        c->SetGrid();

        // =========================================================
        // --- Save
        // =========================================================
        c->Print(Form("%s/c_bar%02d_%s.pdf", outDir.c_str(), bar, sideLabel.c_str()));
        c->Print(Form("%s/c_bar%02d_%s.png", outDir.c_str(), bar, sideLabel.c_str()));

        delete leg;
        delete c;
    }
}



//---------------------------------------------------------------------------------
// --- Plot fit function almplitude walk for the two channels of the REF bar
//---------------------------------------------------------------------------------

inline void DrawOverlayTW_differentCh_REF(
    std::map<int,TF1*>& fitMapL,
    std::map<int,TF1*>& fitMapR,
    const std::string& plotDir
)
{
    for (auto& mapL : fitMapL)
    {
        int index_ref = mapL.first;

        if (fitMapR.find(index_ref) == fitMapR.end()) continue;

        TF1* fL = mapL.second;
        TF1* fR = fitMapR[index_ref];

        if (!fL || !fR) continue;

        int vov_int = index_ref / 10000;
        float Vov = vov_int / 100.;

        int vth = (index_ref % 10000) / 100;

        TCanvas* c = new TCanvas(
            Form("c_overlay_REF_LR_Vov%.2f_th%02d",Vov,vth),
            Form("c_overlay_REF_LR_Vov%.2f_th%02d",Vov,vth),
            800,600
        );

        fL->SetLineColor(kBlue);
        fL->SetLineWidth(2);

        fR->SetLineColor(kRed);
        fR->SetLineWidth(2);

        fL->SetTitle(
            Form("REF TW correction - Vov %.2f th %d;Energy [a.u.];TW correction [ps]",
            Vov,vth)
        );

        fL->Draw();
        fR->Draw("same");

        TLegend* leg = new TLegend(0.60,0.75,0.88,0.88);
        leg->AddEntry(fL,"REF channel L","l");
        leg->AddEntry(fR,"REF channel R","l");
        leg->Draw();

        c->SaveAs(Form("%s/Overlay_REF_LR_Vov%.2f_th%02d.png",
                       plotDir.c_str(),Vov,vth));

        c->SaveAs(Form("%s/Overlay_REF_LR_Vov%.2f_th%02d.pdf",
                       plotDir.c_str(),Vov,vth));

        delete c;
    }
}








//--------------------------------------------------------------------------------------------------
// --- Plot fit function amplitude walk for the same channel of the REF bar and different vth values
//--------------------------------------------------------------------------------------------------
inline void DrawOverlayTW_same_ch_different_Vth_REF(
    std::map<int,TF1*>& fitMap,
    const std::string& plotDir,
    const std::string& channelLabel
)
{
    TCanvas* c = new TCanvas(
        Form("c_overlay_REF_%s_vsVth",channelLabel.c_str()),
        Form("c_overlay_REF_%s_vsVth",channelLabel.c_str()),
        900,700
    );

    TLegend* leg = new TLegend(0.60,0.60,0.88,0.88);

    bool first = true;

    int color = 1;

    for (auto& mapIt : fitMap)
    {
        int index_ref = mapIt.first;

        TF1* func = mapIt.second;

        if (!func) continue;

        int vov_int = index_ref / 10000;
        float Vov = vov_int / 100.;

        int vth = (index_ref % 10000) / 100;

        func->SetLineColor(color);
        func->SetLineWidth(2);

        func->SetTitle(
            Form("REF channel %s TW correction;Energy [a.u.];TW correction [ps]",
            channelLabel.c_str())
        );

        if (first)
        {
            func->Draw();
            first = false;
        }
        else
        {
            func->Draw("same");
        }

        leg->AddEntry(func,
            Form("th %d",vth),
            "l");

        color++;

        if (color == 5) color++;
    }

    leg->Draw();

    c->SaveAs(Form("%s/Overlay_REF_%s_vsVth.png",
                   plotDir.c_str(),channelLabel.c_str()));

    c->SaveAs(Form("%s/Overlay_REF_%s_vsVth.pdf",
                   plotDir.c_str(),channelLabel.c_str()));

    delete c;
}







//------------------------------------------------------------------------
// --- draw and fit of TProfile with a (maybe) more robust fit procedure:
//------------------------------------------------------------------------



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



#include <utility> // per std::pair
#include <cmath>



//--------------------------------------------------------------------------------------
// --- From a TProfile, determination of the mean value of the points in a certain range:
//     In this case we are considering the simple aritmetic mean not considering the weigth 
//     of the statistic per bin
//--------------------------------------------------------------------------------------

inline std::pair<double,double> ComputeProfileMeanInRange(
    TProfile* prof,
    double xMin,
    double xMax,
    int nMinEntries = 5   // filtro statistico per bin
)
{
    if (!prof) return std::make_pair(0., 0.);

    double sumY = 0.;
    double sumY2 = 0.;
    int nPoints = 0;

    for (int i = 1; i <= prof->GetNbinsX(); ++i)
    {
        double x = prof->GetBinCenter(i);
        if (x < xMin || x > xMax) continue;

        double y   = prof->GetBinContent(i);
        double err = prof->GetBinError(i);
        double n   = prof->GetBinEntries(i);

        if (n < nMinEntries) continue;
        if (err <= 0) continue;

        sumY  += y;
        sumY2 += y*y;
        nPoints++;
    }

    if (nPoints == 0) return std::make_pair(0., 0.);

    double mean = sumY / nPoints;

    // RMS dei punti
    double variance = (sumY2 / nPoints) - (mean * mean);
    double rms = (variance > 0) ? std::sqrt(variance) : 0.;

    // errore sulla media
    double errMean = rms / std::sqrt(nPoints);

    return std::make_pair(mean, errMean);
}



//---------------------------------------------------------------------------------------
// --- determiation of the aritmetic mean and weithed mean for a sample of TProfile bin :
//--------------------------------------------------------------------------------------


inline std::tuple<double,double,double> ComputeProfileMeansComparison(
    TProfile* prof,
    double xMin,
    double xMax,
    int nMinEntries = 5
)
{
    if (!prof) return std::make_tuple(0.,0.,0.);

    // --- media semplice
    double sumY = 0.;
    double sumY2 = 0.;
    int nPoints = 0;

    // --- media pesata (peso = nEntries)
    double sumWeightedY = 0.;
    double sumWeights   = 0.;

    for (int i = 1; i <= prof->GetNbinsX(); ++i)
    {
        double x = prof->GetBinCenter(i);
        if (x < xMin || x > xMax) continue;

        double y   = prof->GetBinContent(i);
        double err = prof->GetBinError(i);
        double n   = prof->GetBinEntries(i);

        if (n < nMinEntries) continue;
        if (err <= 0) continue;
        if (std::isnan(y) || std::isinf(y)) continue;

        // --- media semplice
        sumY  += y;
        sumY2 += y*y;
        nPoints++;

        // --- media pesata
        sumWeightedY += y * n;
        sumWeights   += n;
    }

    if (nPoints == 0 || sumWeights == 0)
        return std::make_tuple(0.,0.,0.);

    // ============================
    // --- media semplice
    // ============================
    double mean_simple = sumY / nPoints;

    double variance = (sumY2 / nPoints) - (mean_simple * mean_simple);
    double rms = (variance > 0) ? std::sqrt(variance) : 0.;
    double errMean_simple = rms / std::sqrt(nPoints);

    // ============================
    // --- media pesata
    // ============================
    double mean_weighted = sumWeightedY / sumWeights;

    // (errore pesato opzionale, qui non lo calcolo per semplicità)

    // ============================
    // --- differenza
    // ============================
    double diff = mean_weighted - mean_simple;

    return std::make_tuple(mean_simple, mean_weighted, diff);
}











//------------------------------------------------------------
// --- Save diff(mean_simple - mean_weighted) per vth
//------------------------------------------------------------
inline void SaveProfileComparisonToTxt(
    const std::map<std::tuple<int,float,int,std::string>, std::pair<double,double>>& inputMap,
    const std::string& outDir
)
{
    // mappa: vth -> ofstream*
    std::map<int, std::ofstream*> fileMap;

    for (const auto& it : inputMap)
    {
        int bar;
        float vov;
        int vth;
        std::string side;

        std::tie(bar, vov, vth, side) = it.first;

        double value = it.second.first;
	double valueErr = it.second.second;

        // --- crea file se non esiste
        if (fileMap.find(vth) == fileMap.end())
        {
            std::string fileName = Form("%s/profileComparison_vth%d.txt", outDir.c_str(), vth);
            fileMap[vth] = new std::ofstream(fileName);

            // header
            (*fileMap[vth]) << "# bar  side  Value ValueErr\n";
        }

        // --- scrivi riga
        (*fileMap[vth]) << bar << "  "
                        << side << "  "
                        << value << "  "
	                << valueErr << "\n";
    }

    // --- chiudi tutti i file
    for (auto& it : fileMap)
    {
        it.second->close();
        delete it.second;
    }
}




//-----------------------------------------------------------------------
// --- Save txt File from map <int,float,int> -> <double,double>  per vth
//-----------------------------------------------------------------------
inline void SaveProfileComparisonToTxt_noSide(
    const std::map<std::tuple<int,float,int>, std::pair<double,double>>& inputMap,
    const std::string& outDir
)
{
    // mappa: vth -> ofstream*
    std::map<int, std::ofstream*> fileMap;

    for (const auto& it : inputMap)
    {
        int bar;
        float vov;
        int vth;

        std::tie(bar, vov, vth) = it.first;

        double value = it.second.first;
        double valueErr = it.second.second;

        // --- crea file se non esiste
        if (fileMap.find(vth) == fileMap.end())
        {
            std::string fileName = Form("%s/outfileFile_vth%d.txt", outDir.c_str(), vth);
            fileMap[vth] = new std::ofstream(fileName);

            // header
            (*fileMap[vth]) << "# bar  side  diff_mean\n";
        }

        // --- scrivi riga
        (*fileMap[vth]) << bar << "  "
                        << value << "  "
                        << valueErr << "\n";
    }

    // --- chiudi tutti i file
    for (auto& it : fileMap)
    {
        it.second->close();
        delete it.second;
    }
}



//-----------------------------------------------------------------------------------
// --- Save T1F function in a root file:
//-----------------------------------------------------------------------------------

inline void SaveTF1Map(
    TFile* outFile,
    const std::map<int,TF1*>& funcMap,
    const std::string& funcLabel
)
{
    outFile->cd();

    for(const auto& map_func : funcMap)
    {
        int index_ref = map_func.first;

        int vov_int = index_ref / 10000;
        float Vov = vov_int / 100.;

        int vth = (index_ref % 10000) / 100;

        std::string dirName =
            Form("Vov%.2f_th%02d",Vov,vth);

        if( !outFile->GetDirectory(dirName.c_str()) )
            outFile->mkdir(dirName.c_str());

        outFile->cd(dirName.c_str());

        TF1* f = (TF1*)(map_func.second->Clone());

        f->SetName(
            Form("%s_Vov%.2f_th%02d",
            funcLabel.c_str(),
            Vov,
            vth)
        );

        f->Write("",TObject::kOverwrite);
	delete f;
    }
}



//---------------------------------------------------------------------------------------
// --- Load the T1F function saved from a root file and re-build the map [index_ref]->TF1
//---------------------------------------------------------------------------------------

inline std::map<int,TF1*> LoadTF1Map(
    const std::string& fileName,
    const std::string& funcLabel
)
{
    std::map<int,TF1*> funcMap;

    TFile* inFile = TFile::Open(fileName.c_str());

    if(!inFile || inFile->IsZombie())
    {
        std::cout << "ERROR: cannot open file "
                  << fileName << std::endl;

        return funcMap;
    }

    TIter nextDir(inFile->GetListOfKeys());
    TKey* keyDir;

    while( (keyDir = (TKey*)nextDir()) )
    {
        TObject* obj = keyDir->ReadObj();

        if( !obj->InheritsFrom("TDirectory") )
            continue;

        TDirectory* dir = (TDirectory*)obj;

        std::string dirName = dir->GetName();

        float Vov;
        int vth;

        sscanf(
            dirName.c_str(),
            "Vov%f_th%d",
            &Vov,
            &vth
        );

        int index_ref =
            (10000*int(Vov*100.))
            + (100*vth);

        std::string funcName =
            Form("%s_Vov%.2f_th%02d",
            funcLabel.c_str(),
            Vov,
            vth);

        TF1* f = (TF1*)dir->Get(funcName.c_str());

        if(!f)
        {
            std::cout << "WARNING: function "
                      << funcName
                      << " not found"
                      << std::endl;

            continue;
        }

        funcMap[index_ref] =
            (TF1*)f->Clone();
    }

    inFile->Close();

    return funcMap;
}





//-------------------------------------------------------------------------------------------------------------------------
// --- Funzione che prende i valori della Sigma(Delta T) Raw e Sigma(Delta T) Corrected per i due canali della barra di REF
//     e ne fa un TGraphErros in funzione del valore di vth.
//------------------------------------------------------------------------------------------------------------------------

inline void DrawSigmaComparisonVsVth_REF(
    std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> >& sigmaRawMap,
    std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> >& sigmaCorrMap,
    const std::string& plotDir,
    const std::string& outName,
    int iBar,
    float Vov,
    float ymin,
    float ymax
)
{
    TGraphErrors* gL_raw  = new TGraphErrors();
    TGraphErrors* gL_corr = new TGraphErrors();

    TGraphErrors* gR_raw  = new TGraphErrors();
    TGraphErrors* gR_corr = new TGraphErrors();

    int pL_raw  = 0;
    int pL_corr = 0;
    int pR_raw  = 0;
    int pR_corr = 0;

    // ---------------------------------------------------
    // --- RAW
    // ---------------------------------------------------

    for(auto mapIt : sigmaRawMap)
    {
        int bar          = std::get<0>(mapIt.first);
        float Vov_map    = std::get<1>(mapIt.first);
        int vth          = std::get<2>(mapIt.first);
        std::string side = std::get<3>(mapIt.first);

        if(bar != iBar) continue;
        if(std::abs(Vov_map - Vov) > 0.001) continue;

        double sigma    = mapIt.second.first;
        double sigmaErr = mapIt.second.second;

        if(side == "L")
        {
            gL_raw->SetPoint(pL_raw,vth,sigma);
            gL_raw->SetPointError(pL_raw,0,sigmaErr);
            ++pL_raw;
        }

        if(side == "R")
        {
            gR_raw->SetPoint(pR_raw,vth,sigma);
            gR_raw->SetPointError(pR_raw,0,sigmaErr);
            ++pR_raw;
        }
    }

    // ---------------------------------------------------
    // --- CORRECTED
    // ---------------------------------------------------

    for(auto mapIt : sigmaCorrMap)
    {
        int bar          = std::get<0>(mapIt.first);
        float Vov_map    = std::get<1>(mapIt.first);
        int vth          = std::get<2>(mapIt.first);
        std::string side = std::get<3>(mapIt.first);

        if(bar != iBar) continue;
        if(std::abs(Vov_map - Vov) > 0.001) continue;

        double sigma    = mapIt.second.first;
        double sigmaErr = mapIt.second.second;

        if(side == "L")
        {
            gL_corr->SetPoint(pL_corr,vth,sigma);
            gL_corr->SetPointError(pL_corr,0,sigmaErr);
            ++pL_corr;
        }

        if(side == "R")
        {
            gR_corr->SetPoint(pR_corr,vth,sigma);
            gR_corr->SetPointError(pR_corr,0,sigmaErr);
            ++pR_corr;
        }
    }

    // ---------------------------------------------------
    // --- STYLE
    // ---------------------------------------------------

    // --- L raw
    gL_raw->SetMarkerStyle(24); 
    gL_raw->SetMarkerColor(kRed);
    gL_raw->SetLineColor(kRed);
    gL_raw->SetLineWidth(2);
    gL_raw->SetLineStyle(2);

    // --- L corrected
    gL_corr->SetMarkerStyle(21); 
    gL_corr->SetMarkerColor(kRed);
    gL_corr->SetLineColor(kRed);
    gL_corr->SetLineWidth(2);

    // --- R raw
    gR_raw->SetMarkerStyle(25); 
    gR_raw->SetMarkerColor(kBlue);
    gR_raw->SetLineColor(kBlue);
    gR_raw->SetLineWidth(2);
    gR_raw->SetLineStyle(2);

    // --- R corrected
    gR_corr->SetMarkerStyle(20); 
    gR_corr->SetMarkerColor(kBlue);
    gR_corr->SetLineColor(kBlue);
    gR_corr->SetLineWidth(2);

    // ---------------------------------------------------
    // --- DRAW
    // ---------------------------------------------------

    TCanvas* c = new TCanvas(Form("c_%s",outName.c_str()),Form("c_%s",outName.c_str()),1200,700);

    gL_raw->SetTitle(";Threshold [DAC];#sigma(#DeltaT) [ps]");

    gL_raw->GetYaxis()->SetRangeUser(ymin, ymax);

    gL_raw->Draw("APL");
    gL_corr->Draw("PL SAME");

    gR_raw->Draw("PL SAME");
    gR_corr->Draw("PL SAME");

    // ---------------------------------------------------
    // --- LEGEND
    // ---------------------------------------------------
    double mean_L_raw, rms_L_raw, mean_R_raw, rms_R_raw;
    double mean_L_corr, rms_L_corr, mean_R_corr, rms_R_corr;
    computeMeanRMSTGraphErrors(gL_raw, mean_L_raw, rms_L_raw);
    computeMeanRMSTGraphErrors(gR_raw, mean_R_raw, rms_R_raw);
    computeMeanRMSTGraphErrors(gL_corr, mean_L_corr, rms_L_corr);
    computeMeanRMSTGraphErrors(gR_corr, mean_R_corr, rms_R_corr);

    TLegend* leg = new TLegend(0.35,0.66,0.77,0.94);
    leg->SetHeader("V_{OV} = 3.0 V, REF bar = 07 ", "C");

    leg->AddEntry(gL_raw ,Form("#sigma(#Delta T) raw - L (Mean=%.4f, RMS=%.4f)", mean_L_raw, rms_L_raw),"pl");
    leg->AddEntry(gR_raw ,Form("#sigma(#Delta T) raw - R (Mean=%.4f, RMS=%.4f)", mean_R_raw, rms_R_raw),"pl");
    leg->AddEntry(gL_corr,Form("#sigma(#Delta T) corr - L (Mean=%.4f, RMS=%.4f)", mean_L_corr, rms_L_corr),"pl");
    leg->AddEntry(gR_corr,Form("#sigma(#Delta T) corr - R (Mean=%.4f, RMS=%.4f)", mean_R_corr, rms_R_corr),"pl");

    leg->SetBorderSize(1);
    leg->Draw();

    c->SetGridy();
    c->SetTickx();
    c->SetTicky();

    c->SaveAs(Form("%s/%s.png",plotDir.c_str(),outName.c_str()));
    c->SaveAs(Form("%s/%s.pdf",plotDir.c_str(),outName.c_str()));

    delete leg;
    delete c;
}






//-------------------------------------------------------------------------------------------------
// --- Prende due mappe di sigma Delta T e ne fa il plot combinato nella stessa canvas vs bar index
//-------------------------------------------------------------------------------------------------
inline void DrawSigmaComparisonVsBar(
    std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> >& SigmaRawRaw,
    std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> >& SigmaRaw,
    const std::vector<int>& barList,
    const std::vector<std::string>& stepLabels,
    std::map<std::string,float>& map_Vovs,
    std::map<std::string,float>& map_ths,
    const std::string& plotDir,
    const std::string& outName
)
{
    for(auto stepLabel : stepLabels)
    {
        float Vov = map_Vovs[stepLabel];
        int vth   = map_ths[stepLabel];

        TGraphErrors* gL_raw_raw = new TGraphErrors();
        TGraphErrors* gL_raw    = new TGraphErrors();

        TGraphErrors* gR_raw_raw = new TGraphErrors();
        TGraphErrors* gR_raw    = new TGraphErrors();

        int pL1=0, pL2=0, pR1=0, pR2=0;

        double ymin = 999999.;
        double ymax = -999999.;

        // ---------------------------------------------------
        // --- LOOP OVER BARS
        // ---------------------------------------------------

        for(auto iBar : barList)
        {
            auto keyL = std::make_tuple(iBar,Vov,vth,"L");
            auto keyR = std::make_tuple(iBar,Vov,vth,"R");

            // ---------------- L raw raw ----------------
            if( SigmaRawRaw.count(keyL) )
            {
                double sigma = SigmaRawRaw[keyL].first;
                double err   = SigmaRawRaw[keyL].second;

                gL_raw_raw->SetPoint(pL1, iBar, sigma);
                gL_raw_raw->SetPointError(pL1, 0, err);

                ymin = std::min(ymin, sigma);
                ymax = std::max(ymax, sigma);

                ++pL1;
            }

            // ---------------- L raw ----------------
            if( SigmaRaw.count(keyL) )
            {
                double sigma = SigmaRaw[keyL].first;
                double err   = SigmaRaw[keyL].second;

                gL_raw->SetPoint(pL2, iBar, sigma);
                gL_raw->SetPointError(pL2, 0, err);

                ymin = std::min(ymin, sigma);
                ymax = std::max(ymax, sigma);

                ++pL2;
            }

            // ---------------- R raw raw ----------------
            if( SigmaRawRaw.count(keyR) )
            {
                double sigma = SigmaRawRaw[keyR].first;
                double err   = SigmaRawRaw[keyR].second;

                gR_raw_raw->SetPoint(pR1, iBar, sigma);
                gR_raw_raw->SetPointError(pR1, 0, err);

                ymin = std::min(ymin, sigma);
                ymax = std::max(ymax, sigma);

                ++pR1;
            }

            // ---------------- R raw ----------------
            if( SigmaRaw.count(keyR) )
            {
                double sigma = SigmaRaw[keyR].first;
                double err   = SigmaRaw[keyR].second;

                gR_raw->SetPoint(pR2, iBar, sigma);
                gR_raw->SetPointError(pR2, 0, err);

                ymin = std::min(ymin, sigma);
                ymax = std::max(ymax, sigma);

                ++pR2;
            }
        }

        ymin *= 0.90;
        ymax *= 1.10;

        // ---------------------------------------------------
        // --- STYLE
        // ---------------------------------------------------

        // --- L raw raw
        gL_raw_raw->SetMarkerStyle(24);
        gL_raw_raw->SetMarkerColor(kRed);
        gL_raw_raw->SetLineColor(kRed);
        gL_raw_raw->SetLineWidth(2);
        gL_raw_raw->SetLineStyle(2);

        // --- L raw
        gL_raw->SetMarkerStyle(21);
        gL_raw->SetMarkerColor(kRed);
        gL_raw->SetLineColor(kRed);
        gL_raw->SetLineWidth(2);

        // --- R raw raw
        gR_raw_raw->SetMarkerStyle(25);
        gR_raw_raw->SetMarkerColor(kBlue);
        gR_raw_raw->SetLineColor(kBlue);
        gR_raw_raw->SetLineWidth(2);
        gR_raw_raw->SetLineStyle(2);

        // --- R raw
        gR_raw->SetMarkerStyle(20);
        gR_raw->SetMarkerColor(kBlue);
        gR_raw->SetLineColor(kBlue);
        gR_raw->SetLineWidth(2);

        // ---------------------------------------------------
        // --- DRAW
        // ---------------------------------------------------

        TCanvas* c = new TCanvas(
            Form("c_%s_Vov%.2f_th%02d",outName.c_str(),Vov,vth),
            Form("c_%s_Vov%.2f_th%02d",outName.c_str(),Vov,vth),
            1200,700
        );

        gL_raw_raw->SetTitle(
            Form("V_{OV} = %.2f V, threshold = %d DAC;Bar index;#sigma(#DeltaT) [ps]",
            Vov,vth)
        );

        gL_raw_raw->GetYaxis()->SetRangeUser(50.,130);
	gL_raw_raw->GetXaxis()->SetRangeUser(-1.,16.);
        gL_raw_raw->GetXaxis()->SetNdivisions(17, 0, 0, kTRUE);

        gL_raw_raw->Draw("APL");
        gL_raw->Draw("PL SAME");

        gR_raw_raw->Draw("PL SAME");
        gR_raw->Draw("PL SAME");

        // ---------------------------------------------------
        // --- LEGEND
        // ---------------------------------------------------
        double mean_L_raw, rms_L_raw, mean_R_raw, rms_R_raw;
        double mean_L_raw_raw, rms_L_raw_raw, mean_R_raw_raw, rms_R_raw_raw;
        computeMeanRMSTGraphErrors(gL_raw, mean_L_raw, rms_L_raw);
        computeMeanRMSTGraphErrors(gR_raw, mean_R_raw, rms_R_raw);
        computeMeanRMSTGraphErrors(gL_raw_raw, mean_L_raw_raw, rms_L_raw_raw);
        computeMeanRMSTGraphErrors(gR_raw_raw, mean_R_raw_raw, rms_R_raw_raw);

        TLegend* leg = new TLegend(0.35,0.66,0.77,0.94);
        leg->SetHeader(Form("V_{OV} = %.2f V | vth = %d DAC", Vov, vth), "C");

        leg->AddEntry(gL_raw ,Form("#sigma(#Delta T) raw - L (Mean=%.4f, RMS=%.4f)", mean_L_raw, rms_L_raw),"pl");
        leg->AddEntry(gR_raw ,Form("#sigma(#Delta T) raw - R (Mean=%.4f, RMS=%.4f)", mean_R_raw, rms_R_raw),"pl");
        leg->AddEntry(gL_raw_raw,Form("#sigma(#Delta T) raw raw - L (Mean=%.4f, RMS=%.4f)", mean_L_raw_raw, rms_L_raw_raw),"pl");
        leg->AddEntry(gR_raw_raw,Form("#sigma(#Delta T) raw raw - R (Mean=%.4f, RMS=%.4f)", mean_R_raw_raw, rms_R_raw_raw),"pl");

        leg->SetBorderSize(1);
        leg->Draw();

        c->SetGridx();
        c->SetGridy();
	c->SetTickx();
	c->SetTicky();

        c->Print(Form("%s/%s_Vov%.2f_th%02d.png",plotDir.c_str(),outName.c_str(),Vov,vth));
        c->Print(Form("%s/%s_Vov%.2f_th%02d.pdf",plotDir.c_str(),outName.c_str(),Vov,vth));

        delete leg;
        delete c;
    }
}






//----------------------------------------------------------------------------------
// --- Comparison sigma delta T (raw, Ene corr, Ene+phase corr) vs bar index
//----------------------------------------------------------------------------------
inline void DrawSigmaEvolutionVsBar(
    std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> >& SigmaRaw,
    std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> >& SigmaEneCorr,
    std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> >& SigmaEnePhaseCorr,
    const std::vector<int>& barList,
    const std::vector<std::string>& stepLabels,
    std::map<std::string,float>& map_Vovs,
    std::map<std::string,float>& map_ths,
    const std::string& plotDir,
    const std::string& outName
)
{
    for(auto stepLabel : stepLabels)
    {
        float Vov = map_Vovs[stepLabel];
        int vth   = map_ths[stepLabel];

        for(auto side : {"L","R"})
        {
            TGraphErrors* gRaw           = new TGraphErrors();
            TGraphErrors* gEneCorr       = new TGraphErrors();
            TGraphErrors* gEnePhaseCorr  = new TGraphErrors();

            int p0=0, p1=0, p2=0;

            double ymin = 999999.;
            double ymax = -999999.;

            // ---------------------------------------------------
            // --- LOOP OVER BARS
            // ---------------------------------------------------

            for(auto iBar : barList)
            {
                auto key = std::make_tuple(iBar,Vov,vth,std::string(side));

                // ---------------- RAW ----------------

                if( SigmaRaw.count(key) )
                {
                    double sigma = SigmaRaw[key].first;
                    double err   = SigmaRaw[key].second;

                    gRaw->SetPoint(p0, iBar, sigma);
                    gRaw->SetPointError(p0, 0, err);

                    ymin = std::min(ymin, sigma);
                    ymax = std::max(ymax, sigma);

                    ++p0;
                }

                // ---------------- ENERGY CORR ----------------

                if( SigmaEneCorr.count(key) )
                {
                    double sigma = SigmaEneCorr[key].first;
                    double err   = SigmaEneCorr[key].second;

                    gEneCorr->SetPoint(p1, iBar, sigma);
                    gEneCorr->SetPointError(p1, 0, err);

                    ymin = std::min(ymin, sigma);
                    ymax = std::max(ymax, sigma);

                    ++p1;
                }

                // ---------------- ENERGY + PHASE CORR ----------------

                if( SigmaEnePhaseCorr.count(key) )
                {
                    double sigma = SigmaEnePhaseCorr[key].first;
                    double err   = SigmaEnePhaseCorr[key].second;

                    gEnePhaseCorr->SetPoint(p2, iBar, sigma);
                    gEnePhaseCorr->SetPointError(p2, 0, err);

                    ymin = std::min(ymin, sigma);
                    ymax = std::max(ymax, sigma);

                    ++p2;
                }
            }

            ymin *= 0.90;
            ymax *= 1.10;

            // ---------------------------------------------------
            // --- STYLE
            // ---------------------------------------------------

            // --- RAW
            gRaw->SetMarkerStyle(24);
            gRaw->SetMarkerColor(kBlack);
            gRaw->SetLineColor(kBlack);
            gRaw->SetLineWidth(2);
            gRaw->SetLineStyle(2);

            // --- ENERGY CORR
            gEneCorr->SetMarkerStyle(20);
            gEneCorr->SetMarkerColor(kBlue);
            gEneCorr->SetLineColor(kBlue);
            gEneCorr->SetLineWidth(2);

            // --- ENERGY + PHASE CORR
            gEnePhaseCorr->SetMarkerStyle(21);
            gEnePhaseCorr->SetMarkerColor(kRed);
            gEnePhaseCorr->SetLineColor(kRed);
            gEnePhaseCorr->SetLineWidth(2);

            // ---------------------------------------------------
            // --- DRAW
            // ---------------------------------------------------
            TCanvas* c = new TCanvas(Form("c_%s_%s_Vov%.2f_th%02d",outName.c_str(),side,Vov,vth),Form("c_%s_%s_Vov%.2f_th%02d",outName.c_str(),side,Vov,vth),1200,700);

            gRaw->SetTitle(Form("Side %s - V_{OV}=%.2f V - threshold=%d DAC;Bar index;#sigma(#DeltaT) [ps]",side,Vov,vth));
            gRaw->GetYaxis()->SetRangeUser(40.,100.);
	    gRaw->GetXaxis()->SetRangeUser(-1.,16.);
    	    gRaw->GetXaxis()->SetNdivisions(17, 0, 0, kTRUE);

            gRaw->Draw("APL");
            gEneCorr->Draw("PL SAME");
            gEnePhaseCorr->Draw("PL SAME");

            // ---------------------------------------------------
            // --- LEGEND
            // ---------------------------------------------------
            double mean_raw, rms_raw;
            double mean_EneCorr, rms_EneCorr;
	    double mean_EnePhaseCorr, rms_EnePhaseCorr;
            computeMeanRMSTGraphErrors(gRaw, mean_raw, rms_raw);
            computeMeanRMSTGraphErrors(gEneCorr, mean_EneCorr, rms_EneCorr);
	    computeMeanRMSTGraphErrors(gEnePhaseCorr, mean_EnePhaseCorr, rms_EnePhaseCorr);

            TLegend* leg = new TLegend(0.35,0.66,0.77,0.94);
            leg->SetHeader(Form("V_{OV} = %.2f V | vth = %d DAC | Side %s", Vov, vth, side), "C");

            leg->AddEntry(gRaw ,Form("#sigma(#Delta T) Raw (Mean=%.4f, RMS=%.4f)", mean_raw, rms_raw),"pl");
            leg->AddEntry(gEneCorr ,Form("#sigma(#Delta T) TW corr (Mean=%.4f, RMS=%.4f)", mean_EneCorr, rms_EneCorr),"pl");
            leg->AddEntry(gEnePhaseCorr,Form("#sigma(#Delta T) TW+Phase corr (Mean=%.4f, RMS=%.4f)", mean_EnePhaseCorr, rms_EnePhaseCorr),"pl");

            leg->SetBorderSize(1);
            leg->Draw();

            c->SetGridx();
            c->SetGridy();
	    c->SetTickx();
	    c->SetTicky();

            c->Print(Form("%s/%s_%s_Vov%.2f_th%02d.png",plotDir.c_str(),outName.c_str(),side,Vov,vth));
            c->Print(Form("%s/%s_%s_Vov%.2f_th%02d.pdf",plotDir.c_str(),outName.c_str(),side,Vov,vth));

            delete leg;
            delete c;
        }
    }
}







//----------------------------------------------------------------------------------
// --- Comparison sigma delta T (raw, Ene corr, Ene+phase corr) same channels vs vth
//----------------------------------------------------------------------------------
inline void DrawSigmaEvolutionVsThreshold(
    std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> >& SigmaRaw,
    std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> >& SigmaEneCorr,
    std::map< std::tuple<int,float,int,std::string>, std::pair<double,double> >& SigmaEnePhaseCorr,
    const std::vector<int>& barList,
    const std::vector<std::string>& stepLabels,
    std::map<std::string,float>& map_Vovs,
    std::map<std::string,float>& map_ths,
    const std::string& plotDir,
    const std::string& outName
)
{
    for(auto side : {"L","R"})
    {
        for(auto iBar : barList)
        {
            TGraphErrors* gRaw          = new TGraphErrors();
            TGraphErrors* gEneCorr      = new TGraphErrors();
            TGraphErrors* gEnePhaseCorr = new TGraphErrors();

            int p0=0, p1=0, p2=0;

            double ymin = 999999.;
            double ymax = -999999.;

            // ---------------------------------------------------
            // --- LOOP OVER THRESHOLDS
            // ---------------------------------------------------

            for(auto stepLabel : stepLabels)
            {
                float Vov = map_Vovs[stepLabel];
                int vth   = map_ths[stepLabel];

                auto key = std::make_tuple(iBar,Vov,vth,std::string(side));

                // ---------------- RAW ----------------

                if( SigmaRaw.count(key) )
                {
                    double sigma = SigmaRaw[key].first;
                    double err   = SigmaRaw[key].second;

                    gRaw->SetPoint(p0, vth, sigma);
                    gRaw->SetPointError(p0, 0, err);

                    ymin = std::min(ymin, sigma);
                    ymax = std::max(ymax, sigma);

                    ++p0;
                }

                // ---------------- ENERGY CORR ----------------

                if( SigmaEneCorr.count(key) )
                {
                    double sigma = SigmaEneCorr[key].first;
                    double err   = SigmaEneCorr[key].second;

                    gEneCorr->SetPoint(p1, vth, sigma);
                    gEneCorr->SetPointError(p1, 0, err);

                    ymin = std::min(ymin, sigma);
                    ymax = std::max(ymax, sigma);

                    ++p1;
                }

                // ---------------- ENERGY + PHASE CORR ----------------

                if( SigmaEnePhaseCorr.count(key) )
                {
                    double sigma = SigmaEnePhaseCorr[key].first;
                    double err   = SigmaEnePhaseCorr[key].second;

                    gEnePhaseCorr->SetPoint(p2, vth, sigma);
                    gEnePhaseCorr->SetPointError(p2, 0, err);

                    ymin = std::min(ymin, sigma);
                    ymax = std::max(ymax, sigma);

                    ++p2;
                }
            }

            ymin *= 0.90;
            ymax *= 1.10;

            // ---------------------------------------------------
            // --- STYLE
            // ---------------------------------------------------

            // --- RAW
            gRaw->SetMarkerStyle(24);
            gRaw->SetMarkerColor(kBlack);
            gRaw->SetLineColor(kBlack);
            gRaw->SetLineWidth(2);
            gRaw->SetLineStyle(2);

            // --- ENERGY CORR
            gEneCorr->SetMarkerStyle(20);
            gEneCorr->SetMarkerColor(kBlue);
            gEneCorr->SetLineColor(kBlue);
            gEneCorr->SetLineWidth(2);

            // --- ENERGY + PHASE CORR
            gEnePhaseCorr->SetMarkerStyle(21);
            gEnePhaseCorr->SetMarkerColor(kRed);
            gEnePhaseCorr->SetLineColor(kRed);
            gEnePhaseCorr->SetLineWidth(2);

            // ---------------------------------------------------
            // --- DRAW
            // ---------------------------------------------------

            TCanvas* c = new TCanvas(Form("c_%s_bar%02d_%s",outName.c_str(),iBar,side),Form("c_%s_bar%02d_%s",outName.c_str(),iBar,side),1200,700);

            gRaw->SetTitle(Form("Bar %02d - Side %s;Threshold [DAC];#sigma(#DeltaT) [ps]",iBar,side));
            gRaw->GetYaxis()->SetRangeUser(40.,100.);

            gRaw->Draw("APL");
            gEneCorr->Draw("PL SAME");
            gEnePhaseCorr->Draw("PL SAME");

            // ---------------------------------------------------
            // --- LEGEND
            // ---------------------------------------------------
            double mean_raw, rms_raw;
            double mean_EneCorr, rms_EneCorr;
            double mean_EnePhaseCorr, rms_EnePhaseCorr;
            computeMeanRMSTGraphErrors(gRaw, mean_raw, rms_raw);
            computeMeanRMSTGraphErrors(gEneCorr, mean_EneCorr, rms_EneCorr);
            computeMeanRMSTGraphErrors(gEnePhaseCorr, mean_EnePhaseCorr, rms_EnePhaseCorr);

            TLegend* leg = new TLegend(0.35,0.66,0.77,0.94);
            leg->SetHeader(Form(" Bar DUT = %02d | Side %s", iBar, side), "C");

            leg->AddEntry(gRaw ,Form("#sigma(#Delta T) Raw (Mean=%.4f, RMS=%.4f)", mean_raw, rms_raw),"pl");
            leg->AddEntry(gEneCorr ,Form("#sigma(#Delta T) TW corr (Mean=%.4f, RMS=%.4f)", mean_EneCorr, rms_EneCorr),"pl");
            leg->AddEntry(gEnePhaseCorr,Form("#sigma(#Delta T) TW+Phase corr (Mean=%.4f, RMS=%.4f)", mean_EnePhaseCorr, rms_EnePhaseCorr),"pl");

            leg->SetBorderSize(1);
            leg->Draw();

            c->SetGridx();
            c->SetGridy();
	    c->SetTickx();
	    c->SetTicky();

            c->Print(Form("%s/%s_bar%02d_%s.png",plotDir.c_str(),outName.c_str(),iBar,side));
            c->Print(Form("%s/%s_bar%02d_%s.pdf",plotDir.c_str(),outName.c_str(),iBar,side));

            delete leg;
            delete c;
        }
    }
}



#endif
