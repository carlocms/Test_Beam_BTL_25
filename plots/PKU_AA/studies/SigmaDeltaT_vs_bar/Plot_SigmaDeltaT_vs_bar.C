#include <iostream>
#include <fstream>
#include <map>
#include <cmath>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

//------------------------------------------------------------------------
// Compute mean and RMS from TGraph
//------------------------------------------------------------------------
void computeMeanRMS(TGraph* g, double& mean, double& rms)
{
    int n = g->GetN();
    if(n == 0){
        mean = 0.0;
        rms  = 0.0;
        return;
    }

    double sum = 0.0, sum2 = 0.0;

    for(int i = 0; i < n; ++i){
        double x, y;
        g->GetPoint(i, x, y);
        sum  += y;
        sum2 += y*y;
    }

    mean = sum / n;
    rms  = std::sqrt(sum2/n - mean*mean);
}

//------------------------------------------------------------------------
// Main macro
//------------------------------------------------------------------------
void Plot_SigmaDeltaT_vs_bar(int vth)
{
    // -------------------------------
    // --- Hardcoded base path
    // -------------------------------
    std::string baseDir = "SigmaDeltaT_raw/";
    std::string file = baseDir + Form("profileComparison_vth%d.txt", vth);

    std::ifstream fin(file);

    if(!fin){
        std::cout << "Error opening file: " << file << std::endl;
        return;
    }

    // -------------------------------
    // --- Maps
    // -------------------------------
    std::map<int,double> map_L;
    std::map<int,double> map_R;

    int bar;
    std::string side;
    double val, valErr;

    // -------------------------------
    // --- Read file
    // -------------------------------
    while(fin >> bar >> side >> val >> valErr)
    {
        if(side == "L") map_L[bar] = valErr;
        else if(side == "R") map_R[bar] = valErr;
    }

    fin.close();

    // -------------------------------
    // --- Graphs
    // -------------------------------
    TGraph* gL = new TGraph();
    TGraph* gR = new TGraph();

    int ipL = 0;
    for(const auto& [b, v] : map_L){
        gL->SetPoint(ipL, b, v);
        ipL++;
    }

    int ipR = 0;
    for(const auto& [b, v] : map_R){
        gR->SetPoint(ipR, b, v);
        ipR++;
    }

    // -------------------------------
    // --- Style
    // -------------------------------
    gL->SetMarkerStyle(20);
    gL->SetMarkerColor(kRed);
    gL->SetLineColor(kRed);

    gR->SetMarkerStyle(21);
    gR->SetMarkerColor(kBlue);
    gR->SetLineColor(kBlue);

    // -------------------------------
    // --- Canvas
    // -------------------------------
    TCanvas* c = new TCanvas("c_sigma","Sigma DeltaT vs bar",1000,600);

    gL->SetTitle(";Bar Index;#sigma(#Delta T_{raw}) [ps]");
    gL->GetXaxis()->SetLimits(-1,16);
    gL->GetYaxis()->SetRangeUser(70.,100.);
    gL->GetXaxis()->SetNdivisions(17,0,0,kTRUE);

    gL->GetXaxis()->SetTitleSize(0.05);
    gL->GetYaxis()->SetTitleSize(0.05);
    gL->GetYaxis()->SetTitleOffset(0.95);

    gL->Draw("APL");
    gR->Draw("PL SAME");

    // -------------------------------
    // --- Stats
    // -------------------------------
    double meanL, rmsL;
    double meanR, rmsR;

    computeMeanRMS(gL, meanL, rmsL);
    computeMeanRMS(gR, meanR, rmsR);

    // -------------------------------
    // --- Legend
    // -------------------------------
    TLegend* leg = new TLegend(0.35,0.65,0.7,0.85);

    leg->SetHeader(Form("#bf{V_{OV} = 3.00 V | vth2 = %d DAC}", vth), "C");

    leg->AddEntry(gL, Form("Side L (Mean=%.2f, RMS=%.2f)", meanL, rmsL), "lp");
    leg->AddEntry(gR, Form("Side R (Mean=%.2f, RMS=%.2f)", meanR, rmsR), "lp");

    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.03);
    leg->SetTextFont(42);

    leg->Draw();

    // -------------------------------
    // --- Grid
    // -------------------------------
    c->SetGridx();
    c->SetGridy();
    c->SetTickx();
    c->SetTicky();

    // -------------------------------
    // --- Save
    // -------------------------------
    std::string outName = Form("plot_Sigma_DeltaT_Raw/Sigma_DeltaT_Raw_vs_bar_vth%d", vth);

    c->SaveAs((outName + ".pdf").c_str());
//    c->SaveAs((outName + ".png").c_str());
}
