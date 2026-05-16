//============================================================
// Plot_DeltaTdiff.C
//============================================================

#include <iostream>
#include <fstream>
#include <map>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"

void Plot_FileValue(std::string file_L, std::string file_R)
{
    TGraph* gL = new TGraph();
    TGraph* gR = new TGraph();

    // -------------------------------
    // --- Read txt files
    // -------------------------------
    std::map<int,double> Side_L;
    std::map<int,double> Side_R;

    std::ifstream finL(file_L);
    std::ifstream finR(file_R);

    int bar;
    double val, valErr;

    while(finL >> bar >> val >> valErr)
        Side_L[bar] = val;

    while(finR >> bar >> val >> valErr)
        Side_R[bar] = val;

    finL.close();
    finR.close();

    // -------------------------------
    // --- Reference values (bar 0)
    // -------------------------------
    if (Side_L.find(0) == Side_L.end() || Side_R.find(0) == Side_R.end()) {
        std::cerr << "ERROR: bar 0 not found in one of the files!" << std::endl;
        return;
    }

    double refL = Side_L[0];
    double refR = Side_R[0];

    // -------------------------------
    // --- Fill graphs
    // -------------------------------
    int ipL = 0;
    for (const auto& [b, v] : Side_L) {
        double delta = v - refL;
        gL->SetPoint(ipL, b, delta);
        ipL++;
    }

    int ipR = 0;
    for (const auto& [b, v] : Side_R) {
        double delta = v - refR;
        gR->SetPoint(ipR, b, delta);
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
    TCanvas* c = new TCanvas("c_relative_diff","Relative DeltaT diff",1000,600);

    gL->SetTitle(";Bar Index;#Delta mean (w.r.t bar 0) [ps]");

    gL->Draw("AP");
    gR->Draw("P SAME");

    // -------------------------------
    // --- Legend
    // -------------------------------
    TLegend* leg = new TLegend(0.20,0.75,0.45,0.88);

    leg->AddEntry(gL,"Side L","p");
    leg->AddEntry(gR,"Side R","p");

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    c->SetGrid();

    // -------------------------------
    // --- Save
    // -------------------------------
    c->SaveAs("LO_vs_DeltaTdiff.pdf");
    c->SaveAs("LO_vs_DeltaTdiff.png");
}
