#include <iostream>
#include <fstream>
#include <map>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"

std::map<int,double> readFile(const std::string& filename)
{
    std::map<int,double> data;
    std::ifstream fin(filename);

    int bar;
    double val, valErr;

    while(fin >> bar >> val >> valErr)
        data[bar] = val;

    fin.close();
    return data;
}

void Plot_DeltaDelta(
    std::string file_L1,
    std::string file_R1,
    std::string file_L2,
    std::string file_R2
)
{
    // -------------------------------
    // --- Read files
    // -------------------------------
    auto L1 = readFile(file_L1);
    auto L2 = readFile(file_L2);
    auto R1 = readFile(file_R1);
    auto R2 = readFile(file_R2);

    // -------------------------------
    // --- Check bar 0 exists
    // -------------------------------
    if (!L1.count(0) || !L2.count(0) || !R1.count(0) || !R2.count(0)) {
        std::cerr << "ERROR: bar 0 missing in one of the files!" << std::endl;
        return;
    }

    double refL1 = L1[0];
    double refL2 = L2[0];
    double refR1 = R1[0];
    double refR2 = R2[0];

    // -------------------------------
    // --- Graphs
    // -------------------------------
    TGraph* gL = new TGraph();
    TGraph* gR = new TGraph();

    int ipL = 0;
    for (const auto& [b, v1] : L1) {
        if (!L2.count(b)) continue;

        double d1 = v1 - refL1;
        double d2 = L2[b] - refL2;

        double deltaDelta = d1 - d2;

        gL->SetPoint(ipL, b, deltaDelta);
        ipL++;
    }

    int ipR = 0;
    for (const auto& [b, v1] : R1) {
        if (!R2.count(b)) continue;

        double d1 = v1 - refR1;
        double d2 = R2[b] - refR2;

        double deltaDelta = d1 - d2;

        gR->SetPoint(ipR, b, deltaDelta);
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
    TCanvas* c = new TCanvas("c_deltadelta","DeltaDelta comparison",800,600);

    gL->SetTitle(";Bar Index;#Delta#Delta mean [ps]");

    gL->Draw("AP");
    gR->Draw("P SAME");

    // -------------------------------
    // --- Legend
    // -------------------------------
    TLegend* leg = new TLegend(0.20,0.75,0.45,0.88);

    leg->AddEntry(gL,"Side L (#Delta#Delta)","p");
    leg->AddEntry(gR,"Side R (#Delta#Delta)","p");

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    c->SetGrid();

    // -------------------------------
    // --- Save
    // -------------------------------
    c->SaveAs("DeltaDelta_LR.pdf");
    c->SaveAs("DeltaDelta_LR.png");
}
