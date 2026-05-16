//============================================================
// Plot_LO_vs_DeltaTdiff_fromTwoFiles.C
//============================================================

#include <iostream>
#include <fstream>
#include <map>
#include <tuple>
#include <sstream>

void Plot_LO_vs_DeltaTdiff_fromTwoFiles(
    std::string fitFile,
    std::string profileFile,
    std::string LO_L_file,
    std::string LO_R_file
)
{
    // --------------------------------------------------
    // --- MAPPE
    // --------------------------------------------------
    std::map<std::pair<int,std::string>, double> fitMap;
    std::map<std::pair<int,std::string>, double> profMap;

    std::map<int,double> LO_L;
    std::map<int,double> LO_R;

    // --------------------------------------------------
    // --- READ FIT FILE
    // --------------------------------------------------
    std::ifstream finFit(fitFile);
    int bar;
    std::string side;
    double val;

    while(finFit >> bar >> side >> val)
    {
        fitMap[{bar, side}] = val;
    }
    finFit.close();

    // --------------------------------------------------
    // --- READ PROFILE FILE
    // --------------------------------------------------
    std::ifstream finProf(profileFile);

    while(finProf >> bar >> side >> val)
    {
        profMap[{bar, side}] = val;
    }
    finProf.close();

    // --------------------------------------------------
    // --- READ LO FILES
    // --------------------------------------------------
    std::ifstream finL(LO_L_file);
    std::ifstream finR(LO_R_file);

    while(finL >> bar >> val) LO_L[bar] = val;
    while(finR >> bar >> val) LO_R[bar] = val;

    finL.close();
    finR.close();

    // --------------------------------------------------
    // --- GRAPH
    // --------------------------------------------------
    TGraph* gL = new TGraph();
    TGraph* gR = new TGraph();

    int pL = 0;
    int pR = 0;

    // --------------------------------------------------
    // --- BUILD DIFFERENCE + GRAPH
    // --------------------------------------------------
    for(const auto& it : fitMap)
    {
        int bar = it.first.first;
        std::string side = it.first.second;

        // controlla esistenza profile
        if(profMap.find(it.first) == profMap.end()) continue;

        double mean_fit  = it.second;
        double mean_prof = profMap[it.first];

        double diff = mean_fit - mean_prof;

        if(side == "L")
        {
            if(LO_L.find(bar) == LO_L.end()) continue;

            double LO = LO_L[bar];
            gL->SetPoint(pL++, LO, diff);
        }
        else if(side == "R")
        {
            if(LO_R.find(bar) == LO_R.end()) continue;

            double LO = LO_R[bar];
            gR->SetPoint(pR++, LO, diff);
        }
    }

    // --------------------------------------------------
    // --- STYLE
    // --------------------------------------------------
    gL->SetMarkerStyle(20);
    gL->SetMarkerColor(kRed);
    gL->SetLineColor(kRed);

    gR->SetMarkerStyle(21);
    gR->SetMarkerColor(kBlue);
    gR->SetLineColor(kBlue);

    // --------------------------------------------------
    // --- CANVAS
    // --------------------------------------------------
    TCanvas* c = new TCanvas("c_LO_vs_diff","LO vs DeltaT diff",800,600);

    gL->SetTitle(";LO [a.u.];Mean_{fit} - Mean_{profile} [ps]");

    gL->Draw("AP");
    gR->Draw("P SAME");

    // --------------------------------------------------
    // --- LEGEND
    // --------------------------------------------------
    TLegend* leg = new TLegend(0.20,0.75,0.45,0.88);

    leg->AddEntry(gL,"Side L","p");
    leg->AddEntry(gR,"Side R","p");

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    c->SetGrid();

    // --------------------------------------------------
    // --- OPTIONAL: FIT LINEARE
    // --------------------------------------------------
    gL->Fit("pol1","Q");
    gR->Fit("pol1","Q+");

    // --------------------------------------------------
    // --- SAVE
    // --------------------------------------------------
    c->Print("LO_vs_DeltaTdiff_fromTwoFiles.pdf");
    c->Print("LO_vs_DeltaTdiff_fromTwoFiles.png");
}
