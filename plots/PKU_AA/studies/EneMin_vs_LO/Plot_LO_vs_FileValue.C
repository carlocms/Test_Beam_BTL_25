//============================================================
// Plot_LO_vs_DeltaTdiff.C
//============================================================

#include <iostream>
#include <fstream>
#include <map>

void Plot_LO_vs_FileValue(
    std::string FileValue,
    std::string LO_L_file,
    std::string LO_R_file
)
{
    // -------------------------------
    // --- Read LO files
    // -------------------------------
    std::map<int,double> LO_L;
    std::map<int,double> LO_R;

    std::ifstream finL(LO_L_file);
    std::ifstream finR(LO_R_file);

    int bar;
    double val;

    while(finL >> bar >> val)
        LO_L[bar] = val;

    while(finR >> bar >> val)
        LO_R[bar] = val;

    finL.close();
    finR.close();

    // -------------------------------
    // --- Graphs
    // -------------------------------
    TGraph* gL = new TGraph();
    TGraph* gR = new TGraph();

    int pL = 0;
    int pR = 0;

    // -------------------------------
    // --- Read diff file
    // -------------------------------
    std::ifstream fin(FileValue);

    int barFile;
    std::string side;
    double Value;

    std::string line;

    while(std::getline(fin, line))
    {
        if(line[0] == '#') continue;

        std::stringstream ss(line);
        ss >> barFile >> side >> Value;

        if(side == "L")
        {
            if(LO_L.find(barFile) == LO_L.end()) continue;

            double LO = LO_L[barFile];
            gL->SetPoint(pL++, LO, Value);
        }
        else if(side == "R")
        {
            if(LO_R.find(barFile) == LO_R.end()) continue;

            double LO = LO_R[barFile];
            gR->SetPoint(pR++, LO, Value);
        }
    }

    fin.close();

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
    TCanvas* c = new TCanvas("c_LO_vs_diff","LO vs DeltaT diff",800,600);

    gL->SetTitle(";LO [a.u.];#Delta(mean) [ps]");

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
    c->Print("LO_vs_DeltaTdiff.pdf");
    c->Print("LO_vs_DeltaTdiff.png");
}
