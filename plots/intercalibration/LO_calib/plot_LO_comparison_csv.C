#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>

#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <vector>

// struttura per contenere dati
struct LOData {
    std::vector<double> bar_L, calib_L;
    std::vector<double> bar_R, calib_R;
};

// -----------------------------------------------------
// funzione per leggere CSV
// -----------------------------------------------------
LOData readCSV(const std::string& filename) {
    LOData data;
    std::ifstream file(filename);

    std::string line;
    std::getline(file, line); // skip header

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string item;

        int bar;
        std::string side;
        double calib;

        std::getline(ss, item, ',');
        bar = std::stoi(item);

        std::getline(ss, side, ',');

        std::getline(ss, item, ',');
        calib = std::stod(item);

        if (side == "L") {
            data.bar_L.push_back(bar);
            data.calib_L.push_back(calib);
        } else if (side == "R") {
            data.bar_R.push_back(bar);
            data.calib_R.push_back(calib);
        }
    }

    return data;
}

// -----------------------------------------------------
// MAIN MACRO
// -----------------------------------------------------
void plot_LO_comparison() {

    gStyle->SetOptStat(0);

    // -------- file paths --------
    std::string file_global = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/module_32110020004436_LO_calibration_factors_global.csv";
    std::string file_sep    = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/module_32110020004436_LO_calibration_factors_separated_side.csv";

    std::string outdir = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/intercalibration/LO_calib/";

    // -------- read data --------
    LOData global = readCSV(file_global);
    LOData sep    = readCSV(file_sep);

    // =====================================================
    // 🔹 PLOT 1: GLOBAL
    // =====================================================
    TCanvas* c1 = new TCanvas("c1","Global",800,600);

    TGraph* gL = new TGraph(global.bar_L.size(), &global.bar_L[0], &global.calib_L[0]);
    TGraph* gR = new TGraph(global.bar_R.size(), &global.bar_R[0], &global.calib_R[0]);

    gL->SetMarkerStyle(20);
    gL->SetMarkerColor(kRed);
    gL->SetLineColor(kRed);

    gR->SetMarkerStyle(21);
    gR->SetMarkerColor(kBlue);
    gR->SetLineColor(kBlue);

    gL->SetTitle("LO calibration (global);Bar;Calibration factor");
    gL->Draw("AP");
    gR->Draw("P SAME");

    c1->SaveAs((outdir + "LO_global.pdf").c_str());

    // =====================================================
    // 🔹 PLOT 2: SEPARATED SIDE
    // =====================================================
    TCanvas* c2 = new TCanvas("c2","Separated",800,600);

    TGraph* gL2 = new TGraph(sep.bar_L.size(), &sep.bar_L[0], &sep.calib_L[0]);
    TGraph* gR2 = new TGraph(sep.bar_R.size(), &sep.bar_R[0], &sep.calib_R[0]);

    gL2->SetMarkerStyle(20);
    gL2->SetMarkerColor(kRed);
    gL2->SetLineColor(kRed);

    gR2->SetMarkerStyle(21);
    gR2->SetMarkerColor(kBlue);
    gR2->SetLineColor(kBlue);

    gL2->SetTitle("LO calibration (separated side);Bar;Calibration factor");
    gL2->Draw("AP");
    gR2->Draw("P SAME");

    c2->SaveAs((outdir + "LO_separated.pdf").c_str());

    // =====================================================
    // 🔹 PLOT 3: DIFFERENCE (global - separated)
    // =====================================================
    TCanvas* c3 = new TCanvas("c3","Difference",800,600);

    std::vector<double> diff_L, diff_R;

    for (size_t i = 0; i < global.calib_L.size(); ++i) {
        diff_L.push_back(global.calib_L[i] - sep.calib_L[i]);
        diff_R.push_back(global.calib_R[i] - sep.calib_R[i]);
    }

    TGraph* gDiffL = new TGraph(global.bar_L.size(), &global.bar_L[0], &diff_L[0]);
    TGraph* gDiffR = new TGraph(global.bar_R.size(), &global.bar_R[0], &diff_R[0]);

    gDiffL->SetMarkerStyle(20);
    gDiffL->SetMarkerColor(kRed);
    gDiffL->SetLineColor(kRed);

    gDiffR->SetMarkerStyle(21);
    gDiffR->SetMarkerColor(kBlue);
    gDiffR->SetLineColor(kBlue);

    gDiffL->SetTitle("Difference (global - separated);Bar;#Delta calibration");
    gDiffL->Draw("AP");
    gDiffR->Draw("P SAME");

    c3->SaveAs((outdir + "LO_difference.pdf").c_str());

    std::cout << "Plots saved in: " << outdir << std::endl;
}
