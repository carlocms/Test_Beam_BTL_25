#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TStyle.h>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>

// struttura dati
struct LOData {
    std::vector<double> bar_L, calib_L;
    std::vector<double> bar_R, calib_R;
};

// -----------------------------------------------------
// lettura TXT (space-separated)
// -----------------------------------------------------
LOData readTXT(const std::string& filename) {
    LOData data;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "ERROR: cannot open file " << filename << std::endl;
        return data;
    }

    std::string line;

    while (std::getline(file, line)) {

        std::stringstream ss(line);

        int bar;
        std::string side;
        double calib;

        ss >> bar >> side >> calib;

        if (ss.fail()) {
            std::cout << "Skipping line: " << line << std::endl;
            continue;
        }

        if (side == "L") {
            data.bar_L.push_back(bar);
            data.calib_L.push_back(calib);
        } else if (side == "R") {
            data.bar_R.push_back(bar);
            data.calib_R.push_back(calib);
        }
    }

    std::cout << "Read file: " << filename << std::endl;
    std::cout << "  L points: " << data.bar_L.size() << std::endl;
    std::cout << "  R points: " << data.bar_R.size() << std::endl;

    return data;
}

void computeStats(const std::vector<double>& v, double& mean, double& rms) {
    mean = 0;
    rms  = 0;

    if (v.empty()) return;

    for (auto val : v) mean += val;
    mean /= v.size();

    for (auto val : v) rms += (val - mean)*(val - mean);
    rms = std::sqrt(rms / v.size());
}


// -----------------------------------------------------
// MAIN
// -----------------------------------------------------
void plot_LO_comparison_txt() {

    gStyle->SetOptStat(0);

    // -------- file paths --------
    std::string file_LO_L = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/Mean_DeltaT_vs_LO/LO_L.txt";
    std::string file_Mean_L    = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/Mean_DeltaT_vs_LO/Mean_DeltaT_L_vth12.txt";
    std::string file_LO_R = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/Mean_DeltaT_vs_LO/LO_R.txt";
    std::string file_Mean_R    = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/Mean_DeltaT_vs_LO/Mean_DeltaT_R_vth12.txt";
    std::string outdir = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/Mean_DeltaT_vs_LO/";

    // -------- read data --------
    LOData LO_L = readTXT(file_LO_L);
    LOData Mean_L    = readTXT(file_Mean_L);
    LOData LO_R = readTXT(file_LO_R);
    LOData Mean_R    = readTXT(file_Mean_R);

    // =====================================================
    // 🔹 PLOT 1: GLOBAL
    // =====================================================
    TCanvas* c1 = new TCanvas("c1","Global",800,600);

    TGraph* gL = new TGraph(LO_L.Value_LO_L.size(), &LO_L.Value_LO_L[0], &Mean_L.Values_Mean_L[0]);
    TGraph* gR = new TGraph(LO_R.Value_LO_R.size(), &LO_R.Value_LO_R[0], &Mean_R.Values_Mean_R[0]);

    gL->SetMarkerStyle(20);
    gL->SetMarkerColor(kRed);
    gL->SetLineColor(kRed);

    gR->SetMarkerStyle(21);
    gR->SetMarkerColor(kBlue);
    gR->SetLineColor(kBlue);

    // grid
    c1->SetGridy();
    c1->SetTickx();
    c1->SetTicky();

    gL->SetTitle("LO Values (p.e. /MeV) vs Mean #Delta T; LO Values (p.e. /MeV); Mean #Delta T [ps]");
    gL->GetYaxis()->SetRangeUser(0.93, 1.09);
    gL->Draw("APL");
    gR->Draw("PL SAME");
 
//    double mean_L, rms_L, mean_R, rms_R;
//computeStats(global.calib_L, mean_L, rms_L);
//computeStats(global.calib_R, mean_R, rms_R);

    // legenda
    TLegend* leg1 = new TLegend(0.27,0.70,0.65,0.9);
    leg1->SetHeader(Form("LO calibration - global"), "C");
    leg1->AddEntry(gL, Form("Side L","lp");
    leg1->AddEntry(gR, Form("Side R","lp");
    leg1->SetFillStyle(0);
    leg1->SetBorderSize(0);
    leg1->SetTextSize(0.03);
    leg1->SetTextFont(42);
    leg1->Draw();




    c1->SaveAs((outdir + "LO_factors_global.pdf").c_str());

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
    gL2->GetYaxis()->SetRangeUser(0.93, 1.09); 
    gL2->Draw("APL");
    gR2->Draw("PL SAME");

    double mean_L2, rms_L2, mean_R2, rms_R2;
computeStats(sep.calib_L, mean_L2, rms_L2);
computeStats(sep.calib_R, mean_R2, rms_R2);

TLegend* leg2 = new TLegend(0.27,0.70,0.65,0.9);
leg2->SetHeader("LO calibration - separated", "C");

leg2->AddEntry(gL2, Form("Side L (Mean=%.4f, RMS=%.4f)", mean_L2, rms_L2), "lp");
leg2->AddEntry(gR2, Form("Side R (Mean=%.4f, RMS=%.4f)", mean_R2, rms_R2), "lp");

leg2->SetFillStyle(0);
leg2->SetBorderSize(0);
leg2->SetTextSize(0.03);
leg2->Draw();

    c2->SaveAs((outdir + "LO_separated.pdf").c_str());

    // =====================================================
    // 🔹 PLOT 3: DIFFERENCE
    // =====================================================
    TCanvas* c3 = new TCanvas("c3","Difference",800,600);

    std::vector<double> diff_L, diff_R;

    for (size_t i = 0; i < global.calib_L.size(); ++i) {
        diff_L.push_back(global.calib_L[i] - sep.calib_L[i]);
        diff_R.push_back(global.calib_R[i] - sep.calib_R[i]);

	std::cout << "Bar index: " << i+1 << " diff_L: " << global.calib_L[i] - sep.calib_L[i] << " diff_R: " << global.calib_R[i] - sep.calib_R[i] << std::endl; 
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
    gDiffL->GetYaxis()->SetRangeUser(-0.0015, 0.0015);
    gDiffL->Draw("APL");
    gDiffR->Draw("PL SAME");


    double mean_dL, rms_dL, mean_dR, rms_dR;
computeStats(diff_L, mean_dL, rms_dL);
computeStats(diff_R, mean_dR, rms_dR);

TLegend* leg3 = new TLegend(0.27,0.70,0.65,0.9);
leg3->SetHeader("Difference (global - separated)", "C");

leg3->AddEntry(gDiffL, Form("Side L (Mean=%.5f, RMS=%.5f)", mean_dL, rms_dL), "lp");
leg3->AddEntry(gDiffR, Form("Side R (Mean=%.5f, RMS=%.5f)", mean_dR, rms_dR), "lp");

leg3->SetFillStyle(0);
leg3->SetBorderSize(0);
leg3->SetTextSize(0.03);
leg3->Draw();

    c3->SaveAs((outdir + "LO_difference.pdf").c_str());

    std::cout << "Plots saved in: " << outdir << std::endl;
}
