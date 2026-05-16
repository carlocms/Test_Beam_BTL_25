#include <iostream>
#include <fstream>
#include <map>
#include <iomanip>  // per padding

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TF1.h"

void Plot_slope_vs_calibFactors(int vth)
{
    // =========================================================
    // --- PATH
    // =========================================================
    std::string basePath = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/";  // <-- metti qui la tua directory

    std::string fileSlope = basePath + Form("/TProfile_fit_slope/E_l350_E_h500/profileComparison_vth%d.txt", vth);

    std::cout << "fileSlope: " << fileSlope << std::endl;
    std::string fileCalib = basePath + Form("/intercalibration/LO_calib/module_32110020004436_LO_calibration_factors_global.txt");

    std::ifstream finSlope(fileSlope);
    std::ifstream finCalib(fileCalib);


    if(!finSlope || !finCalib){
        std::cout << "Error opening input files!" << std::endl;
        return;
    }

    // =========================================================
    // --- MAPPE
    // =========================================================
    std::map<int,double> slope_L, slope_R;
    std::map<int,double> calib_L, calib_R;

    int bar;
    std::string side;
    double val;

    // --- Read slope
    while(finSlope >> bar >> side >> val)
    {
        if(side == "L") slope_L[bar] = val;
        if(side == "R") slope_R[bar] = val;
    }

    // --- Read calibration
    while(finCalib >> bar >> side >> val)
    {
        if(side == "L") calib_L[bar] = val;
        if(side == "R") calib_R[bar] = val;
    }

    finSlope.close();
    finCalib.close();

    // =========================================================
    // --- GRAPH L
    // =========================================================
    TGraph* gL = new TGraph();
    int ipL = 0;

    for(const auto& [b, slope] : slope_L)
    {
        if(calib_L.find(b) == calib_L.end()) continue;

        double calib = calib_L[b];

        gL->SetPoint(ipL, calib, slope);
        ipL++;
    }

    // =========================================================
    // --- GRAPH R
    // =========================================================
    TGraph* gR = new TGraph();
    int ipR = 0;

    for(const auto& [b, slope] : slope_R)
    {
        if(calib_R.find(b) == calib_R.end()) continue;

        double calib = calib_R[b];

        gR->SetPoint(ipR, calib, slope);
        ipR++;
    }

    // =========================================================
    // --- STYLE
    // =========================================================
    gL->SetMarkerStyle(20);
    gL->SetMarkerColor(kRed);
    gL->SetLineColor(kRed);

    gR->SetMarkerStyle(21);
    gR->SetMarkerColor(kBlue);
    gR->SetLineColor(kBlue);

    // =========================================================
    // --- CANVAS L
    // =========================================================
    TCanvas* cL = new TCanvas("cL","Slope vs LO Calibration L",850,600);

    gL->SetTitle(";LO Calibration factor;Slope [ps/a.u.]");
    gL->GetXaxis()->SetTitleSize(0.05);
    gL->GetYaxis()->SetTitleSize(0.05);
    gL->GetYaxis()->SetTitleOffset(0.8);
    gL->GetYaxis()->SetRangeUser(-2.,1.);

    gL->Draw("AP");

    // =========================================================
    // --- FIT LINEARE
    // =========================================================
    TF1* fitL = new TF1("fitL","pol1",0.7,1.4);

    fitL->SetLineColor(kBlack);
    fitL->SetLineWidth(2);

    gL->Fit(fitL,"RQ");



    TLegend* legL = new TLegend(0.20,0.75,0.60,0.88);
    legL->SetHeader(Form("#bf{Side L | vth = %d DAC}", vth),"C");
    //legL->AddEntry(gL,"Slope vs TOFHIR calibration factor","pl");
    legL->AddEntry(gL,"Data","pl");

    legL->AddEntry(fitL,Form("Linear fit: y = %.3f x %+ .3f",fitL->GetParameter(1),fitL->GetParameter(0)),"l");

    legL->SetFillColor(kWhite);
    legL->SetFillStyle(1001);
    legL->SetBorderSize(1);
    legL->SetTextSize(0.03);
    legL->Draw();

    cL->SetGrid();
    cL->SaveAs(Form("plot_slope_vs_LO_factors/Slope_vs_LO_Calibration_L_vth%d.pdf", vth));

    // =========================================================
    // --- CANVAS R
    // =========================================================
    TCanvas* cR = new TCanvas("cR","Slope vs LO Calibration R",850,600);

    gR->SetTitle(";LO Calibration factor;Slope [ps/a.u.]");
    gR->GetXaxis()->SetTitleSize(0.05);
    gR->GetYaxis()->SetTitleSize(0.05);
    gR->GetYaxis()->SetTitleOffset(0.8);
    gR->GetYaxis()->SetRangeUser(-2.,1.);

    gR->Draw("AP");
 
    // =========================================================
    // --- FIT LINEARE
    // =========================================================
    TF1* fitR = new TF1("fitR","pol1",0.7,1.4);

    fitR->SetLineColor(kBlack);
    fitR->SetLineWidth(2);

    gR->Fit(fitR,"RQ");


    TLegend* legR = new TLegend(0.20,0.75,0.60,0.88);
    legR->SetHeader(Form("#bf{Side R | vth = %d DAC}", vth),"C");
    //legR->AddEntry(gR,"Slope vs LO calibration factor","pl");

    legR->AddEntry(gR,"Data","pl");

    legR->AddEntry(fitR,Form("Linear fit: y = %.3f x %+ .3f",fitR->GetParameter(1),fitR->GetParameter(0)),"l");

//legR->AddEntry(
//    (TObject*)0,
//    Form("#chi^{2}/NDF = %.2f / %d",
//         fitR->GetChisquare(),
//         fitR->GetNDF()),
//    ""
//);

    legR->SetFillColor(kWhite);
    legR->SetFillStyle(1001);
    legR->SetBorderSize(1);
    legR->SetTextSize(0.03);
    legR->Draw();

    cR->SetGrid();
    cR->SaveAs(Form("plot_slope_vs_LO_factors/Slope_vs_LO_Calibration_R_vth%d.pdf", vth));

    // cleanup (opzionale)
//    delete legL;
//    delete legR;
//    delete cL;
//    delete cR;
//    delete fitL;
//    delete fitR;
}
