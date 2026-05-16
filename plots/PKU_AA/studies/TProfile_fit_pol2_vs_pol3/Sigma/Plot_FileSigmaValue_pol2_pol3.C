#include <iostream>
#include <fstream>
#include <map>

#include "TGraph.h"
#include "TCanvas.h"
#include "TLegend.h"

//------------------------------------------------------------------------
// Compute the RMS and the mean from a TGraph plot
//------------------------------------------------------------------------
void computeMeanRMS(TGraph* g, double& mean, double& rms)
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



void Plot_FileSigmaValue_pol2_pol3()
{
    // -------------------------------
    // --- Hardcoded paths
    // -------------------------------
    std::string file_1 = "pol2/profileComparison_vth32.txt";
    std::string file_2 = "pol3/profileComparison_vth32.txt";

    TGraphErrors* g1_L = new TGraphErrors();
    TGraphErrors* g2_L = new TGraphErrors();
    TGraphErrors* g1_R = new TGraphErrors();
    TGraphErrors* g2_R = new TGraphErrors();


    std::ifstream fin1(file_1);
    std::ifstream fin2(file_2);

    if(!fin1 || !fin2){
        std::cout << "Error opening files!" << std::endl;
        return;
    }

    // -------------------------------
    // --- Read files
    // -------------------------------
    int bar;
    std::string side;
    double val;

    std::map<int,double> map1_L, map2_L, map1_R, map2_R;

    while(fin1 >> bar >> side >> val){
	    if(side == "L") map1_L[bar] = val;
	    if(side == "R") map1_R[bar] = val;
    }

    while(fin2 >> bar >> side >> val){
            if(side == "L") map2_L[bar] = val;
	    if(side == "R") map2_R[bar] = val;
    }

    int ip_L = 0;
    for(const auto& [b, v1] : map1_L)
    {
    if(map2_L.find(b) == map2_L.end()) continue;

    g1_L->SetPoint(ip_L, b, v1);
    g2_L->SetPoint(ip_L, b, map2_L[b]);
    ip_L++;
}

    int ip_R = 0;
    for(const auto& [b, v1] : map1_R)
    {
    if(map2_R.find(b) == map2_R.end()) continue;

    g1_R->SetPoint(ip_R, b, v1);
    g2_R->SetPoint(ip_R, b, map2_R[b]);
    ip_R++;
}

    fin1.close();
    fin2.close();

    // -------------------------------
    // --- Style
    // -------------------------------
    g1_L->SetMarkerStyle(20);
    g1_L->SetMarkerColor(kRed);
    g1_L->SetLineColor(kRed);

    g1_R->SetMarkerStyle(20);
    g1_R->SetMarkerColor(kRed);
    g1_R->SetLineColor(kRed);

    g2_L->SetMarkerStyle(21);
    g2_L->SetMarkerColor(kBlue);
    g2_L->SetLineColor(kBlue);

    g2_R->SetMarkerStyle(21);
    g2_R->SetMarkerColor(kBlue);
    g2_R->SetLineColor(kBlue);

    // -------------------------------
    // --- Canvas
    // -------------------------------
    TCanvas* c_L = new TCanvas("c_Sigma_Values_L","#sigma (#Delta T) with Amplitude walk correction",1000,600);

    g1_L->SetTitle(";Bar Index;#sigma(#Delta T corr)");
    g1_L->GetXaxis()->SetTitleSize(0.05);
    g1_L->GetYaxis()->SetTitleSize(0.05);
    g1_L->GetYaxis()->SetRangeUser(60.,85.);
    g1_L->GetXaxis()->SetLimits(-1, 16);
    g1_L->GetXaxis()->SetNdivisions(17, 0, 0, kTRUE);

    g1_L->Draw("APL");
    g2_L->Draw("PL SAME");

    // -------------------------------
    // --- Legend
    // -------------------------------
	
    double mean_1_L, rms_1_L;    
    double mean_2_L, rms_2_L;

    computeMeanRMS(g1_L, mean_1_L, rms_1_L);
    computeMeanRMS(g2_L, mean_2_L, rms_2_L);
    
    TLegend* leg_L = new TLegend(0.35,0.65,0.85,0.85);

    // Header in grassetto (metodo 1: latex-style ROOT)
    leg_L->SetHeader(Form("#bf{V_{OV} = 3.00 V | vth2 = 32 DAC | Side L }"), "C");

    // Entries
    leg_L->AddEntry(g1_L, Form("Fit function pol2 (Mean=%.4f, RMS=%.4f)", mean_1_L, rms_1_L), "lp");
    leg_L->AddEntry(g2_L, Form("Fit function pol3 (Mean=%.4f, RMS=%.4f)", mean_2_L, rms_2_L), "lp");

    // Stile legenda
    leg_L->SetFillColor(kWhite);   // sfondo bianco
    leg_L->SetFillStyle(1001);     // pieno (non trasparente)
    leg_L->SetBorderSize(1);       // bordo visibile (metti 0 se non lo vuoi)
    leg_L->SetTextSize(0.03);
    leg_L->SetTextFont(42);        // font standard

    leg_L->Draw();
    
    c_L->SetGridx();
    c_L->SetGridy();
    c_L->SetTickx();
    c_L->SetTicky();
 

    // -------------------------------
    // --- Save
    // -------------------------------
    c_L->SaveAs("plot_pol2_vs_pol3/Sigma_DeltaT_Corr_pol2_vs_pol3_vth32_L.pdf");

 
    // -------------------------------
    // --- Canvas
    // -------------------------------
    TCanvas* c_R = new TCanvas("c_Sigma_Values_R","#sigma (#Delta T) with Amplitude walk correction",1000,600);
 
    g1_R->SetTitle(";Bar Index;#sigma(#Delta T corr)");
    g1_R->GetXaxis()->SetTitleSize(0.05);
    g1_R->GetYaxis()->SetTitleSize(0.05);
    g1_R->GetYaxis()->SetRangeUser(60.,85.);
    g1_R->GetXaxis()->SetLimits(-1, 16);
    g1_R->GetXaxis()->SetNdivisions(17, 0, 0, kTRUE);

    g1_R->Draw("APL");
    g2_R->Draw("PL SAME");

    // -------------------------------
    // --- Legend
    // -------------------------------
   
    double mean_1_R, rms_1_R;
    double mean_2_R, rms_2_R;

    computeMeanRMS(g1_R, mean_1_R, rms_1_R);
    computeMeanRMS(g2_R, mean_2_R, rms_2_R);

    TLegend* leg_R = new TLegend(0.35,0.65,0.85,0.85);

    // Header in grassetto (metodo 1: latex-style ROOT)
    leg_R->SetHeader(Form("#bf{V_{OV} = 3.00 V | vth2 = 32 DAC | Side R }"), "C");

    // Entries
    leg_R->AddEntry(g1_R, Form("Fit function pol2 (Mean=%.4f, RMS=%.4f)", mean_1_R, rms_1_R), "lp");
    leg_R->AddEntry(g2_R, Form("Fit function pol3 (Mean=%.4f, RMS=%.4f)", mean_2_R, rms_2_R), "lp");

    // Stile legenda
    leg_R->SetFillColor(kWhite);   // sfondo bianco
    leg_R->SetFillStyle(1001);     // pieno (non trasparente)
    leg_R->SetBorderSize(1);       // bordo visibile (metti 0 se non lo vuoi)
    leg_R->SetTextSize(0.03);
    leg_R->SetTextFont(42);        // font standard

    leg_R->Draw();

    c_R->SetGridx();
    c_R->SetGridy();
    c_R->SetTickx();
    c_R->SetTicky();
    

    // -------------------------------
    // --- Save
    // -------------------------------
    c_R->SaveAs("plot_pol2_vs_pol3/Sigma_DeltaT_Corr_pol2_vs_pol3_vth32_R.pdf");
    
    }

