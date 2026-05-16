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



void Plot_FileChi2Value_pol2_pol3()
{
    // -------------------------------
    // --- Hardcoded paths
    // -------------------------------
    std::string file_1 = "pol2/profileComparison_vth7.txt";
    std::string file_2 = "pol3/profileComparison_vth7.txt";

    TGraph* g1 = new TGraph();
    TGraph* g2 = new TGraph();

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

    int ip1 = 0;
    int ip2 = 0;

    std::map<int,double> map1, map2;

    while(fin1 >> bar >> side >> val){
	    if(side == "L") map1[bar] = val;
    }

    while(fin2 >> bar >> side >> val){
            if(side == "L") map2[bar] = val;
    }

    int ip = 0;
    for(const auto& [b, v1] : map1)
    {
    if(map2.find(b) == map2.end()) continue;

    g1->SetPoint(ip, b, v1);
    g2->SetPoint(ip, b, map2[b]);
    ip++;
}


    fin1.close();
    fin2.close();

    // -------------------------------
    // --- Style
    // -------------------------------
    g1->SetMarkerStyle(20);
    g1->SetMarkerColor(kRed);
    g1->SetLineColor(kRed);

    g2->SetMarkerStyle(21);
    g2->SetMarkerColor(kBlue);
    g2->SetLineColor(kBlue);

    // -------------------------------
    // --- Canvas
    // -------------------------------
    TCanvas* c = new TCanvas("c_Chi2_TProfileFit","Chi2 TProfile values",1000,600);

    g1->SetTitle(";Bar Index;#chi^{2}/NDF");
    g1->GetXaxis()->SetTitleSize(0.05);
    g1->GetYaxis()->SetTitleSize(0.05);
    g1->GetYaxis()->SetRangeUser(0.,6.);
    g1->GetXaxis()->SetLimits(-1, 16);
    g1->GetXaxis()->SetNdivisions(17, 0, 0, kTRUE);
    //g1->GetXaxis()->CenterLabels(true);

    g1->Draw("APL");
    g2->Draw("PL SAME");

    // -------------------------------
    // --- Legend
    // -------------------------------
	
    double mean_1, rms_1;    
    double mean_2, rms_2;

    computeMeanRMS(g1, mean_1, rms_1);
    computeMeanRMS(g2, mean_2, rms_2);
    
    TLegend* leg = new TLegend(0.35,0.65,0.85,0.85);

    // Header in grassetto (metodo 1: latex-style ROOT)
    leg->SetHeader(Form("#bf{V_{OV} = 3.00 V | vth2 = 7 DAC | Side L }"), "C");

    // Entries
    leg->AddEntry(g1, Form("Fit function pol2 (Mean=%.4f, RMS=%.4f)", mean_1, rms_1), "lp");
    leg->AddEntry(g2, Form("Fit function pol3 (Mean=%.4f, RMS=%.4f)", mean_2, rms_2), "lp");

    // Stile legenda
    leg->SetFillColor(kWhite);   // sfondo bianco
    leg->SetFillStyle(1001);     // pieno (non trasparente)
    leg->SetBorderSize(1);       // bordo visibile (metti 0 se non lo vuoi)
    leg->SetTextSize(0.03);
    leg->SetTextFont(42);        // font standard

    leg->Draw();
    
    c->SetGridx();
    c->SetGridy();
    c->SetTickx();
    c->SetTicky();
 

    // -------------------------------
    // --- Save
    // -------------------------------
    c->SaveAs("plot_pol2_vs_pol3/Chi2_pol2_vs_pol3_vth7_L.pdf");
}

