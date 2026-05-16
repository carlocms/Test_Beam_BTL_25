#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TSystem.h"

void Plot_MeanDeltaT_vs_REF_forDifferentVth(int dutBar = 7)
{
    std::string basePath = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/DifferentREF_studies";

    // REF bars disponibili
    std::vector<int> refBars = {4,5,6,7,8,9,10,11};

    // Threshold da confrontare
    std::vector<int> vthList = {7,8,9,10,11,12,13,18,25,32};

    // colori
    std::vector<int> colors = {
        kRed, kBlue, kGreen+2, kMagenta, kCyan+2,
        kOrange+7, kBlack, kViolet, kTeal+2, kPink+7, kAzure+2
    };

    // canvas
    TCanvas* c = new TCanvas("c_vthScan","DeltaT vs REF (vth scan)",1200,700);

    TLegend* leg = new TLegend(0.55,0.65,0.88,0.88);
    leg->SetNColumns(2);
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.03);
    leg->SetHeader(Form("#bf{DUT bar = %d}", dutBar), "C");

    bool first = true;

    // ============================================
    // LOOP sui vth
    // ============================================
    for (size_t ivth = 0; ivth < vthList.size(); ++ivth)
    {
        int vth = vthList[ivth];

        TGraphErrors* g = new TGraphErrors();
        int ip = 0;

        // ----------------------------------------
        // LOOP sui REF bars
        // ----------------------------------------
        for (int ref : refBars)
        {
            std::string filePath = Form("%s/REF_Bar%d/outfileFile_vth%d.txt",
                                        basePath.c_str(), ref, vth);

            std::ifstream fin(filePath);
            if (!fin.is_open()) {
                std::cout << "Missing file: " << filePath << std::endl;
                continue;
            }

            int bar;
            double mean, sigma;

            std::string line;
            std::getline(fin, line); // skip header

            while (fin >> bar >> mean >> sigma)
            {
                if (bar != dutBar) continue;

                g->SetPoint(ip, ref, mean);
                g->SetPointError(ip, 0., sigma);
                ip++;
            }

            fin.close();
        }

        if (g->GetN() == 0) continue;

        // stile
        g->SetMarkerStyle(20 + ivth % 5);
        g->SetMarkerColor(colors[ivth % colors.size()]);
        g->SetLineColor(colors[ivth % colors.size()]);

        if (first)
        {
            g->SetTitle(";REF bar;Mean #Delta T (L-R) [ps]");
            g->GetXaxis()->SetLimits(3, 12);
            g->GetXaxis()->SetNdivisions(9,0,0,kTRUE);
            g->GetYaxis()->SetTitleSize(0.05);
            g->GetXaxis()->SetTitleSize(0.05);
	    g->GetYaxis()->SetTitleOffset(0.95);

            g->Draw("APL");
            first = false;
        }
        else
        {
            g->Draw("PL SAME");
        }

        leg->AddEntry(g, Form("vth = %d", vth), "pl");
    }

    leg->Draw();

    c->SetGrid();
    c->SaveAs(Form("MeanDeltaT_vs_REF_vthScan_bar%d.pdf", dutBar));
}
