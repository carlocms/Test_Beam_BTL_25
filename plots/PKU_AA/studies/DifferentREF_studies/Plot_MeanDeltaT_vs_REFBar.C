#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <string>

#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TAxis.h"

void Plot_MeanDeltaT_vs_REFBar(int vth)
{
    // -------------------------------
    // --- Base directory
    // -------------------------------
    std::string baseDir = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/DifferentREF_studies/";

    // REF bars disponibili
    std::vector<int> refBars = {4,5,6,7,8,9,10,11};

    // -------------------------------
    // --- Map: DUT bar -> graph
    // -------------------------------
    std::map<int, TGraphErrors*> graphs;

    for(int i = 0; i < 16; i++)
        graphs[i] = new TGraphErrors();

    // -------------------------------
    // --- Loop sulle REF bars
    // -------------------------------
    for(int refBar : refBars)
    {
        std::string fileName = baseDir +
            "REF_Bar" + std::to_string(refBar) +
            "/outfileFile_vth" + std::to_string(vth) + ".txt";

        std::ifstream fin(fileName);

        if(!fin.is_open())
        {
            std::cout << "Missing file: " << fileName << std::endl;
            continue;
        }

        std::cout << "Reading: " << fileName << std::endl;

        int bar;
        double mean, sigma;

        // skip header
        std::string line;
        std::getline(fin, line);

        while(fin >> bar >> mean >> sigma)
        {
            int ip = graphs[bar]->GetN();

            graphs[bar]->SetPoint(ip, refBar, mean);
            graphs[bar]->SetPointError(ip, 0., sigma);
        }

        fin.close();
    }

    // -------------------------------
    // --- Canvas
    // -------------------------------
    TCanvas* c = new TCanvas("c_REFstudy","Mean #DeltaT vs REF bar",1200,800);

    // colori
    int colors[] = {
        kRed, kBlue, kGreen+2, kMagenta, kCyan+2,
        kOrange+1, kBlack, kViolet, kAzure+2,
        kPink+2, kTeal+2, kSpring+5,
        kGray+2, kYellow+2, kBlue+3, kRed+3
    };

    bool first = true;

    // -------------------------------
    // --- Draw
    // -------------------------------
    for(int bar = 0; bar < 16; bar++)
    {
        auto g = graphs[bar];

        g->SetMarkerStyle(20);
        g->SetMarkerColor(colors[bar % 16]);
        g->SetLineColor(colors[bar % 16]);

        if(first)
        {
            g->SetTitle(";REF bar;Mean #Delta T (L-R) [ps]");
            g->GetXaxis()->SetLimits(3,12);
            g->GetXaxis()->SetNdivisions(9, 0, 0, kTRUE);	    
            g->GetYaxis()->SetTitleSize(0.05);
            g->GetXaxis()->SetTitleSize(0.05);
	    g->GetYaxis()->SetRangeUser(-1500.,1500.);
	    g->GetYaxis()->SetTitleOffset(0.95);

            g->Draw("APL");
            first = false;
        }
        else
        {
            g->Draw("PL SAME");
        }
    }

    // -------------------------------
    // --- Legend
    // -------------------------------
    TLegend* leg = new TLegend(0.65, 0.60, 0.90, 0.90); // alto dx

    leg->SetNColumns(2);  // <-- due colonne

    //leg->SetHeader(Form("#bf{vth = %d}", vth), "C");

    // stile
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(1001);
    leg->SetBorderSize(1);
    leg->SetTextSize(0.025);
    leg->SetTextFont(42);

    // entries
    for(int bar = 0; bar < 16; bar++)
    {
       leg->AddEntry(graphs[bar], Form("DUT bar %d", bar), "pl");
    }

    leg->Draw();
    
    c->SetGrid();

    // -------------------------------
    // --- Save
    // -------------------------------
    c->SaveAs(Form("MeanDeltaT_vs_REFBar_vth%d.pdf", vth));
}

