void plot_two_files_comparison() {

    // ===== FILES =====
    std::string file1 = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/intercalibration/macros/Global_Calib_factors_Carlo.txt";
    std::string file2 = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/intercalibration/macros/Global_Calib_factors_Simona.txt";

    const int nBars = 16;

    // ===== FUNZIONE PER LEGGERE FILE =====
    auto readFile = [&](std::string filename,
                        double* xL, double* yL,
                        double* xR, double* yR) {

        std::ifstream fin(filename);

        if (!fin.is_open()) {
            std::cerr << "ERROR opening " << filename << std::endl;
            return;
        }

        std::string line;

        while (std::getline(fin, line)) {

            if (line.empty()) continue;

            // rimuove eventuali \r
            line.erase(std::remove(line.begin(), line.end(), '\r'), line.end());

            std::stringstream ss(line);

            int bar;
            std::string side_str;
            double val;

            if (!(ss >> bar >> side_str >> val)) {
                std::cerr << "Parse error: " << line << std::endl;
                continue;
            }

            char side = side_str[0];

            if (bar < 0 || bar >= nBars) continue;

            if (side == 'L') {
                xL[bar] = bar;
                yL[bar] = val;
            }
            else if (side == 'R') {
                xR[bar] = bar;
                yR[bar] = val;
            }

            std::cout << "[" << filename << "] bar=" << bar
                      << " side=" << side
                      << " val=" << val << std::endl;
        }

        fin.close();
    };

    // ===== ARRAYS =====
    double xL1[nBars], yL1[nBars], xR1[nBars], yR1[nBars];
    double xL2[nBars], yL2[nBars], xR2[nBars], yR2[nBars];

    // ===== READ FILES =====
    readFile(file1, xL1, yL1, xR1, yR1);
    readFile(file2, xL2, yL2, xR2, yR2);
/*
    // ===== FUNZIONE PER PLOT =====
    auto makePlot = [&](double* xL, double* yL,
                        double* xR, double* yR,
                        std::string canvasName,
                        std::string outName) {

        TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
        c->SetGridy();

        TGraph* gL = new TGraph(nBars, xL, yL);
        TGraph* gR = new TGraph(nBars, xR, yR);

        // stile
        gL->SetMarkerStyle(20);
        gL->SetMarkerColor(kRed);
        gL->SetLineColor(kRed);
        gL->SetLineWidth(2);

        gR->SetMarkerStyle(21);
        gR->SetMarkerColor(kBlue);
        gR->SetLineColor(kBlue);
        gR->SetLineWidth(2);

        // assi
        gL->SetTitle(";Bar index;Calibration value");
        gL->GetXaxis()->SetLimits(-0.5, nBars - 0.5);
        gL->GetXaxis()->SetNdivisions(nBars, kFALSE);

        gL->Draw("ALP");
        gR->Draw("LP SAME");

        // legenda
        TLegend* leg = new TLegend(0.7, 0.75, 0.88, 0.88);
        leg->SetTextSize(0.04);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);

        leg->AddEntry(gL, "Side L", "lp");
        leg->AddEntry(gR, "Side R", "lp");

        leg->Draw();

        // linea opzionale a y=1
        TLine* line = new TLine(0, 1.0, nBars-1, 1.0);
        line->SetLineStyle(2);
        line->Draw("SAME");

        // salva
        c->SaveAs(outName.c_str());

        c->Update();
    };
*/

    auto makePlot = [&](double* xL, double* yL,
                    double* xR, double* yR,
                    std::string canvasName,
                    std::string outName) {

    TCanvas* c = new TCanvas(canvasName.c_str(), canvasName.c_str(), 800, 600);
    c->SetGridy();

    TGraph* gL = new TGraph(nBars, xL, yL);
    TGraph* gR = new TGraph(nBars, xR, yR);

    // ===== STATS =====
    auto computeStats = [&](double* y, double& mean, double& rms) {
        mean = 0;
        rms  = 0;

        for (int i = 0; i < nBars; ++i) mean += y[i];
        mean /= nBars;

        for (int i = 0; i < nBars; ++i)
            rms += (y[i] - mean)*(y[i] - mean);

        rms = std::sqrt(rms / nBars);
    };

    double meanL, rmsL, meanR, rmsR;
    computeStats(yL, meanL, rmsL);
    computeStats(yR, meanR, rmsR);

    // ===== STYLE =====
    gL->SetMarkerStyle(20);
    gL->SetMarkerColor(kRed);
    gL->SetLineColor(kRed);
    gL->SetLineWidth(2);

    gR->SetMarkerStyle(21);
    gR->SetMarkerColor(kBlue);
    gR->SetLineColor(kBlue);
    gR->SetLineWidth(2);

    // ===== AXES =====
    gL->SetTitle(";Bar index;TOFHIR #times LO calib ");
    gL->SetMinimum(0.5);
    gL->SetMaximum(1.5);
    //gL->GetXaxis()->SetLimits(-0.5, nBars - 0.5);
    //gL->GetXaxis()->SetNdivisions(nBars, kFALSE);

    gL->Draw("ALP");
    gR->Draw("LP SAME");

    // ===== LEGEND =====
    TLegend* leg = new TLegend(0.55, 0.7, 0.88, 0.88);

    leg->SetTextSize(0.04);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);

    leg->AddEntry(gL, Form("Side L (#mu=%.3f, RMS=%.3f)", meanL, rmsL), "lp");
    leg->AddEntry(gR, Form("Side R (#mu=%.3f, RMS=%.3f)", meanR, rmsR), "lp");

    leg->Draw();

    // ===== LINEA DI RIFERIMENTO =====
    TLine* line = new TLine(0, 1.0, nBars-1, 1.0);
    line->SetLineStyle(2);
    line->Draw("SAME");

    // ===== SAVE =====
    c->SaveAs(outName.c_str());

    c->Update();
};


    // ===== CREA I DUE PLOT =====
    makePlot(xL1, yL1, xR1, yR1, "file1_plot", "file1_plot.pdf");
    makePlot(xL2, yL2, xR2, yR2, "file2_plot", "file2_plot.pdf");
}
