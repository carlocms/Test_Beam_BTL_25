void plot_calibration_diff() {

    // ===== PATH =====
    std::string path1 = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/intercalibration/module_intercalib/";
    std::string path2 = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/intercalibration/macros/";

    std::string file1 = path1 + "TOFHIR_LO_calibration_factors_Vov3.00_th12_global.txt";
    std::string file2 = path2 + "TOFHIR_intercalibration_factors_SIMONA.txt";

    // ===== MAPS =====
    std::map<std::pair<int,char>, double> calib1;
    std::map<std::pair<int,char>, double> calib2;

    std::string line;

    // ===== READ FILE 1 =====
    std::ifstream fin1(file1);
    while (std::getline(fin1, line)) {
        if (line.empty() || line[0] == '#') continue;

        std::stringstream ss(line);
        int bar;
        char side;
        double val;

        ss >> bar >> side >> val;

        calib1[{bar, side}] = val;

        std::cout << "[FILE1] bar=" << bar << " side=" << side << " val=" << val << std::endl;
    }
    fin1.close();


    // ===== READ FILE 2 =====
std::ifstream fin2(file2);

if (!fin2.is_open()) {
    std::cerr << "ERROR: cannot open file2: " << file2 << std::endl;
    return;
}

std::string line2;

while (std::getline(fin2, line2)) {

    if (line2.empty()) continue;

    std::stringstream ss(line2);

    int bar;
    std::string side_str;
    double vov, th, val;

    // parsing robusto
    if (!(ss >> bar >> side_str >> val)) {
        std::cerr << "ERROR parsing line2: " << line2 << std::endl;
        continue;
    }

    // sicurezza sul lato
    if (side_str.size() == 0) {
        std::cerr << "ERROR: empty side in line2: " << line2 << std::endl;
        continue;
    }

    char side = side_str[0];

    // ===== DEBUG COMPLETO =====
    std::cout << "[FILE2 RAW] " << line2 << std::endl;
    std::cout << "[FILE2 PARSED] bar=" << bar
              << " side=" << side
              << " val=" << val << std::endl;

    // ===== EVENTUALE FILTER =====
    //if (std::abs(th - 12.0) > 1e-6) continue;

    calib2[{bar, side}] = val;
}

fin2.close();

// ===== CHECK FINALE =====
std::cout << "\n[SUMMARY] Entries in calib2 = " << calib2.size() << std::endl;


    // ===== CREATE GRAPHS =====
    const int nBars = 16;

    double xL[nBars], yL[nBars];
    double xR[nBars], yR[nBars];

    std::vector<double> vecL, vecR;

    for (int i = 0; i < nBars; ++i) {

        xL[i] = i;
        xR[i] = i;

        if (!calib1.count({i,'L'}) || !calib2.count({i,'L'}) ||
            !calib1.count({i,'R'}) || !calib2.count({i,'R'})) {

            std::cout << "WARNING missing entry for bar " << i << std::endl;
            continue;
        }

        double c1_L = calib1[{i,'L'}];
        double c2_L = calib2[{i,'L'}];

        double c1_R = calib1[{i,'R'}];
        double c2_R = calib2[{i,'R'}];

        yL[i] = c2_L - c1_L;
        yR[i] = c2_R - c1_R;

        vecL.push_back(yL[i]);
        vecR.push_back(yR[i]);

        // ===== DEBUG DIFFERENZE =====
        std::cout << "[DIFF] bar=" << i
                  << " L: " << c2_L << " - " << c1_L << " = " << yL[i]
                  << " | R: " << c2_R << " - " << c1_R << " = " << yR[i]
                  << std::endl;
    }

    // ===== FUNCTION MEAN / RMS =====
    auto computeStats = [](std::vector<double>& v, double& mean, double& rms) {
        mean = 0;
        rms = 0;

        for (auto val : v) mean += val;
        mean /= v.size();

        for (auto val : v) rms += (val - mean)*(val - mean);
        rms = std::sqrt(rms / v.size());
    };

    double meanL, rmsL, meanR, rmsR;
    computeStats(vecL, meanL, rmsL);
    computeStats(vecR, meanR, rmsR);

    std::cout << "\n=== STATS ===" << std::endl;
    std::cout << "L: mean=" << meanL << " RMS=" << rmsL << std::endl;
    std::cout << "R: mean=" << meanR << " RMS=" << rmsR << std::endl;

    // ===== GRAPHS =====
    TGraph* gL = new TGraph(nBars, xL, yL);
    TGraph* gR = new TGraph(nBars, xR, yR);

    gL->SetMarkerStyle(20);
    gL->SetMarkerColor(kRed);
    gL->SetLineColor(kRed);

    gL->SetMinimum(-0.03);
    gL->SetMaximum(0.03);

    gL->GetXaxis()->SetRangeUser(-0.5, 15.5);

    gR->SetMarkerStyle(21);
    gR->SetMarkerColor(kBlue);
    gR->SetLineColor(kBlue);

    // ===== CANVAS =====
    TCanvas* c = new TCanvas("c", "Calibration differences", 800, 600);
    c->SetGridy();

    gL->SetTitle("Calibration difference;Bar index;#Delta TOFHIR calib ");
    gL->Draw("ALP");
    gR->Draw("LP SAME");

    gL->SetLineWidth(2);
    gR->SetLineWidth(2);

    // ===== AXIS SETTINGS =====
    // forza tick solo su interi (0,1,2,...)
    //gL->GetXaxis()->SetLimits(-0.5, nBars - 0.5);
    //gL->GetXaxis()->SetNdivisions(nBars, kFALSE);  // solo interi

    // ===== LEGEND CON STATS =====
    TLegend* leg = new TLegend(0.6, 0.7, 0.88, 0.88);
    leg->SetTextSize(0.03);
    leg->AddEntry(gL, Form("Side L (mean=%.5f, RMS=%.5f)", meanL, rmsL), "lp");
    leg->AddEntry(gR, Form("Side R (mean=%.5f, RMS=%.5f)", meanR, rmsR), "lp");

    leg->Draw();

    // ===== SAVE PDF =====
    c->SaveAs("TOFHIR_calibration_diff.pdf");

    c->Update();
}
