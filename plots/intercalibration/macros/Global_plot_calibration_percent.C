void plot_calibration_diff_percent() {

    // ===== PATH =====
    std::string path = "/eos/home-c/cgiraldi/MTD/TB_CERN_Sep25/Lab5015Analysis/plots/intercalibration/macros/";

    std::string file1 = path + "Global_Calib_factors_Carlo.txt";
    std::string file2 = path + "Global_Calib_factors_Simona.txt";

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
    double val;

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





    // ===== DATA =====
    const int nBars = 16;

    double xL[nBars], yL[nBars];
    double xR[nBars], yR[nBars];

    std::vector<double> vecL, vecR;

    for (int i = 0; i < nBars; ++i) {

        if (!calib1.count({i,'L'}) || !calib2.count({i,'L'}) ||
            !calib1.count({i,'R'}) || !calib2.count({i,'R'})) {
            std::cout << "Missing bar " << i << std::endl;
            continue;
        }

        double c1_L = calib1[{i,'L'}];
        double c2_L = calib2[{i,'L'}];

        double c1_R = calib1[{i,'R'}];
        double c2_R = calib2[{i,'R'}];

        xL[i] = i;
        xR[i] = i;

        // ===== DELTA PERCENTUALE =====
        yL[i] = (c2_L - c1_L) / c1_L * 100.0;
        yR[i] = (c2_R - c1_R) / c1_R * 100.0;

        vecL.push_back(yL[i]);
        vecR.push_back(yR[i]);
    }

    // ===== STATS =====
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

    // ===== GRAPHS =====
    TGraph* gL = new TGraph(nBars, xL, yL);
    TGraph* gR = new TGraph(nBars, xR, yR);

    gL->SetMarkerStyle(20);
    gL->SetMarkerColor(kRed+1);
    gL->SetLineColor(kRed+1);
    gL->SetLineWidth(2);

    gR->SetMarkerStyle(21);
    gR->SetMarkerColor(kBlue+1);
    gR->SetLineColor(kBlue+1);
    gR->SetLineWidth(2);

    // ===== CANVAS =====
    TCanvas* c = new TCanvas("c", "Delta % TOFHIR Calib", 900, 700);
    c->SetGridy();

    gL->SetTitle(";Bar index;#Delta TOFHIR Calib [%]");

    // range consigliato
    gL->SetMinimum(-3);
    gL->SetMaximum(3);

    gL->Draw("ALP");
    gR->Draw("LP SAME");

    // ===== LINEA ZERO =====
    TLine* linezero = new TLine(0, 0, nBars-1, 0);
    linezero->SetLineStyle(2);
    linezero->SetLineWidth(2);
    linezero->Draw("SAME");

    // ===== LEGEND =====
    TLegend* leg = new TLegend(0.6, 0.75, 0.88, 0.88);
    leg->SetTextSize(0.03);
    leg->AddEntry(gL, Form("Side L (mean=%.2f%%)", meanL), "lp");
    leg->AddEntry(gR, Form("Side R (mean=%.2f%%)", meanR), "lp");

    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->Draw();

    c->Update();
}
