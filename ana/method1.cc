//++++++++++++++++++ How to Use +++++++++++++++++++++
//
// compare beta loss with bethe calculation and geant data, and calculate initial P(mom at spectrometer).
//
// 2024/05/09 K.Amemiya
//+++++++++++++++++++++++++++++++++++++++++++++++++++

double c_ns = 0.299792458; // [m/ns]
double m = 0.93827203;     // [GeV/c^2]
double resgtime[18];       // [ns]
double gtime[18];          // [ns]
double adc[18];            // [mV]
double l1 = 7.495;         // spectrometer -> Al length [m]
double l3 = 7.440;         // Scintillator -> CVC length [m]
int RUN_NUM = 22;          // run number

double Bethe(double beta, double Z, double A, double charge, double density, double I, double C0, double a, double n, double X1, double X0) // calculate energy loss from bethe formula
{
    // +++++++++++++++++++++++++++++++++++++++++++++++++
    //     calculate energy loss rate [GeV/cm]
    //     arguments： beta, Z, A, charge, density[g/cm^3], I[eV], C0, a, n, X1, X0 (cf. LEO p26)
    // +++++++++++++++++++++++++++++++++++++++++++++++++
    double EnergyLoss_rate = 0.;
    double gamma = 1 / sqrt(1 - pow(beta, 2));
    double eta = beta * gamma;
    double C = 0.;
    if (eta > 0.1)
    {
        C = (0.4223 * pow(eta, -2) + 0.03040 * pow(eta, -4) - 0.0003810 * pow(eta, -6)) * pow(10, -6) * pow(I, 2) + (3.850 * pow(eta, -2) - 0.1667 * pow(eta, -4) + 0.001579 * pow(eta, -6)) * pow(10, -9) * pow(I, 3);
    }
    double X = log10(eta);
    double delta = 0.;
    if (X0 < X && X < X1)
    {
        delta = 4.6052 * X + C0 + a * pow((X1 - X), n);
    }
    if (X > X1)
    {
        delta = 4.6052 * X + C0;
    }
    double meM = 0.000511 / m;                                                                                                                                                                                                   // me/M
    double Tmax = 1.02 * pow(10, -3) * pow(gamma, 2) * pow(beta, 2) / (1 + 2 * gamma * meM + pow(meM, 2));                                                                                                                       // [GeV]
    EnergyLoss_rate = density * 1.535 * pow(10, -4) * (Z / A) * pow(charge / beta, 2) * (std::log(1.02 * pow(10, -3) * pow(gamma, 2) * pow(beta, 2) * Tmax / pow((I * pow(10, -9)), 2)) - 2 * pow(beta, 2) - delta - 2 * C / Z); // [GeV/cm]

    return EnergyLoss_rate; // [GeV/cm]
}

void method1()
{
    TFile *file = new TFile(Form("/home/had/kohki/work/geant/e73/rootfiles/E73simulator_run%03d.root", RUN_NUM));

    TTree *tree = (TTree *)file->Get("tree");

    TCanvas *c1 = new TCanvas("c1", "deltaE distribution");
    TCanvas *c2 = new TCanvas("c2", "P1 distribution");
    c1->SetGrid();
    c2->SetGrid();

    int nEntries = tree->GetEntries();

    // beta drop from VP information(h1) and Bethe calculation(h2), and corrected p1(h3)

    TH1D *hene1 = new TH1D("hene", "energy loss", 1500, 0, 15);            // bin width=0.1
    TH1D *hene2 = new TH1D("hene", "energy loss", 1500, 0, 15);            // bin width=0.1
    TH1D *hene1_vp = new TH1D("hene_vp", "energy loss", 1500, 0, 15);      // bin width=0.1
    TH1D *hene2_vp = new TH1D("hene_vp", "energy loss", 1500, 0, 15);      // bin width=0.1
    TH1D *hmom = new TH1D("hmom", "reconstructed P_{1} ", 6000, 0.7, 1.3); // bin width=1MeV/c
    TLegend *legend = new TLegend(0.78, 0.64, 0.98, 0.77);

    tree->SetBranchAddress("resgtime", resgtime);
    tree->SetBranchAddress("gtime", gtime);
    tree->SetBranchAddress("adc", adc);

    // calculate momentum drop in 2 ways and fill histograms
    for (int j = 0; j < nEntries; ++j)
    {
        tree->GetEntry(j);
        if (gtime[0] == 0 || gtime[5] == 0 || gtime[6] == 0) // exclude no hit events
            continue;

        // double beta1_vp = l1 / (c_ns * gtime[0]);
        // double beta2_vp = l3 / (c_ns * (gtime[6] - gtime[5]));
        // double delta_beta_vp = beta1_vp - beta2_vp;

        double beta_input = (l1 + l3) / (c_ns * resgtime[6]); // beta for bethe input
        // double beta_input = beta1_vp;                                             // test
        double p_input = m * beta_input / sqrt(1 - pow(beta_input, 2));                                      // mom at bethe input
        double deltaE1 = Bethe(beta_input, 13, 26.98, -1, 2.70, 166, -4.24, 0.0802, 3.63, 3.01, 0.1708) * 1; // energy loss at 1cm Al [GeV]
        double deltaP1 = deltaE1 * sqrt((pow(m, 2)) + pow(p_input, 2)) / p_input;                            // [GeV/c]
        double delta_beta1 = deltaP1 * pow(m, 2) / pow((pow(m, 2) + pow(p_input, 2)), 3 / 2);
        double p2_ = p_input - deltaP1; //[GeV/c]
        double beta2 = beta_input - delta_beta1;
        double deltaE2 = Bethe(beta2, 5.575, 11.07, -1, 1.032, 64.7, -3.20, 0.1610, 3.24, 2.49, 0.1464) * 3; // energy loss at 3cm Scintillator [GeV]
        double deltaP2 = deltaE2 * sqrt((pow(m, 2)) + pow(p2_, 2)) / p2_;
        double delta_beta2 = deltaP2 * pow(m, 2) / pow((pow(m, 2) + pow(p2_, 2)), 3 / 2);
        double delta_beta = delta_beta1 + delta_beta2; // total beta loss at Al and Scintillator

        double tof = resgtime[6]; //[ns]
        double a_ = c_ns * tof;
        double b_ = c_ns * tof * delta_beta + l1 + l3;
        double c_ = l1 * delta_beta;
        double beta1 = (b_ + sqrt(pow(b_, 2) - 4 * a_ * c_)) / (2 * a_);
        double p1 = m * beta1 / sqrt(1 - pow(beta1, 2)); // [GeV/c]
        // std::cout << "delta_beta_vp: " << delta_beta_vp << " delta_beta: " << delta_beta << " p1: " << p1 << std::endl;

        hene1->Fill(deltaE1 * 1000); //[MeV]
        hene2->Fill(deltaE2 * 1000); //[MeV]
        hene1_vp->Fill(adc[1]);
        hene2_vp->Fill(adc[4]);
        hmom->Fill(p1);
        // std::cout << "deltaE1: " << deltaE1 << ", deltaE2: " << deltaE2 << std::endl;
        // std::cout << "deltaE1_vp: " << adc[1] * 0.001 << ", deltaE2_vp: " << adc[4] * 0.001 << std::endl;
    }
    // fitting and get mpv
    double max_gaus = hmom->GetMean() + 2 * hmom->GetRMS();
    double min_gaus = hmom->GetMean() - 2 * hmom->GetRMS();
    TF1 *gaus = new TF1("gaus", "gaus", min_gaus, max_gaus);
    hmom->Fit("gaus", "", "", min_gaus, max_gaus);

    c1->cd();
    legend->AddEntry(hene1, "bethe", "l");
    hene1->SetXTitle("#Delta E (MeV)");
    hene1->GetXaxis()->CenterTitle();
    hene1->SetYTitle("count");
    hene1->GetYaxis()->CenterTitle();
    hene1->SetLineColor(4);
    hene1->SetLineWidth(2);
    hene1->Draw();
    hene2->SetLineColor(4);
    hene2->SetLineWidth(2);
    hene2->Draw("same");
    legend->AddEntry(hene1_vp, "VP", "l");
    hene1_vp->SetLineColor(1);
    hene1_vp->SetLineWidth(2);
    hene1_vp->Draw("same");
    hene2_vp->SetLineColor(1);
    hene2_vp->SetLineWidth(2);
    hene2_vp->Draw("same");
    legend->Draw();

    c2->cd();
    hmom->SetXTitle("P_{1} (GeV/C)");
    hmom->GetXaxis()->CenterTitle();
    hmom->SetYTitle("count");
    hmom->GetYaxis()->CenterTitle();
    hmom->SetLineWidth(2);
    hmom->Draw();

    std::cout << "deltaE1 mean: " << hene1->GetMean() << ", deltaE2: " << hene2->GetMean() << std::endl;
    std::cout << "deltaE1_vp: " << hene1_vp->GetMean() << ", deltaE2_vp: " << hene2_vp->GetMean() << std::endl;

    std::cout << "***********************************" << std::endl
              << " reconstructed P1 : " << gaus->GetParameter(1) << " GeV/c ± " << gaus->GetParameter(2) << "GeV/c" << std::endl
              << "***********************************" << std::endl;

    c1->SaveAs(Form("../img/bethe_mom/Eloss_%03d.png", RUN_NUM));
    c2->SaveAs(Form("../img/bethe_mom/p1_%03d.png", RUN_NUM));
}
