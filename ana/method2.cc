//++++++++++++++++++ How to Use +++++++++++++++++++++
//
// - Scan beta between 0.600 < beta < 0.850 and find the beta
//   which gives the minimum difference between the calculated ToF
//   and the measured ToF.
//
// - Using this beta, calculate the momentum at the spectrometer.
//
//   2024/06/07 K.Amemiya
//+++++++++++++++++++++++++++++++++++++++++++++++++++

#include <iostream>
#include <limits>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"

double c_ns = 0.299792458; // [m/ns]
double m = 0.93827203;     // [GeV/c^2]
double resgtime[18];       // [ns]
double adc[18];            // [MeV]
int RUN_NUM = 26;          // run number
TFile *file = new TFile(Form("/home/had/kohki/work/geant/e73/rootfiles/E73simulator_run%03d.root", RUN_NUM));
TTree *tree = (TTree *)file->Get("tree");
int nEntries = tree->GetEntries();

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

double ReconstMom()
{
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //    scan beta and find the most likely initial beta at the spectrometer
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    // gROOT->SetStyle("ATLAS");

    TFile *file = new TFile(Form("/home/had/kohki/work/geant/e73/rootfiles/E73simulator_run%03d.root", RUN_NUM));
    TTree *tree = (TTree *)file->Get("tree");
    tree->SetBranchAddress("resgtime", resgtime);

    TH1D *htof = new TH1D("htof", "", 3000, 50, 80); // bin width=10ps

    for (int k = 0; k < nEntries; ++k)
    {
        tree->GetEntry(k);
        if (resgtime[0] == 0 || resgtime[3] == 0 || resgtime[6] == 0) // exclude no hit events
            continue;
        double tof_cvc = resgtime[6]; //[ns]
        htof->Fill(tof_cvc);
    }
    double tof_cvc_mean = htof->GetMean();

    TGraph *graph = new TGraph();
    double diff_min = std::numeric_limits<double>::max(); // initialize with biggest value for <double>
    double beta_init = 0.;                                // estimated input beta

    for (int i = 0; i < 2500; ++i) // scan between 0.6000 < beta < 0.8500
    {
        double beta = 0.8500 - (i * 0.0001);
        double beta_ = 0.8500 - (i * 0.0001); // beta for calculation
        double tof_sum = 0.;

        // calculate sum of ToF
        for (int j = 0; j < 15000; ++j)
        {
            double tof = 0.001 / (beta_ * c_ns); // [ns]
            tof_sum += tof;

            double deltaE = 0.;
            if (j > 7499 && j < 7510)
                deltaE = Bethe(beta_, 13, 26.98, -1, 2.70, 166, -4.24, 0.0802, 3.63, 3.01, 0.1708) * 0.1; // dE at Al [GeV]
            if (j > 7519 && j < 7550)
                deltaE = Bethe(beta_, 5.575, 11.07, -1, 1.032, 64.7, -3.20, 0.1610, 3.24, 2.49, 0.1464) * 0.1; // dE at Scintillator [GeV]
            else
                deltaE = Bethe(beta_, 7.372, 14.80, -1, 0.001205, 85.7, -10.6, 0.1091, 3.40, 4.28, 1.742) * 0.1; // dE at Air [GeV]

            double p = m * beta_ / sqrt(1 - pow(beta_, 2));
            double deltaP = deltaE * sqrt((pow(m, 2)) + pow(p, 2)) / p;
            double delta_beta = deltaP * pow(m, 2) / pow((pow(m, 2) + pow(p, 2)), 3 / 2);
            beta_ = beta_ - delta_beta;
            p = p - deltaP;
        }

        // calculate average ToF difference
        double diff = pow((tof_cvc_mean - tof_sum), 2);
        graph->SetPoint(i, beta, diff);

        if (diff < diff_min)
        {
            diff_min = diff;
            beta_init = beta;
        }
        // std::cout << "Beta: " << beta << ", diff: " << diff << std::endl;
    }

    double p_init = m * beta_init / sqrt(1 - pow(beta_init, 2));
    std::cout << "Minimum diff: " << diff_min << ", Corresponding beta: " << beta_init << ", Reconstructed mom:" << p_init << std::endl;

    TCanvas *c1 = new TCanvas("c1", "diff dependency on beta");
    c1->SetGrid();
    graph->SetMarkerStyle(7);
    graph->SetTitle("#Beta scan");
    graph->GetXaxis()->SetTitle("initial #Beta");
    graph->GetXaxis()->CenterTitle();
    graph->GetYaxis()->SetTitle("Diff");
    graph->GetYaxis()->CenterTitle();
    graph->Draw("AP");

    // c1->SaveAs(Form("../img/ReconstMom/scan_run%03d.png", RUN_NUM));

    return beta_init;
}

void method2()
{
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++
    // check if the bethe calculation is correct
    // by looking at the energy loss at Al and Scintillator
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++

    // CALCULATE dE1 and dE2S
    double beta_input = ReconstMom();
    double deltaE1 = 0.;
    double deltaE2 = 0.;
    for (int j = 0; j < 7550; ++j)
    {
        double deltaE = 0.;
        if (j > 7499 && j < 7510)
        {
            deltaE = Bethe(beta_input, 13, 26.98, -1, 2.70, 166, -4.24, 0.0802, 3.63, 3.01, 0.1708) * 0.1; // dE at Al [GeV]
            deltaE1 += deltaE;
        }
        if (j > 7519 && j < 7550)
        {
            deltaE = Bethe(beta_input, 5.575, 11.07, -1, 1.032, 64.7, -3.20, 0.1610, 3.24, 2.49, 0.1464) * 0.1; // dE at Scintillator [GeV]
            deltaE2 += deltaE;
        }
        else
            deltaE = Bethe(beta_input, 7.372, 14.80, -1, 0.001205, 85.7, -10.6, 0.1091, 3.40, 4.28, 1.742) * 0.1; // dE at Air [GeV]
        // deltaE = 0;

        double p = m * beta_input / sqrt(1 - pow(beta_input, 2));
        double deltaP = deltaE * sqrt((pow(m, 2)) + pow(p, 2)) / p;
        double delta_beta = deltaP * pow(m, 2) / pow((pow(m, 2) + pow(p, 2)), 3 / 2);
        beta_input = beta_input - delta_beta;
        p = p - deltaP;
    }

    // GET dE1 and dE2 FROM VP
    TH1D *hene1_vp = new TH1D("hene1_vp", "", 1500, 0, 15); // bin width=10ps
    TH1D *hene2_vp = new TH1D("hene2_vp", "", 1500, 0, 15); // bin width=10ps
    TLegend *legend = new TLegend(0.78, 0.64, 0.98, 0.77);
    tree->SetBranchAddress("resgtime", resgtime);
    tree->SetBranchAddress("adc", adc);

    for (int k = 0; k < nEntries; ++k)
    {
        tree->GetEntry(k);
        if (resgtime[0] == 0 || resgtime[3] == 0 || resgtime[6] == 0) // exclude no hit events
            continue;
        hene1_vp->Fill(adc[1]);
        hene2_vp->Fill(adc[4]);
        // std::cout << "adc1: " << adc[1] << ", adc2: " << adc[4] << std::endl;
    }

    // Find the maximum values of the histograms
    double max1 = hene1_vp->GetMaximum();
    double max2 = hene2_vp->GetMaximum();
    double maxY = std::max(max1, max2);

    TCanvas *c2 = new TCanvas("c2", "energy loss distribution");
    c2->SetGrid();
    legend->AddEntry(hene1_vp, "VP", "l");
    hene1_vp->SetLineColor(kBlack);
    hene1_vp->SetLineWidth(2);
    hene1_vp->SetTitle("energy loss");
    hene1_vp->GetXaxis()->SetTitle("#Delta E (MeV)");
    hene1_vp->GetXaxis()->CenterTitle();
    hene1_vp->GetYaxis()->SetTitle("counts");
    hene1_vp->GetYaxis()->CenterTitle();
    hene1_vp->Draw();
    hene2_vp->SetLineColor(kBlack);
    hene2_vp->SetLineWidth(2);
    hene2_vp->Draw("same");

    TLine *line1 = new TLine(deltaE1 * 1000, 0, deltaE1 * 1000, maxY);
    line1->SetLineColor(kBlue);
    line1->SetLineStyle(2);
    line1->SetLineWidth(4);
    line1->Draw("same");

    TLine *line2 = new TLine(deltaE2 * 1000, 0, deltaE2 * 1000, maxY);
    line2->SetLineColor(kBlue);
    line2->SetLineStyle(2);
    line2->SetLineWidth(4);
    line2->Draw("same");
    legend->AddEntry(line2, "bethe", "l");

    legend->Draw();

    std::cout << "deltaE1: " << deltaE1 * 1000 << ", deltaE1_vp: " << hene1_vp->GetMean() << std::endl;
    std::cout << "deltaE2: " << deltaE2 * 1000 << ", deltaE2_vp: " << hene2_vp->GetMean() << std::endl;

    // c2->SaveAs(Form("../img/ReconstMom/dE_distribution_%03d.png", RUN_NUM));
}
