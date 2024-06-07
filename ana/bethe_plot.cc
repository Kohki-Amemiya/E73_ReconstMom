//------------------------- HOW TO USE ---------------------------------------------------------
//
// plot dE/dx vs beta*gamma using Bethe formula and check if bethe calculation is correct
//
// 2024/06/06 K.Amemiya
//-----------------------------------------------------------------------------------------------

#include <iostream>
#include <limits>
#include <cmath>
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"

double m = 0.93827203; // proton mass [GeV/c^2]

double Bethe(double beta, double Z, double A, double charge, double I, double C0, double a, double n, double X1, double X0) // calculate energy loss from bethe formula
{
    // +++++++++++++++++++++++++++++++++++++++++++++++++
    //     calculate energy loss rate [GeV/cm]
    //     arguments： beta, Z, A, charge, I[eV], C0, a, n, X1, X0 (cf. LEO p26)
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
    double meM = 0.000511 / m;                                                                             // me/M
    double Tmax = 1.02 * pow(10, -3) * pow(gamma, 2) * pow(beta, 2) / (1 + 2 * gamma * meM + pow(meM, 2)); // [GeV]
    EnergyLoss_rate = 1000 * 1.535 * pow(10, -4) * (Z / A) * pow(charge / beta, 2) * (std::log(1.02 * pow(10, -3) * pow(gamma, 2) * pow(beta, 2) * Tmax / pow((I * pow(10, -9)), 2)) - 2 * pow(beta, 2) - delta - 2 * C / Z);
    return EnergyLoss_rate; // [MeV/g*cm^2]
}

void bethe_plot()
{
    const int nPoints = 1000;
    double etaMin = 0.1;
    double etaMax = 10000.;
    double logEtaMin = log10(etaMin);
    double logEtaMax = log10(etaMax);

    TGraph *graph_Al = new TGraph(nPoints);
    TGraph *graph_Pb = new TGraph(nPoints);

    double minDeltaE_Al = std::numeric_limits<double>::max();
    double etaAtMinDeltaE_Al = etaMin;

    double minDeltaE_Pb = std::numeric_limits<double>::max();
    double etaAtMinDeltaE_Pb = etaMin;

    for (int i = 0; i < nPoints; ++i)
    {
        double logEta = logEtaMin + i * (logEtaMax - logEtaMin) / (nPoints - 1);
        double eta = pow(10, logEta);
        double beta = eta / sqrt(1 + pow(eta, 2));
        double deltaE_Al = Bethe(beta, 13, 26.98, 1, 166, -4.24, 0.0802, 3.63, 3.01, 0.1708); // dE at Al [MeV/g*cm^2]
        double deltaE_Pb = Bethe(beta, 82, 207.2, 1, 823, -6.20, 0.0936, 3.16, 3.81, 0.3776); // dE at Pb [MeV/g*cm^2]

        if (deltaE_Al < minDeltaE_Al)
        {
            minDeltaE_Al = deltaE_Al;
            etaAtMinDeltaE_Al = eta;
        }

        if (deltaE_Pb < minDeltaE_Pb)
        {
            minDeltaE_Pb = deltaE_Pb;
            etaAtMinDeltaE_Pb = eta;
        }

        graph_Al->SetPoint(i, eta, deltaE_Al);
        graph_Pb->SetPoint(i, eta, deltaE_Pb);
    }

    double betaAtEtaMax_Al = 10000 / sqrt(1 + 10000 * 10000);
    double deltaEAtEtaMax_Al = Bethe(betaAtEtaMax_Al, 13, 26.98, 1, 166, -4.24, 0.0802, 3.63, 3.01, 0.1708);

    double betaAtEtaMax_Pb = 10000 / sqrt(1 + 10000 * 10000);
    double deltaEAtEtaMax_Pb = Bethe(betaAtEtaMax_Pb, 82, 207.2, 1, 823, -6.20, 0.0936, 3.16, 3.81, 0.3776);

    TCanvas *canvas = new TCanvas("canvas", "Bethe Test", 800, 600);
    graph_Al->SetTitle("Energy Loss at Al and Pb;#beta #gamma;dE/dx [MeV/g cm^2]");
    graph_Al->GetXaxis()->CenterTitle();
    graph_Al->GetYaxis()->CenterTitle();
    graph_Al->SetMaximum(10);
    graph_Al->SetMinimum(1);
    graph_Al->SetLineColor(kBlue);
    graph_Pb->SetLineColor(kGreen);
    graph_Al->Draw("AL");
    graph_Pb->Draw("same");
    canvas->SetLogx();
    canvas->SetLogy();
    canvas->SetGrid();

    TPaveText *pt = new TPaveText(0.55, 0.65, 0.85, 0.85, "NDC");
    pt->SetFillColor(0);

    // Add Al text in blue
    TText *textAl1 = pt->AddText(Form("Al dE/dx min : %.2f [MeV], eta = %.2f", minDeltaE_Al, etaAtMinDeltaE_Al));
    textAl1->SetTextColor(kBlue);
    TText *textAl2 = pt->AddText(Form("Al dE/dx edge: %.2f [MeV], eta = 10000", deltaEAtEtaMax_Al));
    textAl2->SetTextColor(kBlue);

    // Add Pb text in green
    TText *textPb1 = pt->AddText(Form("Pb dE/dx min : %.2f [MeV], eta = %.2f", minDeltaE_Pb, etaAtMinDeltaE_Pb));
    textPb1->SetTextColor(kGreen);
    TText *textPb2 = pt->AddText(Form("Pb dE/dx edge: %.2f [MeV], eta = 10000", deltaEAtEtaMax_Pb));
    textPb2->SetTextColor(kGreen);

    pt->Draw();

    canvas->Update();
    canvas->SaveAs("Bethe_plot.png");
}
