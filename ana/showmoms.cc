//+++++++++++++++++ How to Use +++++++++++++++++++++
//
// show tof and mom hist at each counters
//
// 2024/05/02 K.Amemiya
//+++++++++++++++++++++++++++++++++++++++++++++++++++

#include <algorithm>

void showmoms()
{
    TFile *file = new TFile("/home/had/kohki/work/geant/e73/rootfiles/E73simulator_run026.root");

    TTree *tree = (TTree *)file->Get("tree");

    TCanvas *ctof = new TCanvas("ctof", "Global Time (TOF) at CVC and NC", 800, 600);
    ctof->SetGrid();
    TCanvas *cmom = new TCanvas("cmom", "momentum at CVC and NC", 800, 600);
    cmom->SetGrid();
    TCanvas *cmom3 = new TCanvas("cmom3", "corrected momentum at CVC and NC", 800, 600);
    cmom3->SetGrid();

    TLegend *legend = new TLegend(0.7, 0.75, 0.9, 0.9);

    double c = 0.299792458;                     // [m/ns]
    double m = 0.93827203;                      // [GeV/c^2]
    double l[4] = {15.00, 15.84, 15.91, 15.98}; // [m] {position: cvc, nc1, nc2, nc3}
    double resgtime[18];                        // [ns]
    double gtime[18];                           // [ns]
    double mom[18];                             // [GeV/c]
    int nEntries = tree->GetEntries();

    // calculate tof/momentum at each cahnnel(CVC and NC) and fill histograms //
    for (int i = 0; i < 4; ++i)
    {

        TH1D *htof = new TH1D(Form("htof_%d", i), "ToF at CVC and NC", 2000, 60, 80);       // bin width=10ps
        TH1D *hmom = new TH1D(Form("hmom_%d", i), "momentum at CVC and NC", 400, 0.8, 1.2); // bin width=1MeV/c

        tree->SetBranchAddress("resgtime", resgtime);
        tree->SetBranchAddress("gtime", gtime);
        tree->SetBranchAddress("mom", mom);

        for (int j = 0; j < nEntries; ++j)
        {
            tree->GetEntry(j);

            if (resgtime[6] == 0 || resgtime[9] == 0 || resgtime[12] == 0 || resgtime[15] == 0) // only take all channel hit events
                continue;

            double tof_counter = gtime[(i + 2) * 3]; // gtime at CVC, nc1, nc2, nc3 (VP)
            // double tof_counter = resgtime[(i + 2) * 3]; // resgtime at CVC, nc1, nc2, nc3 (VP)
            double beta = l[i] / (c * tof_counter);
            // double mom_counter = m * beta / sqrt(1 - pow(beta, 2)); // [GeV/c]
            double mom_counter = mom[(i + 2) * 3] * 0.001; // [GeV] for test
            htof->Fill(tof_counter);
            hmom->Fill(mom_counter);

            // if (i == 0)
            //     std::cout << "channel: " << i << " beta: " << beta
            //               << " Mom: " << mom_counter << " GeV/c" << std::endl; // for test
        }

        // get sigma
        std::cout << "channel: " << i << " ToF sigma: " << htof->GetStdDev() << " ns"
                  << " Mom sigma: " << hmom->GetStdDev() << " GeV/c" << std::endl;
        // get mean
        std::cout << "channel: " << i << " ToF mean: " << htof->GetMean() << " ns"
                  << " Mom mean: " << hmom->GetMean() << " GeV/c" << std::endl;

        // Set histogram style
        htof->SetLineColor(i + 1);
        htof->SetLineWidth(2);
        htof->SetStats(0);
        hmom->SetLineColor(i + 1);
        hmom->SetLineWidth(2);
        hmom->SetStats(0);

        if (i == 0)
        {
            ctof->cd();
            htof->Draw();
            htof->SetXTitle("ToF (ns)");
            htof->SetYTitle("count");
            cmom->cd();
            hmom->Draw();
            hmom->SetXTitle("momentum (GeV/C)");
            hmom->SetYTitle("count");
        }
        else
        {
            ctof->cd();
            htof->Draw("same");
            cmom->cd();
            hmom->Draw("same");
        }

        legend->AddEntry(htof, Form("Channel %d", i), "l");
    }

    // calculate tof/momentum after material and fill histograms //

    TH1D *hmom3 = new TH1D("hmom3", "p_{3}", 600, 0.7, 1.3); // bin width=1MeV/c
    for (int j = 0; j < nEntries; ++j)
    {
        tree->GetEntry(j);
        if (resgtime[5] == 0 || resgtime[6] == 0) // take all channel hit events
            continue;

        double beta3 = 7.45 / (c * (resgtime[6] - resgtime[5])); // corrected beta at cvc
        double mom3 = m * beta3 / sqrt(1 - pow(beta3, 2));
        hmom3->Fill(mom3);
        // if (j == 0)
        //     std::cout << " beta3: " << beta3
        //               << " Mom3: " << mom3 << " GeV/c" << std::endl; // for test
    }
    cmom3->cd();
    hmom3->Draw();
    hmom3->SetXTitle("p_{3}(GeV/C)");
    hmom3->SetYTitle("count");

    // get mean
    std::cout << " Mom3 mean: " << hmom3->GetMean() << " GeV/c "
              << " Mom3 sigma: " << hmom3->GetStdDev() << " GeV/c" << std::endl;

    ctof->cd();
    legend->Draw();
    cmom->cd();
    legend->Draw();
    // ctof->SaveAs("png/tof_run002.png");
    // cmom->SaveAs("png/mom_run002.png");
    // cmom3->SaveAs("png/realmom_run030.png");
}
