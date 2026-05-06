
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

void make_fit(const nlohmann::json& config, TString data_name, float xmin, float xmax) {

    string muon_type = config["muon_type"];
    std::vector<double> pT_bins = config["hists"]["pT_bins"].get<std::vector<double>>();
    float signal_range_min = config["signal_range"]["min"];
    float signal_range_max = config["signal_range"]["max"];
    float background_range_min = config["background_range"]["min"];
    float background_range_max = config["background_range"]["max"];
    
    // Load data
    TString in_file = TString::Format("results/%s_%s/reco/invMassSpektra.root", data_name.Data(), muon_type.c_str());
    TFile* file = TFile::Open(in_file, "READ");
    
    TCanvas *c1 = new TCanvas("c1", "Invariant Mass Spectra by pT", 700, 500);
    
    TString category_str = TString::Format("pT_%.1f_%.1f_invMass", pT_bins[0], pT_bins[1]);
    TH1F *invMass_pt = (TH1F*)file->Get(category_str);
    invMass_pt->SetDirectory(0); // Detach histogram from file

    invMass_pt->SetLineWidth(2);
    invMass_pt->SetLineColor(kRed);
    invMass_pt->SetTitle(category_str);
    invMass_pt->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
    invMass_pt->GetYaxis()->SetTitle("Counts");
    // invMass_pt->Sumw2();
    invMass_pt->Draw();

    double yMinSingle = invMass_pt->GetMinimum();
    double yMaxSingle = invMass_pt->GetMaximum();

    invMass_pt->Sumw2();

    // Draw shaded boxes for background regions
    TBox* box1s = new TBox(background_range_min, yMinSingle, signal_range_min, yMaxSingle);
    box1s->SetFillColorAlpha(kGray, 1); // semi-transparent
    box1s->SetFillStyle(3004);
    box1s->Draw("same");

    TBox* box2s = new TBox(signal_range_max, yMinSingle, background_range_max, yMaxSingle);
    box2s->SetFillColorAlpha(kGray, 1); // semi-transparent
    box2s->SetFillStyle(3004);
    box2s->Draw("same");

    int n_signal = invMass_pt->Integral(invMass_pt->FindBin(signal_range_min), invMass_pt->FindBin(signal_range_max));
    int n_bkg_low = invMass_pt->Integral(invMass_pt->FindBin(background_range_min), invMass_pt->FindBin(signal_range_min));
    int n_bkg_high = invMass_pt->Integral(invMass_pt->FindBin(signal_range_max), invMass_pt->FindBin(background_range_max));
    int subtracted_signal = n_signal - (n_bkg_low + n_bkg_high);

        TF1 *model = nullptr;
        if (data_name == "c3") {
            model = new TF1(
                "model",
                "crystalball(0) + pol2(5)",
                xmin, xmax
            );

            model->SetParameters(
                3000,   // amplitude
                3.1,  // mean
                .2,  // sigma
                0.6,  // alpha
                2.0,  // n
                10,    // p0 background
                -1,   // p1 background
                0       // p2 background
            );

        } else {
            model = new TF1(
                "model",
                "crystalball(0) + expo(5)",
                xmin, xmax
            );
        
            model->SetParameters(
                8000,   // amplitude
                3.1,  // mean
                .1,  // sigma
                0.6,  // alpha
                2.0,  // n
                20.,    // p0 background
                -3.  // p1 background
            );
            // model->FixParameter(0, 10000);
            model->FixParameter(1, 3.1);
            model->FixParameter(2, 0.088);
            model->FixParameter(3, 0.6); 
            model->FixParameter(4, 2.0); 
        }

        invMass_pt->Fit(model, "R");

        TF1 *bkg = nullptr;
        if (data_name == "c3") {
            bkg = new TF1("bkg", "pol2", xmin, xmax);
            bkg->SetParameters(
                model->GetParameter(5),
                model->GetParameter(6),
                model->GetParameter(7)
            );
        } else {
            bkg = new TF1("bkg", "expo", xmin, xmax);
            bkg->SetParameters(
                model->GetParameter(5),
                model->GetParameter(6)
            );
        }

        TH1 *h_sub = (TH1*)invMass_pt->Clone("h_sub");
        h_sub->Add(bkg, -1);

        TF1 *signal = new TF1("signal", "crystalball", xmin, xmax);
        signal->SetParameters(
            model->GetParameter(0),
            model->GetParameter(1),
            model->GetParameter(2),
            model->GetParameter(3),
            model->GetParameter(4)
        );

        double yield = signal->Integral(xmin, xmax) / invMass_pt->GetBinWidth(1);

        model->SetLineColor(kBlack);
        model->Draw("same");
        signal->SetLineColor(kBlack);
        signal->SetLineStyle(kDashed);
        signal->Draw("same");
        bkg->SetLineStyle(kDashed);
        bkg->SetLineColor(kBlack);
        bkg->Draw("same");

        std::cout << data_name << std::endl;
        std::cout << "Yield from fit: " << yield << std::endl;
        std::cout << "Subtracted signal counts: " << subtracted_signal << std::endl << std::endl;

    c1->SaveAs(TString::Format("results/%s_%s/fit_result.png", data_name.Data(), muon_type.c_str()));

}

void analysis_fit_test() {

    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;
    
    TString data_name = "DQ_data";
    TString MC_name = "c3";

    make_fit(config, MC_name, 1.0, 5.0);
    make_fit(config, data_name, 2.0, 5.0);
}