
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TCanvas.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

void fillHist(const int& idx, std::map<int, std::unique_ptr<TH1F>>& hists, double value, int n_bins, double x_min, double x_max, double weight = 1.0) {

    // If histogram doesn't exist yet, create it
    if (hists.find(idx) == hists.end()) {
        hists[idx] = std::make_unique<TH1F>(TString::Format("hist_%d", idx), "Invariant Mass", n_bins, x_min, x_max);
        hists[idx]->SetDirectory(nullptr);
    }

    hists[idx]->Fill(value, weight);
}

void inspect_invMass() {

    ////////////////////////////////////////////////////////////////////
    ////            Load configuration, setup up filenames          ////
    ////////////////////////////////////////////////////////////////////
    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    // Data
    std::string dataset_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    std::string eff_source = config["data_eff_source"];

    TString data_name = TString::Format("%s_%s", dataset_name.c_str(), muon_type.c_str());
    TString MC_name = TString::Format("%s_%s", eff_source.c_str(), muon_type.c_str());

    TString type = "reco";

    bool is_MC = !(dataset_name == "DQ_data");

    // Cuts
    float pT_trigger_min = config["cuts_JPsi"]["pT_JPsi_min"];
    float pT_trigger_max = config["cuts_JPsi"]["pT_JPsi_max"];
    float eta_trigger_min = config["cuts_JPsi"]["eta_JPsi_min"];
    float eta_trigger_max = config["cuts_JPsi"]["eta_JPsi_max"];
    
    float pT_assoc_min = config["cuts_mu"]["pT_mu_min"];
    float pT_assoc_max = config["cuts_mu"]["pT_mu_max"];
    float eta_assoc_min = config["cuts_mu"]["eta_mu_min"];
    float eta_assoc_max = config["cuts_mu"]["eta_mu_max"];

    // Histogram parameters
    int n_bins_mass = config["hists"]["n_bins_mass"];
    int n_bins_pT = config["hists"]["n_bins_pT"];
    float deltaEta_min = config["hists"]["deltaEta_min"];
    float deltaEta_max = config["hists"]["deltaEta_max"];

    float signal_range_min = config["signal_range"]["min"];
    float signal_range_max = config["signal_range"]["max"];
    float background_range_min = config["background_range"]["min"];
    float background_range_max = config["background_range"]["max"];

    ////////////////////////////////////////////////////////////////////
    ////            Open input file tree, setup pT histograms       ////
    ////////////////////////////////////////////////////////////////////
    TString data_file = TString::Format("results/%s/%s/eventmuons_nocuts.root", data_name.Data(), type.Data());
    TFile *file = TFile::Open(data_file);
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file\n";
        return;
    }

    TTree *tree = nullptr;
    file->GetObject("Triggers", tree);
    if (!tree) {
        std::cerr << "Tree 'Triggers' not found\n";
        return;
    }

    std::string *category = nullptr;
    double pT, eta, phi, mass;
    std::vector<double> *pT_assocs = nullptr;
    std::vector<double> *eta_assocs = nullptr;
    std::vector<double> *phi_assocs = nullptr;
    std::vector<int> *MotherPID = nullptr;

    tree->SetBranchAddress("category", &category);
    tree->SetBranchAddress("pT",  &pT);
    tree->SetBranchAddress("eta",  &eta);
    tree->SetBranchAddress("phi",  &phi);
    tree->SetBranchAddress("mass",  &mass);
    tree->SetBranchAddress("pT_assocs",  &pT_assocs);
    tree->SetBranchAddress("eta_assocs",  &eta_assocs);
    tree->SetBranchAddress("phi_assocs",  &phi_assocs);
    tree->SetBranchAddress("MotherPID", &MotherPID);

    TH1F *invMassBeforeCuts = new TH1F("invMassBeforeCuts", "Invariant Mass Before Cuts", n_bins_mass, 1.0, 5.0);
    TH1F *invMassIntermediate = new TH1F("invMassIntermediate", "Invariant Mass Intermediate", n_bins_mass, 1.0, 5.0);
    TH1F *invMassAfterCuts = new TH1F("invMassAfterCuts", "Invariant Mass After Cuts", n_bins_mass, 1.0, 5.0);

    TH2F *invMassPtBeforeCuts = new TH2F("invMassPtBeforeCuts", "Inv mass vs assoc pT Before Cuts", n_bins_mass, 1.0, 5.0, n_bins_pT, 0.0, 20.0);
    TH2F *invMassPtIntermediate = new TH2F("invMassPtIntermediate", "Inv mass vs assoc pT Intermediate", n_bins_mass, 1.0, 5.0, n_bins_pT, 0.0, 20.0);
    TH2F *invMassPtAfterCuts = new TH2F("invMassPtAfterCuts", "Inv mass vs assoc pT After Cuts", n_bins_mass, 1.0, 5.0, n_bins_pT, 0.0, 20.0);

    std::map<int, std::unique_ptr<TH1F>> invMassHistsPairsPtBeforeCuts;
    std::map<int, std::unique_ptr<TH1F>> invMassHistsPairsPtIntermediate;
    std::map<int, std::unique_ptr<TH1F>> invMassHistsPairsPtAfterCuts;

    ////////////////////////////////////////////////////////////////////////
    ////    Loop over J/Psi candidate entries and fill histograms       ////
    ////////////////////////////////////////////////////////////////////////

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {

        tree->GetEntry(i);
        std::cout << "Processing entry " << i+1 << " of " << tree->GetEntries() << "\r" << std::flush;
        if (*category == "All_1") continue; // Skip same-sign pairs for this test

        // Calculate delta eta and phi for associated particles
        for (size_t j = 0; j < eta_assocs->size(); ++j) {
            int bin = static_cast<int>(pT_assocs->at(j) * n_bins_pT / 20.0);

            invMassBeforeCuts->Fill(mass);
            invMassPtBeforeCuts->Fill(mass, pT_assocs->at(j));
            fillHist(bin, invMassHistsPairsPtBeforeCuts, mass, n_bins_mass, 1.0, 5.0);
            
            // First cuts
            if (pT < pT_trigger_min || pT > pT_trigger_max) continue;
            if (eta < eta_trigger_min || eta > eta_trigger_max) continue;
            if (eta_assocs->at(j) < eta_assoc_min || eta_assocs->at(j) > eta_assoc_max) continue;
            
            invMassIntermediate->Fill(mass);
            invMassPtIntermediate->Fill(mass, pT_assocs->at(j));
            fillHist(bin, invMassHistsPairsPtIntermediate, mass, n_bins_mass, 1.0, 5.0);
            
            // pT assoc cut
            if (pT_assocs->at(j) < pT_assoc_min || pT_assocs->at(j) > pT_assoc_max) continue;
            
            invMassAfterCuts->Fill(mass);
            invMassPtAfterCuts->Fill(mass, pT_assocs->at(j));
            fillHist(bin, invMassHistsPairsPtAfterCuts, mass, n_bins_mass, 1.0, 5.0);
        }
    }
    std::cout << std::endl;

    TCanvas *c1 = new TCanvas("c1", "Invariant Mass Distributions", 1200, 400);
    c1->Divide(3, 1);
    c1->cd(1);
    invMassBeforeCuts->Draw();
    c1->cd(2);
    invMassIntermediate->Draw();
    c1->cd(3);
    invMassAfterCuts->Draw();
    c1->SaveAs("invMassTest.png");
    
    TCanvas *c2 = new TCanvas("c2", "Invariant Mass vs pT of Associate", 1200, 400);
    c2->Divide(3, 1);
    c2->cd(1);
    invMassPtBeforeCuts->Draw("COLZ");
    gPad->SetLogz();
    c2->cd(2);
    invMassPtIntermediate->Draw("COLZ");
    gPad->SetLogz();
    c2->cd(3);
    invMassPtAfterCuts->Draw("COLZ");
    gPad->SetLogz();
    c2->SaveAs("invMassPtTest.png");

    TCanvas *c3 = new TCanvas("c3", "Invariant Mass binned by Associate pT", 1200, 400);
    c3->Divide(3, 1);
    
    c3->cd(1);
    if (invMassHistsPairsPtBeforeCuts.find(0) != invMassHistsPairsPtBeforeCuts.end()) {
        invMassHistsPairsPtBeforeCuts[0]->Draw();
    }
    gPad->SetLogy();
    c3->cd(2);
    if (invMassHistsPairsPtIntermediate.find(0) != invMassHistsPairsPtIntermediate.end()) {
        invMassHistsPairsPtIntermediate[0]->Draw();
    }
    gPad->SetLogy();
    c3->cd(3);
    if (invMassHistsPairsPtAfterCuts.find(0) != invMassHistsPairsPtAfterCuts.end()) {
        invMassHistsPairsPtAfterCuts[0]->Draw();
    }
    gPad->SetLogy();
    for (int b = 1; b < n_bins_pT; ++b) {
        c3->cd(1);
        if (invMassHistsPairsPtBeforeCuts.find(b) != invMassHistsPairsPtBeforeCuts.end()) {
            invMassHistsPairsPtBeforeCuts[b]->Draw("same");
        }
        c3->cd(2);
        if (invMassHistsPairsPtIntermediate.find(b) != invMassHistsPairsPtIntermediate.end()) {
            invMassHistsPairsPtIntermediate[b]->Draw("same");
        }
        c3->cd(3);
        if (invMassHistsPairsPtAfterCuts.find(b) != invMassHistsPairsPtAfterCuts.end()) {
            invMassHistsPairsPtAfterCuts[b]->Draw("same");
        }
    }
    c3->SaveAs("invMassPtTestBinned.png");
}