
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <TH1F.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "utils/hists.c"

void analysis_correlations() {

    ////////////////////////////////////////////////////////////////////
    ////            Load configuration, setup up filenames          ////
    ////////////////////////////////////////////////////////////////////
    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    // Data
    std::string dataset_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    TString data_name = TString::Format("%s_%s", dataset_name.c_str(), muon_type.c_str());
    
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
    int n_bins = config["hists"]["n_bins_correlations"];
    float deltaEta_min = config["hists"]["deltaEta_min"];
    float deltaEta_max = config["hists"]["deltaEta_max"];

    float signal_range_min = config["signal_range"]["min"];
    float signal_range_max = config["signal_range"]["max"];
    float background_range_min = config["background_range"]["min"];
    float background_range_max = config["background_range"]["max"];

    std::vector<double> pT_bins = config["hists"]["pT_bins"].get<std::vector<double>>();

    // Set outfile
    TFile* outFile = TFile::Open(TString::Format("results/%s/%s/analysis.root", data_name.Data(), type.Data()), "RECREATE");
    TTree* triggerCounts = new TTree("Correlations", "JPsi correlations");
    
    // Define branches for output tree
    std::string category_out; 
    int count_out;

    triggerCounts->Branch("category", &category_out);
    triggerCounts->Branch("count", &count_out, "count/I");

    ////////////////////////////////////////////////////////////////////
    ////            Open input file tree, setup pT histograms       ////
    ////////////////////////////////////////////////////////////////////
    TString data_file = TString::Format("results/%s/%s/eventmuons.root", data_name.Data(), type.Data());
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

    std::map<TString, int> trigger_counts;
    std::map<TString, std::unique_ptr<TH1F>> invMassHists;
    std::map<TString, std::unique_ptr<TH1F>> deltaEtaHists;
    std::map<TString, std::unique_ptr<TH1F>> deltaPhiHists;
    std::map<TString, std::unique_ptr<TH1F>> deltaPhiHistspT;

    ////////////////////////////////////////////////////////////////////////
    ////    Loop over J/Psi candidate entries and fill histograms       ////
    ////////////////////////////////////////////////////////////////////////

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {

        tree->GetEntry(i);
        std::cout << "Processing entry " << i+1 << " of " << tree->GetEntries() << "\r" << std::flush;

        TString category_str = TString(*category);
        trigger_counts[category_str]++;
        fillHist(category_str + "_invMass", invMassHists, mass, n_bins_mass, 1.0, 5.0);

        // Calculate delta eta and phi for associated particles
        for (size_t j = 0; j < eta_assocs->size(); ++j) {

            double deltaEta = eta - eta_assocs->at(j);
            double deltaPhi = phi - phi_assocs->at(j);
            if (deltaPhi > 3.0 / 2.0 * M_PI) {
                deltaPhi -= 2.0 * M_PI;
            }
            if (deltaPhi < -0.5 * M_PI) {
                deltaPhi += 2.0 * M_PI;
            }

            TString full_category = category_str;

            // pT cuts
            if (pT < pT_trigger_min || pT > pT_trigger_max) continue;
            if (pT_assocs->at(j) < pT_assoc_min || pT_assocs->at(j) > pT_assoc_max) continue;

            // Eta cuts
            if (eta < eta_trigger_min || eta > eta_trigger_max) continue;
            if (eta_assocs->at(j) < eta_assoc_min || eta_assocs->at(j) > eta_assoc_max) continue;

            if (is_MC) {
                if ((std::abs(MotherPID->at(j)) >= 411 && std::abs(MotherPID->at(j)) <= 445) || 
                    (std::abs(MotherPID->at(j)) >= 4101 && std::abs(MotherPID->at(j)) <= 4444)) { // J/psi from charm
                    full_category += "_charm";
                } else { // J/psi not from charm
                    full_category += "_noncharm";
                }
            } else {
                if (mass > signal_range_min && mass < signal_range_max) {
                    full_category += "_signal";
                } else if (mass > background_range_min && mass < signal_range_min || 
                           mass > signal_range_max && mass < background_range_max) {
                    full_category += "_background";
                } else {
                    continue; // Skip if not in signal or background range
                }
            }

            fillHist(full_category + "_deltaEta", deltaEtaHists, deltaEta, n_bins, deltaEta_min, deltaEta_max);
            fillHist(full_category + "_deltaPhi", deltaPhiHists, deltaPhi, n_bins, -0.5 * M_PI, 3.0 / 2.0 * M_PI);

            for (int p = 0; p < 10; ++p) {
                float ptmin = pT_bins[p];
                float ptmax = pT_bins[p + 1];
                if (pT_assocs->at(j) >= ptmin && pT_assocs->at(j) < ptmax) {
                    fillHist(TString::Format("%s_deltaPhi_pT_%.1f_%.1f", full_category.Data(), ptmin, ptmax), deltaPhiHistspT, deltaPhi, n_bins, -0.5 * M_PI, 3.0 / 2.0 * M_PI);
                }
            }
        }
    }
    std::cout << std::endl;

    outFile->cd();
    for (auto& mHist : invMassHists) {
        mHist.second->Write();
    }
    for (auto& dEtaHist : deltaEtaHists) {
        dEtaHist.second->Write();
    }
    for (auto& dPhiHist : deltaPhiHists) {
        dPhiHist.second->Write();
    }
    for (auto& dPhiHist : deltaPhiHistspT) {
        dPhiHist.second->Write();
    }

    for (auto& count : trigger_counts) {
        category_out = count.first;
        count_out = count.second;
        triggerCounts->Fill();
    }

    ////////////////////////////////////////////////////////////////////////////
    ////    Load metadata tree with event count information and propagate   ////
    ////////////////////////////////////////////////////////////////////////////
    TTree *metaTree = nullptr;
    file->GetObject("MetaData", metaTree);
    int *nEvents = nullptr;
    metaTree->SetBranchAddress("nEvents", &nEvents);
    metaTree->GetEntry(0);
    
    int *nEventsOut = nullptr;
    TTree* metaDataTree = new TTree("MetaData", "Metadata about the analysis");
    metaDataTree->Branch("nEvents", &nEventsOut, "nEvents/I");
    nEventsOut = nEvents; // Propagate the number of events to the new tree
    std::cout << "Number of events: " << nEventsOut << std::endl;
    metaDataTree->Fill(); // Fill with the propagated number of events

    outFile->Write();
    outFile->Close();
}