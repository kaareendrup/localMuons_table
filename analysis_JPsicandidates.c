
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
#include <map>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "utils/hists.c"

void analysis_JPsiHists(TString type) {
    
    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;
    
    // Data
    std::string data_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    TString MC_name = TString::Format("%s_%s", data_name.c_str(), muon_type.c_str());

    // Cuts
    float eta_JPsi_min = config["cuts_JPsi"]["eta_JPsi_min"];
    float eta_JPsi_max = config["cuts_JPsi"]["eta_JPsi_max"];
    
    int n_files = config["n_files"];

    // Histogram parameters
    bool scale_by_y = config["scale_by_y"];
    int n_bins_mass = config["hists"]["n_bins_mass"];
    float signal_range_min = config["signal_range"]["min"];
    float signal_range_max = config["signal_range"]["max"];
    float background_range_min = config["background_range"]["min"];
    float background_range_max = config["background_range"]["max"];
    
    std::vector<double> pT_bins = config["hists"]["pT_bins"].get<std::vector<double>>();

    // Set outfile
    TFile* outFile = TFile::Open(TString::Format("results/%s/%s/invMassSpektra.root", MC_name.Data(), type.Data()), "RECREATE");
    
    // Load events
    TString data_file = TString::Format("results/%s/%s/eventmuons.root", MC_name.Data(), type.Data());
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
    
    double pT, eta, phi, mass;
    std::string *category = nullptr;

    tree->SetBranchAddress("pT",  &pT);
    tree->SetBranchAddress("eta",  &eta);
    tree->SetBranchAddress("phi",  &phi);
    tree->SetBranchAddress("mass",  &mass);
    tree->SetBranchAddress("category", &category);
    
    // Create histograms
    std::map<int, std::unique_ptr<TH1F>> invMassHists;
    TH1D* pT_sig = new TH1D("pT_sig", "pT of signal region;p_{T} GeV/c;Counts", pT_bins.size() - 1, pT_bins.data());
    TH1D* pT_bkg_low = new TH1D("pT_bkg_low", "pT of lower background region;p_{T} GeV/c;Counts", pT_bins.size() - 1, pT_bins.data());
    TH1D* pT_bkg_high = new TH1D("pT_bkg_high", "pT of higher background region;p_{T} GeV/c;Counts", pT_bins.size() - 1, pT_bins.data());

    // Loop over entries and create the necessary histograms
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {

        tree->GetEntry(i);
        std::cout << "Processing entry " << i+1 << " of " << tree->GetEntries() << "\r" << std::flush;
        if (type == "gen" && *category != "JPsi") continue; // Only consider J/Psi category for gen-level analysis

        TString category_str;
        for (size_t j = 0; j < pT_bins.size()-1; ++j) {
            if (pT >= pT_bins[j] && pT < pT_bins[j+1]) {
                category_str = TString::Format("pT_%.1f_%.1f", pT_bins[j], pT_bins[j+1]);
                fillHist(j, category_str + "_invMass", invMassHists, mass, n_bins_mass, 1.0, 5.0);   
                break;
            }
        }

        double w = 1.0; 
        if (scale_by_y) {
            w = 1.0 / getDeltaY(pT, eta_JPsi_min, eta_JPsi_max);
        }
        if (mass >= signal_range_min && mass <= signal_range_max) {
            pT_sig->Fill(pT, w);
        } else if (mass >= background_range_min && mass < signal_range_min) {
            pT_bkg_low->Fill(pT, w);
        } else if (mass > signal_range_max && mass <= background_range_max) {
            pT_bkg_high->Fill(pT, w);
        }
    }
    std::cout << std::endl;

    // Write data
    outFile->cd();
    for (auto& mHist : invMassHists) {
        mHist.second->Write();
    }
    pT_sig->Write();
    pT_bkg_low->Write();
    pT_bkg_high->Write();
    
    // Propagate metadata
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

void analysis_JPsicandidates() {

    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;
    std::string data_name = config["data_name"];

    analysis_JPsiHists("reco");
    if (!(data_name == "DQ_data")) {
        analysis_JPsiHists("gen");
    }
}