
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <Math/Vector4D.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>
#include <cmath>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

void analysis_triggers_reco() {
    
    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;
    
    // Data
    std::string data_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    TString MC_name = TString::Format("%s_%s", data_name.c_str(), muon_type.c_str());
    TString type = "reco";

    // Cuts
    float pT_JPsi_min = config["cuts_JPsi"]["pT_JPsi_min"];
    float pT_JPsi_max = config["cuts_JPsi"]["pT_JPsi_max"];

    float eta_JPsi_min = config["cuts_JPsi"]["eta_JPsi_min"];
    float eta_JPsi_max = config["cuts_JPsi"]["eta_JPsi_max"];
    
    float pT_mu_min = config["cuts_mu"]["pT_mu_min"];
    float pT_mu_max = config["cuts_mu"]["pT_mu_max"];
    
    float eta_mu_min = config["cuts_mu"]["eta_mu_min"];
    float eta_mu_max = config["cuts_mu"]["eta_mu_max"];
    
    int n_files = config["n_files"];

    TFile* outFile = TFile::Open(TString::Format("results/%s/%s/eventmuons.root", MC_name.Data(), type.Data()), "RECREATE");
    TTree* outTree = new TTree("Triggers", "JPsi triggers");

    std::string category = "All"; 
    double pT, eta, phi, mass;
    std::vector<double> pT_assocs, eta_assocs, phi_assocs;
    std::vector<int> MotherPID;
    int nEvents = 0;

    outTree->Branch("category", &category);
    outTree->Branch("pT",  &pT,  "pT/D");
    outTree->Branch("eta",  &eta,  "eta/D");
    outTree->Branch("phi",  &phi,  "phi/D");
    outTree->Branch("mass",  &mass,  "mass/D");
    outTree->Branch("pT_assocs",  &pT_assocs);
    outTree->Branch("eta_assocs",  &eta_assocs);
    outTree->Branch("phi_assocs",  &phi_assocs);
    outTree->Branch("MotherPID", &MotherPID);

    TTree* metaDataTree = new TTree("MetaData", "Metadata about the analysis");
    metaDataTree->Branch("nEvents", &nEvents, "nEvents/I");

    // int exceptions[] = {5, 7, 8}; // Files with issues (e.g. missing trees)
    int exceptions[] = {3}; // ONLY FOR DATA FOR NOW
    // int exceptions[] = {-9999}; // Files with issues (e.g. missing trees)

    for (int i = 0; i < n_files; ++i) {

        // Skip files with known issues
        if (std::find(std::begin(exceptions), std::end(exceptions), i) != std::end(exceptions)) {
            std::cout << "Skipping file " << i << " due to known issues.\n";
            continue;
        }

        std::cout << "Processing file " << i << " of " << n_files << std::endl;
        TString data_file = TString::Format("results/%s/%s/muonAOD%d.root", MC_name.Data(), type.Data(), i);

        // Load the dataframe keys
        TFile *file = TFile::Open(data_file);
        TIter nextKey(file->GetListOfKeys());
        TKey* key;

        int dirCount = 0;

        // Loop over dataframes
        while ((key = (TKey*) nextKey())) {

            // Load directory and tree
            TObject* obj = key->ReadObj();
            if (!(obj->InheritsFrom("TDirectory"))) continue;
            
            TDirectory* dir = (TDirectory*) obj;
            TTree *tree = (TTree*)dir->Get("O2dqmuontable");

            std::cout << TString::Format("Reading tracks from dir %d of %d: %s\r", dirCount, file->GetListOfKeys()->GetEntries(), dir->GetName()) << std::flush;

            // Group muons by event index
            std::map<ULong64_t, std::vector<Long64_t>> muon_groups;
            ULong64_t fEventIdx;
            tree->SetBranchAddress("fEventIdx", &fEventIdx);
            
            // First pass: build groups of muons from the same event
            Long64_t n = tree->GetEntries();
            for (Long64_t i = 0; i < n; ++i) {
                tree->GetEntry(i);
                muon_groups[fEventIdx].push_back(i);
            }
            nEvents += muon_groups.size(); // Count unique events for metadata

            // Prepare to label
            Long64_t fGlobalIndexAssoc;
            tree->SetBranchAddress("fGlobalIndexassoc", &fGlobalIndexAssoc);

            // Prepare to read muon kinematics
            float fPt, fPhi, fEta;
            tree->SetBranchAddress("fPtassoc", &fPt);
            tree->SetBranchAddress("fPhiassoc", &fPhi);
            tree->SetBranchAddress("fEtaassoc", &fEta);

            // Second pass: process each group
            for (auto& event : muon_groups) {

                ULong64_t eventID = event.first;
                auto& muon_entries = event.second;

                if (muon_entries.size() < 2) continue; // Needs at least 2 muons to form a pair

                // std::vector<Long64_t> muon_motherIDs, muon_motherPDGs;
                std::vector<ROOT::Math::PtEtaPhiMVector> muon_vectors;

                // Read muon kinematics and build 4-vectors
                for (auto entry : muon_entries) {
                    tree->GetEntry(entry);
                    ROOT::Math::PtEtaPhiMVector muon_vec(fPt, fEta, fPhi, 0.105658); // Muon mass ~105.658 MeV/c^2
                    muon_vectors.push_back(muon_vec);
                }

                // Store the indexes of the best candidate muon pair
                ROOT::Math::PtEtaPhiMVector JPsiCandidate(-9999,0,0,-1); 
                size_t idx_cand_1, idx_cand_2;
                
                // Find the muon pair with invariant mass closest to J/Psi mass
                for (size_t j = 0; j < muon_vectors.size(); ++j) {
                    if (muon_vectors[j].Pt() < pT_mu_min || muon_vectors[j].Pt() > pT_mu_max) continue;
                    if (muon_vectors[j].Eta() < eta_mu_min || muon_vectors[j].Eta() > eta_mu_max) continue;

                    for (size_t k = j + 1; k < muon_vectors.size(); ++k) {
                        if (muon_vectors[k].Pt() < pT_mu_min || muon_vectors[k].Pt() > pT_mu_max) continue;
                        if (muon_vectors[k].Eta() < eta_mu_min || muon_vectors[k].Eta() > eta_mu_max) continue;

                        auto track = muon_vectors[j] + muon_vectors[k];

                        pT_assocs.clear();
                        eta_assocs.clear();
                        phi_assocs.clear();
                        MotherPID.clear();
 
                        // Replace candidate if closer to J/Psi mass
                        // if (std::abs(track.M() - 3.0969) < std::abs(JPsiCandidate.M() - 3.0969)) { // J/Psi mass ~3.0969 GeV/c^2
                        //     JPsiCandidate = track;
                        //     idx_cand_1 = j;
                        //     idx_cand_2 = k;
                        // }
                        JPsiCandidate = track;
                        idx_cand_1 = j;
                        idx_cand_2 = k;
                        if (JPsiCandidate.Eta() < eta_JPsi_min || JPsiCandidate.Eta() > eta_JPsi_max) continue; // Apply J/Psi eta cut
                        if (JPsiCandidate.Pt() < pT_JPsi_min || JPsiCandidate.Pt() > pT_JPsi_max) continue; // Apply J/Psi pT cut

                        pT = JPsiCandidate.Pt();
                        eta = JPsiCandidate.Eta();
                        phi = JPsiCandidate.Phi();
                        mass = JPsiCandidate.M();

                        // Loop over muons again to find assocs
                        for (size_t a = 0; a < muon_vectors.size(); ++a) {
                            if (a == idx_cand_1 || a == idx_cand_2) continue; // Skip the trigger muons
                            auto assocTrack = muon_vectors[a];
                            pT_assocs.push_back(assocTrack.Pt());
                            eta_assocs.push_back(assocTrack.Eta());
                            phi_assocs.push_back(assocTrack.Phi());
                        }

                        outTree->Fill();
                    }
                }
            }

            tree->ResetBranchAddresses();
            delete tree;

            dirCount++;
        }

        std::cout << std::endl;

        delete key;
        file->Close();
        delete file;
    }

    // Fill metadata tree
    metaDataTree->Fill();

    outFile->cd();
    outTree->Write();
    metaDataTree->Write();
    outFile->Close();
}