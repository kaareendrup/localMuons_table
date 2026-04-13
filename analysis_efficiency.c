
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
#include <set>
#include <map>
#include <cmath>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

bool charm_beauty_cut(Long64_t motherPDG) {
    return ((std::abs(motherPDG) == 443) || (std::abs(motherPDG) == 100443) ||
            (std::abs(motherPDG) >= 411 && std::abs(motherPDG) <= 445) || 
            (std::abs(motherPDG) >= 4101 && std::abs(motherPDG) <= 4444) || 
            (std::abs(motherPDG) >= 511 && std::abs(motherPDG) <= 557) || 
            (std::abs(motherPDG) >= 5101 && std::abs(motherPDG) <= 5554));
}

void analysis_efficiency() {

    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    // Data
    std::string data_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    TString MC_name = TString::Format("%s_%s", data_name.c_str(), muon_type.c_str());

    // Cuts
    float pT_JPsi_min = config["cuts_JPsi"]["pT_JPsi_min"];
    float pT_JPsi_max = config["cuts_JPsi"]["pT_JPsi_max"];

    float eta_JPsi_min = config["cuts_JPsi"]["eta_JPsi_min"];
    float eta_JPsi_max = config["cuts_JPsi"]["eta_JPsi_max"];
    
    float pT_mu_min = config["cuts_mu"]["pT_mu_min"];
    float pT_mu_max = config["cuts_mu"]["pT_mu_max"];
    
    float eta_mu_min = config["cuts_mu"]["eta_mu_min"];
    float eta_mu_max = config["cuts_mu"]["eta_mu_max"];
    
    float signal_range_min = config["signal_range"]["min"];
    float signal_range_max = config["signal_range"]["max"];

    int n_files = config["n_files"];
    
    // Setup output file and trees
    TFile* outFile = TFile::Open(TString::Format("results/%s/particles.root", MC_name.Data()), "RECREATE");
    TTree* outTreeMuonsReco = new TTree("MuonsReco", "Reconstructed Muons");
    TTree* outTreeJPsiReco = new TTree("JPsiReco", "Reconstructed JPsi");
    TTree* outTreeMuonsGen = new TTree("MuonsGen", "Generated Muons");
    TTree* outTreeJPsiGen = new TTree("JPsiGen", "Generated JPsi");
    TTree* metaData = new TTree("MetaData", "Event selection metadata");

    double pTJPsiReco, pTJPsiGen, pTJPsiReco_true, pTMuonReco, pTMuonGen, pTMuonReco_true;
    double massJPsiReco, massJPsiGen;
    double etaJPsiReco, etaJPsiGen, etaMuonReco, etaMuonGen; 

    outTreeJPsiReco->Branch("pTJPsi",  &pTJPsiReco,  "pT/D");
    outTreeJPsiReco->Branch("etaJPsi",  &etaJPsiReco,  "eta/D");
    outTreeJPsiReco->Branch("pTJPsi_true",  &pTJPsiReco_true,  "pT/D"); 
    outTreeJPsiReco->Branch("massJPsi",  &massJPsiReco,  "mass/D");
    outTreeJPsiGen->Branch("pTJPsi",  &pTJPsiGen,  "pT/D");
    outTreeJPsiGen->Branch("etaJPsi",  &etaJPsiGen,  "eta/D");
    outTreeJPsiGen->Branch("massJPsi",  &massJPsiGen,  "mass/D");
    outTreeMuonsReco->Branch("pTMuon",  &pTMuonReco,  "pT/D");
    outTreeMuonsReco->Branch("etaMuon",  &etaMuonReco,  "eta/D");
    outTreeMuonsReco->Branch("pTMuon_true",  &pTMuonReco_true,  "pT/D");
    // outTreeMuonsReco->Branch("etaMuon_true",  &etaMuonReco_true,  "eta/D");
    outTreeMuonsGen->Branch("pTMuon",  &pTMuonGen,  "pT/D");
    outTreeMuonsGen->Branch("etaMuon",  &etaMuonGen,  "eta/D");

    // Metadata branches
    std::vector<float> pTCuts = {pT_JPsi_min, pT_JPsi_max, pT_mu_min, pT_mu_max};
    std::vector<float> etaCuts = {eta_JPsi_min, eta_JPsi_max, eta_mu_min, eta_mu_max};
    metaData->Branch("pTCuts", &pTCuts);
    metaData->Branch("etaCuts", &etaCuts);
    metaData->Fill();
    
    // Loop over gen and reco files
    for (int i = 0; i < n_files; ++i) {

        std::cout << "Processing file " << i << " of " << n_files << std::endl;
        TString reco_file = TString::Format("results/%s/reco/muonAOD%d.root", MC_name.Data(), i);
        TString gen_file = TString::Format("results/%s/gen/muonAOD%d.root", MC_name.Data(), i);
    
        TFile *recoFile = TFile::Open(reco_file);
        TFile *genFile = TFile::Open(gen_file);

        // Load the dataframe keys
        TIter nextGenKey(genFile->GetListOfKeys());
        TKey* genKey;

        TIter nextRecoKey(recoFile->GetListOfKeys());
        TKey* recoKey;
        
        int dirCount = 0;
        int JPsiCount = 0;
        int JPsiCount_mu = 0;

        // Loop over GEN dataframes
        while ((genKey = (TKey*) nextGenKey())) {

            // Load directory and tree
            TObject* obj = genKey->ReadObj();
            if (!(obj->InheritsFrom("TDirectory"))) continue;

            TDirectory* genDir = (TDirectory*) obj;
            TTree *jpsiGenTree = (TTree*)genDir->Get("O2dqjpsitable");
            TTree *muonGenTree = (TTree*)genDir->Get("O2dqmuontable");

            std::cout << TString::Format("Reading tracks from dir %d of %d: %s\r", dirCount, recoFile->GetListOfKeys()->GetEntries(), genDir->GetName()) << std::flush;
            
            float fPtJPsiGen, fPtMuonGen, fEtaJPsiGen, fEtaMuonGen;
            Long64_t fMotherPDG, fGrandmotherPDG, fTrackPDGJPsiGen, fGlobalIndexMCtrack, fMotherID;
            std::set<Long64_t> MCMuonMothers;
            
            muonGenTree->SetBranchAddress("fPtassoc", &fPtMuonGen);
            muonGenTree->SetBranchAddress("fEtaassoc", &fEtaMuonGen);
            muonGenTree->SetBranchAddress("fMotherPDG", &fMotherPDG);
            muonGenTree->SetBranchAddress("fGrandmotherPDG", &fGrandmotherPDG);
            muonGenTree->SetBranchAddress("fMotherID", &fMotherID);
            for (Long64_t i = 0; i < muonGenTree->GetEntries(); ++i) {
                muonGenTree->GetEntry(i);
                
                if (!(charm_beauty_cut(fMotherPDG) || charm_beauty_cut(fGrandmotherPDG))) continue;
                MCMuonMothers.insert(fMotherID);
                if (fEtaMuonGen < eta_mu_min || fEtaMuonGen > eta_mu_max) continue; // Apply eta cut on muons
                if (fPtMuonGen < pT_mu_min || fPtMuonGen > pT_mu_max) continue; // Apply pT cut on muons
                pTMuonGen = fPtMuonGen;
                etaMuonGen = fEtaMuonGen;
                outTreeMuonsGen->Fill();
            }
            
            jpsiGenTree->SetBranchAddress("fPtassoc", &fPtJPsiGen);
            jpsiGenTree->SetBranchAddress("fEtaassoc", &fEtaJPsiGen);
            jpsiGenTree->SetBranchAddress("fTrackPDG", &fTrackPDGJPsiGen);
            jpsiGenTree->SetBranchAddress("fGlobalIndexMCtrack", &fGlobalIndexMCtrack);
            for (Long64_t i = 0; i < jpsiGenTree->GetEntries(); ++i) {
                jpsiGenTree->GetEntry(i);
                if (std::abs(fTrackPDGJPsiGen) != 443) continue; // Ensure we're looking at J/Psi candidates for now
                if (fEtaJPsiGen < eta_JPsi_min || fEtaJPsiGen > eta_JPsi_max) continue; // Apply eta cut on J/Psi
                if (fPtJPsiGen < pT_JPsi_min || fPtJPsiGen > pT_JPsi_max) continue; // Apply pT cut on J/Psi
                JPsiCount++;
                if (!MCMuonMothers.count(fGlobalIndexMCtrack)) continue; 
                pTJPsiGen = fPtJPsiGen;
                etaJPsiGen = fEtaJPsiGen;
                outTreeJPsiGen->Fill();
                JPsiCount_mu++;
            }

            dirCount++;
        }
        std::cout << std::endl << TString::Format("%d J/Psis passed the eta and pT cuts, %d with muon daughters, which is %f%%", JPsiCount, JPsiCount_mu, (double)JPsiCount_mu / JPsiCount * 100) << std::endl;

        std::cout << std::endl;
        dirCount = 0;

        // Loop over RECO dataframes
        while ((recoKey = (TKey*) nextRecoKey())) {

            // Load directory and tree
            TObject* obj = recoKey->ReadObj();
            if (!(obj->InheritsFrom("TDirectory"))) continue;

            TDirectory* recoDir = (TDirectory*) obj;
            TTree *muonRecoTree = (TTree*)recoDir->Get("O2dqmuontable");

            std::cout << TString::Format("Reading tracks from dir %d of %d: %s\r", dirCount, recoFile->GetListOfKeys()->GetEntries(), recoDir->GetName()) << std::flush;

            // Group muons by event index
            std::map<ULong64_t, std::vector<Long64_t>> muon_groups;
            ULong64_t fEventIdx;
            Long64_t fMotherPDG, fMotherID, fGrandmotherPDG, fGlobalIndexMCtrack;
            muonRecoTree->SetBranchAddress("fEventIdx", &fEventIdx);
            muonRecoTree->SetBranchAddress("fMotherPDG", &fMotherPDG);
            muonRecoTree->SetBranchAddress("fMotherID", &fMotherID);
            muonRecoTree->SetBranchAddress("fGrandmotherPDG", &fGrandmotherPDG);
            muonRecoTree->SetBranchAddress("fGlobalIndexMCtrack", &fGlobalIndexMCtrack);

            // Prepare to read muon kinematics
            float fPt, fPt_true, fPt_mother, fPhi, fPhi_mother, fEta, fEta_mother;
            muonRecoTree->SetBranchAddress("fPtassoc", &fPt);
            muonRecoTree->SetBranchAddress("fPtassoctrue", &fPt_true);
            muonRecoTree->SetBranchAddress("fPtmother", &fPt_mother);
            muonRecoTree->SetBranchAddress("fPhiassoc", &fPhi);
            muonRecoTree->SetBranchAddress("fPhimother", &fPhi_mother);
            muonRecoTree->SetBranchAddress("fEtaassoc", &fEta);
            muonRecoTree->SetBranchAddress("fEtamother", &fEta_mother);
            std::unordered_set<Long64_t> seenMCMuons;

            // First pass: build groups of muons from the same event
            for (Long64_t i = 0; i < muonRecoTree->GetEntries(); ++i) {
                muonRecoTree->GetEntry(i);

                // if (seenMCMuons.count(fGlobalIndexMCtrack)) continue; // Skip if we've already seen this MC track index
                if (!(charm_beauty_cut(fMotherPDG) || charm_beauty_cut(fGrandmotherPDG))) continue;
                if (fPt < pT_mu_min || fPt > pT_mu_max) continue;
                if (fEta < eta_mu_min || fEta > eta_mu_max) continue;

                // Save muon
                muon_groups[fEventIdx].push_back(i);
                seenMCMuons.insert(fGlobalIndexMCtrack);
                pTMuonReco = fPt;
                pTMuonReco_true = fPt_true;
                etaMuonReco = fEta;
                outTreeMuonsReco->Fill();
            }

            // Prepare to label
            Long64_t fGlobalIndexAssoc;
            muonRecoTree->SetBranchAddress("fGlobalIndexassoc", &fGlobalIndexAssoc);

            // Second pass: process each group
            for (auto& event : muon_groups) {

                ULong64_t eventID = event.first;
                auto& muon_entries = event.second;

                if (muon_entries.size() < 2) continue; // Needs at least 2 muons to form a pair

                std::vector<ROOT::Math::PtEtaPhiMVector> muon_vectors;
                std::vector<ROOT::Math::PtEtaPhiMVector> muon_mother_vectors;
                std::vector<Long64_t> motherIDs;

                // Read muon kinematics and build 4-vectors
                for (auto entry : muon_entries) {
                    muonRecoTree->GetEntry(entry);
                    if (!(fMotherPDG == 443)) continue; // Only consider muons from J/Psi for the JPsi tree
                    ROOT::Math::PtEtaPhiMVector muon_vec(fPt, fEta, fPhi, 0.105658); // Muon mass ~105.658 MeV/c^2
                    muon_vectors.push_back(muon_vec);
                    ROOT::Math::PtEtaPhiMVector muon_mother_vec(fPt_mother, fEta_mother, fPhi_mother, 3.096916); // J/Psi mass ~3.096916 GeV/c^2
                    muon_mother_vectors.push_back(muon_mother_vec);
                    motherIDs.push_back(fMotherID);
                }

                // Find muon pairs with invariant mass closest to J/Psi mass
                for (size_t j = 0; j < muon_vectors.size(); ++j) {
                    for (size_t k = j + 1; k < muon_vectors.size(); ++k) {

                        if (!(motherIDs[j] == motherIDs[k])) continue; // Check same mother
                        auto track = muon_vectors[j] + muon_vectors[k];
                        auto track_true = muon_mother_vectors[j];

                        // if (!((track.M() > signal_range_min) && (track.M() < signal_range_max))) continue; // Apply invariant mass cut to select J/Psi candidates
                        if (!((track.Eta() > eta_JPsi_min && track.Eta() < eta_JPsi_max) && (track.Pt() > pT_JPsi_min && track.Pt() < pT_JPsi_max))) continue; // Cuts
                        // if (!((track_true.Eta() > eta_JPsi_min && track_true.Eta() < eta_JPsi_max) && (track_true.Pt() > pT_JPsi_min && track_true.Pt() < pT_JPsi_max))) continue; // Cuts
                        pTJPsiReco = track.Pt();
                        pTJPsiReco_true = track_true.Pt();
                        etaJPsiReco = track.Eta ();
                        massJPsiReco = track.M();
                        outTreeJPsiReco->Fill();
                    }
                }
            }

            dirCount++;
        }

        std::cout << std::endl;

        delete genKey;
        genFile->Close();
        delete genFile;

        delete recoKey;
        recoFile->Close();
        delete recoFile;
    }

    outFile->cd();
    outTreeJPsiGen->Write();
    outTreeMuonsGen->Write();
    outTreeJPsiReco->Write();
    outTreeMuonsReco->Write();
    metaData->Write();
    outFile->Close();
}