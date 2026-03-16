
bool charm_beauty_cut(Long64_t motherPDG) {
    return ((std::abs(motherPDG) == 443) || (std::abs(motherPDG) == 100443) ||
            (std::abs(motherPDG) >= 411 && std::abs(motherPDG) <= 445) || 
            (std::abs(motherPDG) >= 4101 && std::abs(motherPDG) <= 4444) || 
            (std::abs(motherPDG) >= 511 && std::abs(motherPDG) <= 557) || 
            (std::abs(motherPDG) >= 5101 && std::abs(motherPDG) <= 5554));
}

void analysis_efficiency() {

    TString MC_name = "c3_global";
    // TString MC_name = "c3_standalone";
    // TString MC_name = "c3_global_temp";

    // float eta_trigger_min = -3.6;
    // float eta_trigger_max = -2.5;
    float eta_trigger_min = -10;
    float eta_trigger_max = 10;

    // Cuts
    // float pT_trigger_leg_min = 0.7;
    // float pT_trigger_leg_max = 20.0;
    float pT_trigger_leg_min = 0.;
    float pT_trigger_leg_max = 20.0;

    float eta_trigger_leg_min = -3.6;
    float eta_trigger_leg_max = -2.5;
    // float eta_trigger_leg_min = -10;
    // float eta_trigger_leg_max = 10;

    // Signal range
    // float signal_range_min = 2.7;
    // float signal_range_max = 3.4;

    int n_files = 25;
    // int n_files = 1;
    
    TFile* outFile = TFile::Open(TString::Format("results/%s/particles.root", MC_name.Data()), "RECREATE");
    TTree* outTreeMuonsReco = new TTree("MuonsReco", "Reconstructed Muons");
    TTree* outTreeJPsiReco = new TTree("JPsiReco", "Reconstructed JPsi");
    TTree* outTreeMuonsGen = new TTree("MuonsGen", "Generated Muons");
    TTree* outTreeJPsiGen = new TTree("JPsiGen", "Generated JPsi");

    double pTJPsiReco, pTJPsiGen, pTMuonReco, pTMuonGen;

    outTreeJPsiReco->Branch("pTJPsi",  &pTJPsiReco,  "pT/D");
    outTreeJPsiGen->Branch("pTJPsi",  &pTJPsiGen,  "pT/D");
    outTreeMuonsReco->Branch("pTMuon",  &pTMuonReco,  "pT/D");
    outTreeMuonsGen->Branch("pTMuon",  &pTMuonGen,  "pT/D");

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

        // Loop over gen dataframes
        while ((genKey = (TKey*) nextGenKey())) {

            // Load directory and tree
            TObject* obj = genKey->ReadObj();
            if (!(obj->InheritsFrom("TDirectory"))) continue;

            TDirectory* genDir = (TDirectory*) obj;
            TTree *jpsiGenTree = (TTree*)genDir->Get("O2dqjpsitable");
            TTree *muonGenTree = (TTree*)genDir->Get("O2dqmuontable");

            std::cout << TString::Format("Reading tracks from dir %d of %d: %s\r", dirCount, recoFile->GetListOfKeys()->GetEntries(), genDir->GetName()) << std::flush;
            
            float fPtJPsiGen, fPtMuonGen, fEtaJPsiGen, fEtaMuonGen;
            Long64_t fMotherPDG, fGrandmotherPDG;

            jpsiGenTree->SetBranchAddress("fPtassoc", &fPtJPsiGen);
            jpsiGenTree->SetBranchAddress("fEtaassoc", &fEtaJPsiGen);
            for (Long64_t i = 0; i < jpsiGenTree->GetEntries(); ++i) {
                jpsiGenTree->GetEntry(i);
                if (fEtaJPsiGen < eta_trigger_min || fEtaJPsiGen > eta_trigger_max) continue; // Apply eta cut on J/Psi
                pTJPsiGen = fPtJPsiGen;
                outTreeJPsiGen->Fill();
            }

            muonGenTree->SetBranchAddress("fPtassoc", &fPtMuonGen);
            muonGenTree->SetBranchAddress("fEtaassoc", &fEtaMuonGen);
            muonGenTree->SetBranchAddress("fMotherPDG", &fMotherPDG);
            muonGenTree->SetBranchAddress("fGrandmotherPDG", &fGrandmotherPDG);
            for (Long64_t i = 0; i < muonGenTree->GetEntries(); ++i) {
                muonGenTree->GetEntry(i);

                if (!(charm_beauty_cut(fMotherPDG) || charm_beauty_cut(fGrandmotherPDG))) continue;
                if (fEtaMuonGen < eta_trigger_leg_min || fEtaMuonGen > eta_trigger_leg_max) continue; // Apply eta cut on muons
                if (fPtMuonGen < pT_trigger_leg_min || fPtMuonGen > pT_trigger_leg_max) continue; // Apply pT cut on muons
                pTMuonGen = fPtMuonGen;
                outTreeMuonsGen->Fill();
            }

            dirCount++;
        }

        std::cout << std::endl;
        dirCount = 0;

        // Loop over reco dataframes
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
            float fPt, fPhi, fEta;
            // muonRecoTree->SetBranchAddress("fPtassoc", &fPt);
            muonRecoTree->SetBranchAddress("fPtassoctrue", &fPt);
            muonRecoTree->SetBranchAddress("fPhiassoc", &fPhi);
            muonRecoTree->SetBranchAddress("fEtaassoc", &fEta);
            
            std::unordered_set<Long64_t> seenMCMuons;

            // First pass: build groups of muons from the same event
            for (Long64_t i = 0; i < muonRecoTree->GetEntries(); ++i) {
                muonRecoTree->GetEntry(i);

                if (seenMCMuons.count(fGlobalIndexMCtrack)) continue; // Skip if we've already seen this MC track index
                if (!(charm_beauty_cut(fMotherPDG) || charm_beauty_cut(fGrandmotherPDG))) continue;
                if (fPt < pT_trigger_leg_min || fPt > pT_trigger_leg_max) continue;
                if (fEta < eta_trigger_leg_min || fEta > eta_trigger_leg_max) continue;

                // Save muon
                muon_groups[fEventIdx].push_back(i);
                seenMCMuons.insert(fGlobalIndexMCtrack);
                pTMuonReco = fPt;
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
                std::vector<Long64_t> motherIDs;

                // Read muon kinematics and build 4-vectors
                for (auto entry : muon_entries) {
                    muonRecoTree->GetEntry(entry);
                    ROOT::Math::PtEtaPhiMVector muon_vec(fPt, fEta, fPhi, 0.105658); // Muon mass ~105.658 MeV/c^2
                    muon_vectors.push_back(muon_vec);
                    motherIDs.push_back(fMotherID);
                }

                // Find muon pairs with invariant mass closest to J/Psi mass
                for (size_t j = 0; j < muon_vectors.size(); ++j) {
                    for (size_t k = j + 1; k < muon_vectors.size(); ++k) {

                        if (!(motherIDs[j] == motherIDs[k])) continue; // Check same mother
                        auto track = muon_vectors[j] + muon_vectors[k];

                        // if (!(track.M() > signal_range_min && track.M() < signal_range_max)) continue; // Check if inv mass is in signal range
                        if (!((track.Eta() > eta_trigger_min && track.Eta() < eta_trigger_max) && (track.Pt() > 0))) continue; // Cuts
                        pTJPsiReco = track.Pt();
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
    outFile->Close();
}