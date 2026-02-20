
void analysis_efficiency() {

    TString MC_name = "DQ";
    // TString MC_name = "HF";
    // TString MC_name = "genpurp";

    float eta_trigger_min = -3.6;
    float eta_trigger_max = -2.5;

    // Cuts
    float pT_trigger_leg_min = 0.7;
    float pT_trigger_leg_max = 20.0;

    float eta_trigger_leg_min = -3.6;
    float eta_trigger_leg_max = -2.5;

    // Signal range
    float signal_range_min = 2.7;
    float signal_range_max = 3.4;

    int n_files = 94;
    
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
        TString reco_file = TString::Format("results/%s/reco/multi/muonAOD%d.root", MC_name.Data(), i);
        TString gen_file = TString::Format("results/%s/gen/multi/muonAOD%d.root", MC_name.Data(), i);
    
        TFile *recoFile = TFile::Open(reco_file);
        TFile *genFile = TFile::Open(gen_file);

        // Load the dataframe keys
        TIter nextGenKey(genFile->GetListOfKeys());
        TKey* genKey;

        TIter nextRecoKey(recoFile->GetListOfKeys());
        TKey* recoKey;
        
        int dirCount = 0;

        while ((genKey = (TKey*) nextGenKey())) {

            // Load directory and tree
            TObject* obj = genKey->ReadObj();
            if (!(obj->InheritsFrom("TDirectory"))) continue;

            TDirectory* genDir = (TDirectory*) obj;
            TTree *jpsiGenTree = (TTree*)genDir->Get("O2dqjpsitable");
            TTree *muonGenTree = (TTree*)genDir->Get("O2dqmuontable");

            std::cout << TString::Format("Reading tracks from dir %d of %d: %s\r", dirCount, recoFile->GetListOfKeys()->GetEntries(), genDir->GetName()) << std::flush;
            
            float fPtJPsiGen, fPtMuonGen, fEtaJPsiGen, fEtaMuonGen;

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
            for (Long64_t i = 0; i < muonGenTree->GetEntries(); ++i) {
                muonGenTree->GetEntry(i);
                if (fEtaMuonGen < eta_trigger_leg_min || fEtaMuonGen > eta_trigger_leg_max) continue; // Apply eta cut on muons
                pTMuonGen = fPtMuonGen;
                outTreeMuonsGen->Fill();
            }

            dirCount++;
        }

        outFile->cd();
        outTreeJPsiGen->Write();
        outTreeMuonsGen->Write();

        std::cout << std::endl;
        dirCount = 0;

        // Loop over dataframes
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
            muonRecoTree->SetBranchAddress("fEventIdx", &fEventIdx);

            // Prepare to read muon kinematics
            float fPt, fPhi, fEta;
            muonRecoTree->SetBranchAddress("fPtassoc", &fPt);
            muonRecoTree->SetBranchAddress("fPhiassoc", &fPhi);
            muonRecoTree->SetBranchAddress("fEtaassoc", &fEta);
            
            // First pass: build groups of muons from the same event
            Long64_t n = muonRecoTree->GetEntries();
            for (Long64_t i = 0; i < n; ++i) {
                muonRecoTree->GetEntry(i);
                if (fPt < pT_trigger_leg_min || fPt > pT_trigger_leg_max) continue;
                if (fEta < eta_trigger_leg_min || fEta > eta_trigger_leg_max) continue;

                muon_groups[fEventIdx].push_back(i);

                // Save muon
                pTMuonReco = fPt;
                outTreeMuonsReco->Fill();
            }

            // Prepare to label
            ULong64_t fGlobalIndexAssoc;
            muonRecoTree->SetBranchAddress("fGlobalIndexassoc", &fGlobalIndexAssoc);

            // Second pass: process each group
            for (auto& event : muon_groups) {

                ULong64_t eventID = event.first;
                auto& muon_entries = event.second;

                if (muon_entries.size() < 2) continue; // Needs at least 2 muons to form a pair

                std::vector<ROOT::Math::PtEtaPhiMVector> muon_vectors;

                // Read muon kinematics and build 4-vectors
                for (auto entry : muon_entries) {
                    muonRecoTree->GetEntry(entry);
                    ROOT::Math::PtEtaPhiMVector muon_vec(fPt, fEta, fPhi, 0.105658); // Muon mass ~105.658 MeV/c^2
                    muon_vectors.push_back(muon_vec);
                }

                // Find muon pairs with invariant mass closest to J/Psi mass
                for (size_t j = 0; j < muon_vectors.size(); ++j) {
                    for (size_t k = j + 1; k < muon_vectors.size(); ++k) {

                        auto track = muon_vectors[j] + muon_vectors[k];

                        if (!(track.M() > signal_range_min && track.M() < signal_range_max)) continue; // Check if inv mass is in signal range
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
    outTreeJPsiReco->Write();
    outTreeMuonsReco->Write();
    outFile->Close();
}