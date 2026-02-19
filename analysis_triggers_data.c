
void analysis_triggers_data() {

    // TString data_name = "DQ";
    TString data_name = "DQ_data";
    // TString data_name = "DQ_data_global";

    TString type = "reco";

    int n_files = 13;
    TString data_file;

    // TString data_file = "results/" + data_name + "/" + type + "/muonAOD.root";
    
    // Cuts
    float pT_trigger_leg_min = 0.7;
    float pT_trigger_leg_max = 20.0;

    float eta_trigger_leg_min = -3.6;
    float eta_trigger_leg_max = -2.5;
  
    TFile* outFile = TFile::Open(TString::Format("results/%s/%s/eventmuons.root", data_name.Data(), type.Data()), "RECREATE");
    TTree* outTree = new TTree("Triggers", "JPsi triggers");

    std::string category = "All"; 
    double pT, eta, phi, mass;
    std::vector<double> pT_assocs, eta_assocs, phi_assocs;
    std::vector<int> MotherPID;

    outTree->Branch("category", &category);
    outTree->Branch("pT",  &pT,  "pT/D");
    outTree->Branch("eta",  &eta,  "eta/D");
    outTree->Branch("phi",  &phi,  "phi/D");
    outTree->Branch("mass",  &mass,  "mass/D");
    outTree->Branch("pT_assocs",  &pT_assocs);
    outTree->Branch("eta_assocs",  &eta_assocs);
    outTree->Branch("phi_assocs",  &phi_assocs);
    outTree->Branch("MotherPID", &MotherPID);

    int exceptions[] = {5, 7, 8}; // Files with issues (e.g. missing trees)

    for (int i = 0; i < n_files; ++i) {

        // Skip files with known issues
        if (std::find(std::begin(exceptions), std::end(exceptions), i) != std::end(exceptions)) {
            std::cout << "Skipping file " << i << " due to known issues.\n";
            continue;
        }

        std::cout << "Processing file " << i << " of " << n_files << std::endl;
        data_file = TString::Format("results/%s/%s/multi/muonAOD%d.root", data_name.Data(), type.Data(), i);

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

            // Prepare to label
            ULong64_t fGlobalIndexAssoc;
            tree->SetBranchAddress("fGlobalIndexassoc", &fGlobalIndexAssoc);

            // Prepare to read muon kinematics
            float fPt, fPhi, fEta;
            tree->SetBranchAddress("fPtassoc", &fPt);
            tree->SetBranchAddress("fPhiassoc", &fPhi);
            tree->SetBranchAddress("fEtaassoc", &fEta);

            // Second pass: process each group
            for (auto& event : muon_groups) {

                pT = -9999; eta = -9999; phi = -9999; mass = -9999;

                pT_assocs.clear();
                eta_assocs.clear();
                phi_assocs.clear();
                MotherPID.clear();

                ULong64_t eventID = event.first;
                auto& muon_entries = event.second;

                if (muon_entries.size() < 2) continue; // Needs at least 2 muons to form a pair

                std::vector<Long64_t> muon_motherIDs, muon_motherPDGs;
                std::vector<ROOT::Math::PtEtaPhiMVector> muon_vectors;

                // Read muon kinematics and build 4-vectors
                for (auto entry : muon_entries) {
                    tree->GetEntry(entry);
                    ROOT::Math::PtEtaPhiMVector muon_vec(fPt, fEta, fPhi, 0.105658); // Muon mass ~105.658 MeV/c^2
                    muon_vectors.push_back(muon_vec);
                }

                // Store the indexes of the best candidate muon pair
                ROOT::Math::PtEtaPhiMVector JPsiCandidate(-9999,0,0,-1); 
                int idx_cand_1, idx_cand_2;

                // Find the muon pair with invariant mass closest to J/Psi mass
                for (size_t j = 0; j < muon_vectors.size(); ++j) {
                    if (muon_vectors[j].Pt() < pT_trigger_leg_min || muon_vectors[j].Pt() > pT_trigger_leg_max) continue;
                    if (muon_vectors[j].Eta() < eta_trigger_leg_min || muon_vectors[j].Eta() > eta_trigger_leg_max) continue;

                    for (size_t k = j + 1; k < muon_vectors.size(); ++k) {
                        if (muon_vectors[k].Pt() < pT_trigger_leg_min || muon_vectors[k].Pt() > pT_trigger_leg_max) continue;
                        if (muon_vectors[k].Eta() < eta_trigger_leg_min || muon_vectors[k].Eta() > eta_trigger_leg_max) continue;

                        auto track = muon_vectors[j] + muon_vectors[k];

                        // Replace candidate if closer to J/Psi mass
                        if (std::abs(track.M() - 3.0969) < std::abs(JPsiCandidate.M() - 3.0969)) { // J/Psi mass ~3.0969 GeV/c^2
                            JPsiCandidate = track;
                            idx_cand_1 = j;
                            idx_cand_2 = k;
                        }
                    }
                }

                if (JPsiCandidate.Pt() < 0) continue; // No valid candidate found, can happen if no pairs pass the pT cut

                // Create best candiddate from vector
                pT = JPsiCandidate.Pt();
                eta = JPsiCandidate.Eta();
                phi = JPsiCandidate.Phi();
                mass = JPsiCandidate.M();

                // Loop over muons again to find assocs
                for (size_t j = 0; j < muon_vectors.size(); ++j) {
                    if (j == idx_cand_1 || j == idx_cand_2) continue; // Skip the trigger muons
                    auto assocTrack = muon_vectors[j];
                    pT_assocs.push_back(assocTrack.Pt());
                    eta_assocs.push_back(assocTrack.Eta());
                    phi_assocs.push_back(assocTrack.Phi());
                }
                
                // Store candidate and associates
                outTree->Fill();
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

    outFile->cd();
    outTree->Write();
    outFile->Close();
}