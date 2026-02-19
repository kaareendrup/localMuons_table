
void analysis_triggers_GEN() {

    TString MC_name = "DQ";
    // TString MC_name = "HF";
    // TString MC_name = "genpurp";

    TString type = "gen";

    float eta_trigger_min = -3.6;
    float eta_trigger_max = -2.5;

    int n_files = 94;
    TString data_file;
    
    int triggers_JPsi;
    int triggers_Psi2S;
    std::vector<double> deltaEta_JPsi, deltaPhi_JPsi;
    std::vector<double> deltaEta_Psi2S, deltaPhi_Psi2S;
    std::vector<double> deltaEta_JPsi_charm, deltaPhi_JPsi_charm;
    std::vector<double> deltaEta_JPsi_noncharm, deltaPhi_JPsi_noncharm;
    
    TFile* outFile = TFile::Open(TString::Format("results/%s/%s/eventmuons.root", MC_name.Data(), type.Data()), "RECREATE");
    TTree* outTree = new TTree("Triggers", "JPsi triggers");

    std::string category; 
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

    for (int i = 0; i < n_files; ++i) {

        std::cout << "Processing file " << i << " of " << n_files << std::endl;
        data_file = TString::Format("results/%s/%s/multi/muonAOD%d.root", MC_name.Data(), type.Data(), i);

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
            TTree *jpsitree = (TTree*)dir->Get("O2dqjpsitable");
            TTree *muontree = (TTree*)dir->Get("O2dqmuontable");

            std::cout << TString::Format("Reading tracks from dir %d of %d: %s\r", dirCount, file->GetListOfKeys()->GetEntries(), dir->GetName()) << std::flush;

            std::map<ULong64_t, std::vector<Long64_t>> jpsi_groups, muon_groups;
            ULong64_t fEventIdx;

            // First loops to match muons and J/psi by event index
            jpsitree->SetBranchAddress("fEventIdx", &fEventIdx);
            for (Long64_t i = 0; i < jpsitree->GetEntries(); ++i) {
                jpsitree->GetEntry(i);
                jpsi_groups[fEventIdx].push_back(i);
            }

            muontree->SetBranchAddress("fEventIdx", &fEventIdx);
            for (Long64_t i = 0; i < muontree->GetEntries(); ++i) {
                muontree->GetEntry(i);
                muon_groups[fEventIdx].push_back(i);
            }

            // Initialize variables to read MC truth info
            Long64_t fGlobalIndexJPsi, fPDGJPsi, fMotherIDMuon, fMotherPDGMuon;
            jpsitree->SetBranchAddress("fGlobalIndexMCtrack",   &fGlobalIndexJPsi);
            jpsitree->SetBranchAddress("fTrackPDG",             &fPDGJPsi);
            muontree->SetBranchAddress("fMotherID",             &fMotherIDMuon);
            muontree->SetBranchAddress("fMotherPDG",            &fMotherPDGMuon);

            float fPtJPsi, fPhiJPsi, fEtaJPsi;
            jpsitree->SetBranchAddress("fPtassoc",  &fPtJPsi);
            jpsitree->SetBranchAddress("fPhiassoc", &fPhiJPsi);
            jpsitree->SetBranchAddress("fEtaassoc", &fEtaJPsi);

            float fPtMuon, fPhiMuon, fEtaMuon;
            muontree->SetBranchAddress("fPtassoc",  &fPtMuon);
            muontree->SetBranchAddress("fPhiassoc", &fPhiMuon);
            muontree->SetBranchAddress("fEtaassoc", &fEtaMuon);

            for (auto& event : jpsi_groups) {

                category = "";
                pT = -9999; eta = -9999; phi = -9999; mass = 3.0969; // Set mass to J/Psi mass

                ULong64_t eventID = event.first;
                auto& jpsi_entries = event.second;
                auto& muon_entries = muon_groups[eventID];
                
                for (auto entry : jpsi_entries) {

                    jpsitree->GetEntry(entry);
                    
                    pT = fPtJPsi;
                    eta = fEtaJPsi;
                    phi = fPhiJPsi;

                    pT_assocs.clear();
                    eta_assocs.clear();
                    phi_assocs.clear();
                    MotherPID.clear();

                    if (eta < eta_trigger_min || eta > eta_trigger_max) continue; // Apply eta cut on J/Psi

                    for (auto muon_entry : muon_entries) {
                        muontree->GetEntry(muon_entry);

                        if (fMotherIDMuon == fGlobalIndexJPsi) continue; // Skip if muon is from trigger J/psi
                        
                        MotherPID.push_back(fMotherPDGMuon);
                        pT_assocs.push_back(fPtMuon);
                        eta_assocs.push_back(fEtaMuon);
                        phi_assocs.push_back(fPhiMuon);
                    }

                    if (fPDGJPsi == 443) {
                        category = "JPsi";
                    } else if (fPDGJPsi == 100443) {
                        category = "Psi2S";
                    }
                    
                    outTree->Fill();
                }
            }

            delete muontree;
            delete jpsitree;
            delete dir;

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