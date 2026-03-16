
void inspect_counts() {

    TString MC_name = "c3_global";

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

    int n_files = 25;
    
    int n_muons_gen = 0;
    int n_muons_reco = 0;
    int n_jpsi_gen = 0;

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

        while ((genKey = (TKey*) nextGenKey())) {

            // Load directory and tree
            TObject* obj = genKey->ReadObj();
            if (!(obj->InheritsFrom("TDirectory"))) continue;

            TDirectory* genDir = (TDirectory*) obj;
            TTree *jpsiGenTree = (TTree*)genDir->Get("O2dqjpsitable");
            TTree *muonGenTree = (TTree*)genDir->Get("O2dqmuontable");

            std::cout << TString::Format("Reading tracks from dir %d of %d: %s\r", dirCount, recoFile->GetListOfKeys()->GetEntries(), genDir->GetName()) << std::flush;
            
            for (Long64_t i = 0; i < jpsiGenTree->GetEntries(); ++i) {
                jpsiGenTree->GetEntry(i);
                n_jpsi_gen++;
            }

            for (Long64_t i = 0; i < muonGenTree->GetEntries(); ++i) {
                muonGenTree->GetEntry(i);
                n_muons_gen++;
            }

            dirCount++;
        }

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

            // First pass: build groups of muons from the same event
            Long64_t n = muonRecoTree->GetEntries();
            for (Long64_t i = 0; i < n; ++i) {
                muonRecoTree->GetEntry(i);
                n_muons_reco++;
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

    std::cout << "Total muons gen: " << n_muons_gen << std::endl;
    std::cout << "Total muons reco: " << n_muons_reco << std::endl;
    std::cout << "Total J/Psi gen: " << n_jpsi_gen << std::endl;
}