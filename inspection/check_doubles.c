
struct EventEntry {
    std::vector<Long64_t> muons;
    std::map<int, std::vector<Long64_t>> neighbour_muons;
};

void check_doubles() {

    // TString MC_name = "c3_global";
    // TString MC_name = "c3_global_temp";
    TString MC_name = "c3_standalone";

    int depth = 10;

    TString reco_file = TString::Format("results/%s/reco/muonAOD0.root", MC_name.Data());
    TFile *recoFile = TFile::Open(reco_file);

    // Setup variables
    ULong64_t fEventIdx;
    Long64_t fMCEventIdx;
    Long64_t fGlobalIndexMCTrack;
    Long64_t fGlobalIndexassoc;
    Long64_t fTrackPDG;

    TIter nextRecoKey(recoFile->GetListOfKeys());
    TKey* recoKey;
    
    while ((recoKey = (TKey*) nextRecoKey())) {
            
        int dirCount = 0;
        int muonCount = 0;
        std::map<ULong64_t, EventEntry> muon_groups; // Map of event ID to list of muon entries 

        // Load directory and tree
        TObject* obj = recoKey->ReadObj();
        if (!(obj->InheritsFrom("TDirectory"))) continue;

        TDirectory* recoDir = (TDirectory*) obj;
        TTree *muonRecoTree = (TTree*)recoDir->Get("O2dqmuontable");

        muonRecoTree->SetBranchAddress("fEventIdx", &fEventIdx);
        muonRecoTree->SetBranchAddress("fMCEventIdx", &fMCEventIdx);
        muonRecoTree->SetBranchAddress("fGlobalIndexMCtrack", &fGlobalIndexMCTrack);
        muonRecoTree->SetBranchAddress("fGlobalIndexassoc", &fGlobalIndexassoc);
        muonRecoTree->SetBranchAddress("fTrackPDG", &fTrackPDG);

        for (Long64_t i = 0; i < muonRecoTree->GetEntries(); ++i) {
            muonRecoTree->GetEntry(i);
            muon_groups[fEventIdx].muons.push_back(fGlobalIndexMCTrack);
            for (int d = 1; d <= depth; ++d) {
                muon_groups[fEventIdx-d].neighbour_muons[d].push_back(fGlobalIndexMCTrack);
                muon_groups[fEventIdx+d].neighbour_muons[-d].push_back(fGlobalIndexMCTrack);
            }
            muonCount++;
        }

        dirCount++;
        std::cout << "Found muon groups for " << muon_groups.size() << " events." << std::endl;
        std::cout << "Total muons found: " << muonCount << std::endl;
        
        int double_count = 0;

        // Second loop over groups
        for (const auto& event : muon_groups) {
            ULong64_t eventIdx = event.first;
            std::vector<Long64_t> muons = event.second.muons;
            
            // Check same event
            if (muons.size() > 1) {
                for (int i = 0; i < muons.size(); ++i) {
                    for (int j = i+1; j < muons.size(); ++j) {
                        if (muons[i] == muons[j]) {
                            std::cout << TString::Format("Found muons with same MC track idx %lld in event %llu", muons[i], eventIdx) << std::endl;
                            double_count++;
                        }
                    }
                }
            }
    
            // Check neighbours
            for (const auto& neighbour : event.second.neighbour_muons) {
                int neighbour_idx = neighbour.first;
                const std::vector<Long64_t>& neighbour_muons = neighbour.second;

                for (const auto& muon : muons) {
                    for (const auto& neighbour_muon : neighbour_muons) {
                        if (muon == neighbour_muon) {
                            if (std::abs(neighbour_idx) > 3) {
                                std::cout << TString::Format("Found muon with same MC track idx %lld in event %llu and neighbour event at depth %d", muon, eventIdx, neighbour_idx) << std::endl;
                            }
                            if (neighbour_idx > 0) {
                                double_count++;
                            }
                        }
                    }
                }
            }
        }

        std::cout << std::endl << TString::Format("   Found %d muons with same MC track idx within neighbour events within depth %d", double_count, depth) << std::endl;
        std::cout << TString::Format("   This gives a rate of %f", float(double_count)/float(muonCount)) << std::endl << std::endl;
    }
}
