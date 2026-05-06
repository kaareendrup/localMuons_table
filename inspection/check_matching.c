
struct eventLookup {
    ULong64_t eventID;
    std::string dirName;
};

eventLookup check_matching_output() {

    // TString MC_name = "c3_global";
    TString MC_name = "c3_global_temp";

    TString reco_file = TString::Format("results/%s/reco/muonAOD0.root", MC_name.Data());
    TString gen_file = TString::Format("results/%s/gen/muonAOD0.root", MC_name.Data());

    TFile *recoFile = TFile::Open(reco_file);
    TFile *genFile = TFile::Open(gen_file);

    // Setup variables
    ULong64_t fEventIdx;
    Long64_t fMCEventIdx;
    Long64_t fGlobalIndexMCTrack;
    Long64_t fGlobalIndexassoc;
    Long64_t fTrackPDG;

    // Load the dataframe keys
    TIter nextGenKey(genFile->GetListOfKeys());
    TKey* genKey;

    std::map<ULong64_t, std::vector<Long64_t>> muon_groups_gen; // Map of event ID to list of muon entries 

    int dirCount = 0;
    std::string target_dir_name;

    // Loop over directories in gen file
    // while ((genKey = (TKey*) nextGenKey())) {
    while ((genKey = (TKey*) nextGenKey()) && dirCount < 1) {

        // Load directory and tree
        TObject* obj = genKey->ReadObj();
        if (!(obj->InheritsFrom("TDirectory"))) continue;

        TDirectory* genDir = (TDirectory*) obj;
        TTree *muonGenTree = (TTree*)genDir->Get("O2dqmuontable");

        std::cout << TString::Format("Reading tracks from dir %d of %d: %s\r", dirCount, genFile->GetListOfKeys()->GetEntries(), genDir->GetName()) << std::flush;
        target_dir_name = genDir->GetName();

        muonGenTree->SetBranchAddress("fEventIdx", &fEventIdx);
        muonGenTree->SetBranchAddress("fGlobalIndexMCtrack", &fGlobalIndexMCTrack);

        for (Long64_t i = 0; i < muonGenTree->GetEntries(); ++i) {
            muonGenTree->GetEntry(i);
            muon_groups_gen[fEventIdx].push_back(fGlobalIndexMCTrack);
        }
    
        dirCount++;
    }
    std::cout << std::endl;

    dirCount = 0;

    TIter nextRecoKey(recoFile->GetListOfKeys());
    TKey* recoKey;

    eventLookup target_event;

    // while ((recoKey = (TKey*) nextRecoKey())) {
    while ((recoKey = (TKey*) nextRecoKey()) && dirCount < 1) {

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
            bool match_found = false;
            for (const auto& idx : muon_groups_gen[fEventIdx]) {
                if (idx == fGlobalIndexMCTrack) {
                    match_found = true;
                    break;
                }
            }
            for (const auto& idx : muon_groups_gen[fEventIdx-1]) {
                if (idx == fGlobalIndexMCTrack) {
                    match_found = true;
                    break;
                }
            }
            for (const auto& idx : muon_groups_gen[fEventIdx+1]) {
                if (idx == fGlobalIndexMCTrack) {
                    match_found = true;
                    break;
                }
            }
            for (const auto& idx : muon_groups_gen[fEventIdx-2]) {
                if (idx == fGlobalIndexMCTrack) {
                    match_found = true;
                    break;
                }
            }
            for (const auto& idx : muon_groups_gen[fEventIdx+2]) {
                if (idx == fGlobalIndexMCTrack) {
                    match_found = true;
                    break;
                }
            }
            if (!match_found) {
                std::cout << "Match not found for event " << fEventIdx << " / " << fMCEventIdx << " and muon " << fGlobalIndexMCTrack << " / " << fGlobalIndexassoc << std::endl;
                std::cout << "Gen muons for this event: ";
                for (const auto& idx : muon_groups_gen[fEventIdx]) {
                    std::cout << idx << " ";
                }
                std::cout << std::endl;
                std::cout << "Muons 1 before: ";
                for (const auto& idx : muon_groups_gen[fEventIdx-1]) {
                    std::cout << idx << " ";
                }
                std::cout << std::endl;
                std::cout << "Muons 1 after: ";
                for (const auto& idx : muon_groups_gen[fEventIdx+1]) {
                    std::cout << idx << " ";
                }
                std::cout << std::endl;
                std::cout << "Muons 2 before: ";
                for (const auto& idx : muon_groups_gen[fEventIdx-2]) {
                    std::cout << idx << " ";
                }
                std::cout << std::endl;
                std::cout << "Muons 2 after: ";
                for (const auto& idx : muon_groups_gen[fEventIdx+2]) {
                    std::cout << idx << " ";
                }
                std::cout << std::endl;

                std::cout << "Muon PID: " << fTrackPDG << std::endl << std::endl;
                target_event = eventLookup(fEventIdx, target_dir_name);
                std::cout << "Target event set to: " << target_event.eventID << " in directory " << target_event.dirName << std::endl;
                break; // Exit loop after first mismatch for this event
            }
        }

        dirCount++;
    }
    return target_event;
}

void check_matching_input(eventLookup target_event) {

    TString data_file = "input_data/c3_global/hy_4205648/AOD/001/AO2D.root";

    TFile *dataFile = TFile::Open(data_file);

    // Setup variables
    Int_t fIndexReducedMCEvents;
    Int_t fIndexReducedMCEventsColl;
    Int_t f;
    Int_t fPdgCode;

    // Load the dataframe keys
    TIter nextDataKey(dataFile->GetListOfKeys());
    TKey* dataKey;

    std::map<Int_t, std::vector<std::tuple<Int_t, Long64_t>>> track_groups_gen; // Map of event ID to list of muon entries 

    int dirCount = 0;

    // Loop over directories in gen file
    while ((dataKey = (TKey*) nextDataKey())) {

        // Load directory and tree
        TObject* obj = dataKey->ReadObj();
        if (!(obj->InheritsFrom("TDirectory"))) continue;

        TDirectory* dataDir = (TDirectory*) obj;
        TTree *MCTrackTree = (TTree*)dataDir->Get("O2reducedmctrack"); 
        TTree *MCTrackTreeColl = (TTree*)dataDir->Get("O2remccollbl"); 

        if (dataDir->GetName() != target_event.dirName) continue;
        std::cout << TString::Format("Reading tracks from dir %d of %d: %s\r", dirCount, dataFile->GetListOfKeys()->GetEntries(), dataDir->GetName()) << std::flush;
        MCTrackTree->SetBranchAddress("fPdgCode", &fPdgCode);
        MCTrackTree->SetBranchAddress("fIndexReducedMCEvents", &fIndexReducedMCEvents);

        MCTrackTreeColl->SetBranchAddress("fIndexReducedMCEvents", &fIndexReducedMCEventsColl);

        MCTrackTreeColl->GetEntry(target_event.eventID-1); // Assuming eventID starts from 1 and corresponds to entry index in the collection tree
        std::cout << std::endl;
        std::cout << "Event ID: " << target_event.eventID << " has fIndexReducedMCEventsColl: " << fIndexReducedMCEventsColl << std::endl;

        for (Long64_t i = 0; i < MCTrackTree->GetEntries(); ++i) {
            MCTrackTree->GetEntry(i);
            if (fIndexReducedMCEvents == fIndexReducedMCEventsColl) {
                track_groups_gen[fIndexReducedMCEvents].push_back(std::make_tuple(fPdgCode,i));
            }
        }
    }
    std::cout << std::endl;

    for (const auto& event : track_groups_gen) {
        Int_t eventID = event.first;
        auto& pdgCodes = event.second;
        std::cout << "Event ID: " << eventID << " has " << pdgCodes.size() << " track with PDG codes and indexes: " << std::endl;
        for (const auto& pdg : pdgCodes) {
            std::cout << std::get<0>(pdg) << " (Index: " << std::get<1>(pdg) << ")" << std::endl;
        }
    }
}

void check_matching() {

    eventLookup target_event = check_matching_output();
    check_matching_input(target_event);

}