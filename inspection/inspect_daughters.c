
void inspect_daughters() {

    TString data_file = "muonAOD0_f4d.root";
    TFile *file = TFile::Open(data_file, "READ");

    // Loop over all keys in the file (directories etc.)
    TIter nextKey(file->GetListOfKeys());
    TKey* key;
    int dirCount = 0;

    int totalJPsis = 0;
    int JPsiWithMuonDaughters = 0;

    while ((key = (TKey*) nextKey())) {
        TObject* obj = key->ReadObj();

        // --- Handle only TDirectoryFile objects ---
        if (obj->InheritsFrom("TDirectory")) {
            TDirectory* dir = (TDirectory*) obj;
            printf("Reading tracks from dir %d of %d: %s\n", dirCount, file->GetListOfKeys()->GetEntries(), dir->GetName());

            TTree* MCJPsis = (TTree*)dir->Get("O2dqjpsitable");

            Long64_t fTrackPDG, fDaughterPDGSum;
            MCJPsis->SetBranchAddress("fTrackPDG", &fTrackPDG);
            MCJPsis->SetBranchAddress("fDaughterPDGSum", &fDaughterPDGSum);

            for (Long64_t i = 0; i < MCJPsis->GetEntries(); ++i) {
                MCJPsis->GetEntry(i);
                
                if (!(std::abs(fTrackPDG) == 443)) continue;
                totalJPsis++;
                if (fDaughterPDGSum == 26) {
                    JPsiWithMuonDaughters++;
                    JPsiWithMuonDaughters++;
                }
                if (fDaughterPDGSum == 13) {
                    JPsiWithMuonDaughters++;
                }
            }

            MCJPsis->ResetBranchAddresses();
            delete MCJPsis;
        }

        dirCount++;
    }

    std::cout << std::endl << TString::Format("%d J/Psis in total, %d with muon daughters, which is %f%%", totalJPsis, JPsiWithMuonDaughters, (double)JPsiWithMuonDaughters / totalJPsis * 100) << std::endl;

    file->Close();
}