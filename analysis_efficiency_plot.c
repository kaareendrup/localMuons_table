void analysis_efficiency_plot() {

    TString MC_name = "c3_global";
    // TString MC_name = "c3_standalone";
    // TString MC_name = "c3_global_temp";

    int n_bins = 20;
    float x_min = 0;
    float x_max = 20;

    // Set outfile
    TFile* outFile = TFile::Open(TString::Format("results/%s/efficiency.root", MC_name.Data()), "RECREATE");
    TTree* pTHists = new TTree("pTHists", "pT histograms");

    TString in_file = TString::Format("results/%s/particles.root", MC_name.Data());
    TFile* inFile = TFile::Open(in_file, "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Cannot open file\n";
        return;
    }

    TTree* muonsReco = nullptr;
    TTree* JPsiReco = nullptr;
    TTree* muonsGen = nullptr;
    TTree* JPsiGen = nullptr;

    inFile->GetObject("MuonsReco", muonsReco);
    inFile->GetObject("JPsiReco", JPsiReco);
    inFile->GetObject("MuonsGen", muonsGen);
    inFile->GetObject("JPsiGen", JPsiGen);

    double pTMuonReco, pTJPsiReco, pTMuonGen, pTJPsiGen, pTMuonReco_true;

    muonsReco->SetBranchAddress("pTMuon",  &pTMuonReco);
    muonsReco->SetBranchAddress("pTMuon_true",  &pTMuonReco_true);
    JPsiReco->SetBranchAddress("pTJPsi",  &pTJPsiReco);
    muonsGen->SetBranchAddress("pTMuon",  &pTMuonGen);
    JPsiGen->SetBranchAddress("pTJPsi",  &pTJPsiGen);

    TH1F *pTMuonRecoHist = new TH1F("pTMuonRecoHist", "pTMuonRecoHist", n_bins, x_min, x_max);
    TH1F *pTMuonRecoTrueHist = new TH1F("pTMuonRecoTrueHist", "pTMuonRecoTrueHist", n_bins, x_min, x_max);
    TH1F *pTJPsiRecoHist = new TH1F("pTJPsiRecoHist", "pTJPsiRecoHist", n_bins, x_min, x_max);
    TH1F *pTMuonGenHist = new TH1F("pTMuonGenHist", "pTMuonGenHist", n_bins, x_min, x_max);
    TH1F *pTJPsiGenHist = new TH1F("pTJPsiGenHist", "pTJPsiGenHist", n_bins, x_min, x_max);
    
    // Loop over entries and create the necessary histograms
    for (Long64_t i = 0; i < muonsReco->GetEntries(); ++i) {
        std::cout << "Processing reco muons entry " << i+1 << " of " << muonsReco->GetEntries() << "\r" << std::flush;
        muonsReco->GetEntry(i);
        pTMuonRecoHist->Fill(pTMuonReco);
        pTMuonRecoTrueHist->Fill(pTMuonReco_true);
    }
    std::cout << std::endl;

    for (Long64_t i = 0; i < JPsiReco->GetEntries(); ++i) {
        std::cout << "Processing reco J/Psi entry " << i+1 << " of " << JPsiReco->GetEntries() << "\r" << std::flush;
        JPsiReco->GetEntry(i);
        pTJPsiRecoHist->Fill(pTJPsiReco);
    }
    std::cout << std::endl;

    for (Long64_t i = 0; i < muonsGen->GetEntries(); ++i) {
        std::cout << "Processing gen muons entry " << i+1 << " of " << muonsGen->GetEntries() << "\r" << std::flush;
        muonsGen->GetEntry(i);
        pTMuonGenHist->Fill(pTMuonGen);
    }
    std::cout << std::endl;

    for (Long64_t i = 0; i < JPsiGen->GetEntries(); ++i) {
        std::cout << "Processing gen J/Psi entry " << i+1 << " of " << JPsiGen->GetEntries() << "\r" << std::flush;        
        JPsiGen->GetEntry(i);
        pTJPsiGenHist->Fill(pTJPsiGen);
    }
    std::cout << std::endl;

    TH1F *muonEffHist = (TH1F*)pTMuonRecoHist->Clone("muonEffHist");
    TH1F *muonEffTrueHist = (TH1F*)pTMuonRecoTrueHist->Clone("muonEffTrueHist");
    TH1F *JPsiEffHist = (TH1F*)pTJPsiRecoHist->Clone("JPsiEffHist");

    muonEffHist->Divide(pTMuonGenHist);
    muonEffTrueHist->Divide(pTMuonGenHist);
    JPsiEffHist->Divide(pTJPsiGenHist);

    outFile->cd();
    pTMuonRecoHist->Write();
    pTMuonRecoTrueHist->Write();
    pTJPsiRecoHist->Write();
    pTMuonGenHist->Write();
    pTJPsiGenHist->Write();

    muonEffHist->Write();
    muonEffTrueHist->Write();
    JPsiEffHist->Write();
    outFile->Close();
}