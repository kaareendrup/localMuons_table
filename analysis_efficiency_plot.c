
#include <nlohmann/json.hpp>
using json = nlohmann::json;

void analysis_efficiency_plot() {

    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    // Data config
    std::string data_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    TString MC_name = TString::Format("%s_%s", data_name.c_str(), muon_type.c_str());

    std::vector<double> pT_bins = config["hists"]["pT_bins"].get<std::vector<double>>();

    // Set outfile
    TFile* outFile = TFile::Open(TString::Format("results/%s/efficiency.root", MC_name.Data()), "RECREATE");
    TTree* metaDataOut = new TTree("MetaData", "Event selection metadata");

    // Set infile
    TString in_file = TString::Format("results/%s/particles.root", MC_name.Data());
    TFile* inFile = TFile::Open(in_file, "READ");
    if (!inFile || inFile->IsZombie()) {
        std::cerr << "Cannot open file\n";
        return;
    }

    // Get input trees and branches
    TTree* muonsReco = nullptr;
    TTree* JPsiReco = nullptr;
    TTree* muonsGen = nullptr;
    TTree* JPsiGen = nullptr;
    TTree* metaData = nullptr;

    inFile->GetObject("MuonsReco", muonsReco);
    inFile->GetObject("JPsiReco", JPsiReco);
    inFile->GetObject("MuonsGen", muonsGen);
    inFile->GetObject("JPsiGen", JPsiGen);
    inFile->GetObject("MetaData", metaData);

    double pTMuonReco, pTJPsiReco, pTMuonGen, pTJPsiGen, pTMuonReco_true, pTJPsiReco_true;

    muonsReco->SetBranchAddress("pTMuon", &pTMuonReco);
    muonsReco->SetBranchAddress("pTMuon_true", &pTMuonReco_true);
    JPsiReco->SetBranchAddress("pTJPsi", &pTJPsiReco);
    JPsiReco->SetBranchAddress("pTJPsi_true", &pTJPsiReco_true);
    muonsGen->SetBranchAddress("pTMuon", &pTMuonGen);
    JPsiGen->SetBranchAddress("pTJPsi", &pTJPsiGen);

    // Propagate metadata
    std::vector<float>* pTCuts = nullptr;
    std::vector<float>* etaCuts = nullptr;
    metaData->SetBranchAddress("pTCuts", &pTCuts);
    metaData->SetBranchAddress("etaCuts", &etaCuts);
    metaData->GetEntry(0);

    std::vector<float>* pTCutsOut = nullptr;
    std::vector<float>* etaCutsOut = nullptr;
    metaDataOut->Branch("pTCuts", &pTCutsOut);
    metaDataOut->Branch("etaCuts", &etaCutsOut);
    pTCutsOut = pTCuts;
    etaCutsOut = etaCuts;
    metaDataOut->Fill();

    // Create output histograms
    TH1F *pTMuonRecoHist = new TH1F("pTMuonRecoHist", "pTMuonRecoHist", pT_bins.size()-1, pT_bins.data());
    TH1F *pTMuonRecoTrueHist = new TH1F("pTMuonRecoTrueHist", "pTMuonRecoTrueHist", pT_bins.size()-1, pT_bins.data());
    TH1F *pTJPsiRecoHist = new TH1F("pTJPsiRecoHist", "pTJPsiRecoHist", pT_bins.size()-1, pT_bins.data());
    TH1F *pTJPsiRecoTrueHist = new TH1F("pTJPsiRecoTrueHist", "pTJPsiRecoTrueHist", pT_bins.size()-1, pT_bins.data());
    TH1F *pTMuonGenHist = new TH1F("pTMuonGenHist", "pTMuonGenHist", pT_bins.size()-1, pT_bins.data());
    TH1F *pTJPsiGenHist = new TH1F("pTJPsiGenHist", "pTJPsiGenHist", pT_bins.size()-1, pT_bins.data());
    
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
        pTJPsiRecoTrueHist->Fill(pTJPsiReco_true);
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
    TH1F *JPsiEffTrueHist = (TH1F*)pTJPsiRecoTrueHist->Clone("JPsiEffTrueHist");

    muonEffHist->Divide(pTMuonGenHist);
    muonEffTrueHist->Divide(pTMuonGenHist);
    JPsiEffHist->Divide(pTJPsiGenHist);
    JPsiEffTrueHist->Divide(pTJPsiGenHist);

    outFile->cd();
    pTMuonRecoHist->Write();
    pTMuonRecoTrueHist->Write();
    pTJPsiRecoHist->Write();
    pTJPsiRecoTrueHist->Write();
    pTMuonGenHist->Write();
    pTJPsiGenHist->Write();

    muonEffHist->Write();
    muonEffTrueHist->Write();
    JPsiEffHist->Write();
    JPsiEffTrueHist->Write();

    outFile->Write();
    outFile->Close();
}