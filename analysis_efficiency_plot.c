
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TEfficiency.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <map>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "utils/hists.c"

void analysis_efficiency_plot() {

    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    // Data config
    std::string data_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    TString MC_name = TString::Format("%s_%s", data_name.c_str(), muon_type.c_str());
    
    float signal_range_min = config["signal_range"]["min"];
    float signal_range_max = config["signal_range"]["max"];
    float background_range_min = config["background_range"]["min"];
    float background_range_max = config["background_range"]["max"];

    float JPsi_branching_ratio = config["config_dataset"][data_name + "_" + muon_type]["branching_ratio"];
    std::vector<double> pT_bins = config["hists"]["pT_bins"].get<std::vector<double>>();
    const int n_pT_bins = pT_bins.size() - 1;

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
    double massJPsiReco, massJPsiGen;

    muonsReco->SetBranchAddress("pTMuon", &pTMuonReco);
    muonsReco->SetBranchAddress("pTMuon_true", &pTMuonReco_true);
    JPsiReco->SetBranchAddress("pTJPsi", &pTJPsiReco);
    JPsiReco->SetBranchAddress("massJPsi", &massJPsiReco);
    JPsiReco->SetBranchAddress("pTJPsi_true", &pTJPsiReco_true);
    muonsGen->SetBranchAddress("pTMuon", &pTMuonGen);
    JPsiGen->SetBranchAddress("pTJPsi", &pTJPsiGen);
    JPsiGen->SetBranchAddress("massJPsi", &massJPsiGen);

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
    TH1F *pTMuonRecoTrueRecoHist = new TH1F("pTMuonRecoTrueRecoHist", "pTMuonRecoTrueRecoHist", 20, 0, .5);
    TH1F *pTJPsiRecoHist = new TH1F("pTJPsiRecoHist", "pTJPsiRecoHist", pT_bins.size()-1, pT_bins.data());
    TH1F *pTJPsiRecoTrueHist = new TH1F("pTJPsiRecoTrueHist", "pTJPsiRecoTrueHist", pT_bins.size()-1, pT_bins.data());
    TH1F *pTMuonGenHist = new TH1F("pTMuonGenHist", "pTMuonGenHist", pT_bins.size()-1, pT_bins.data());
    TH1F *pTJPsiGenHist = new TH1F("pTJPsiGenHist", "pTJPsiGenHist", pT_bins.size()-1, pT_bins.data());

    // Set up for resolution plot
    TH1F* h_ptRes[n_pT_bins];
    TH1F* h_ptResAbs[n_pT_bins];

    for (uint i = 0; i < pT_bins.size()-1; ++i) {
        h_ptRes[i] = new TH1F(
            Form("h_ptRes_bin%d", i),
            Form("pT resolution; (p_{T}^{true} - p_{T}^{reco}) / p_{T}^{true}; Events"),
            100, -0.5, 0.5
        );
        h_ptResAbs[i] = new TH1F(
            Form("h_ptResAbs_bin%d", i),
            Form("Absolute pT resolution; |p_{T}^{true} - p_{T}^{reco}| / p_{T}^{true}; Events"),
            100, 0, 0.5
        );
    }

    // Loop over entries and create the necessary histograms
    for (Long64_t i = 0; i < muonsReco->GetEntries(); ++i) {
        std::cout << "Processing reco muons entry " << i+1 << " of " << muonsReco->GetEntries() << "\r" << std::flush;
        muonsReco->GetEntry(i);
        pTMuonRecoHist->Fill(pTMuonReco);
        pTMuonRecoTrueHist->Fill(pTMuonReco_true);
        if (pTMuonReco_true < 0 || pTMuonReco < 0) {
            std::cerr << "Warning: pTReco or pTGen are negative for entry " << i << " in file " << in_file << "\n";
        }

        pTMuonRecoTrueRecoHist->Fill(std::abs(pTMuonReco_true-pTMuonReco)/pTMuonReco_true);
        for (uint i = 0; i < pT_bins.size()-1; ++i) {
            if (pTMuonReco_true >= pT_bins[i] && pTMuonReco_true < pT_bins[i+1]) {
                double ptRes = (pTMuonReco_true - pTMuonReco) / pTMuonReco_true;
                double ptResAbs =  std::abs(ptRes);
                h_ptRes[i]->Fill(ptRes);
                h_ptResAbs[i]->Fill(ptResAbs);
                break;
            }
        }
    }
    std::cout << std::endl;
    std::cout << "Entries in reco pT histogram: " << pTMuonRecoHist->GetEntries() << ", Entries in reco true pT histogram: " << pTMuonRecoTrueHist->GetEntries() << std::endl;

    TGraphErrors* gr_res = new TGraphErrors(pT_bins.size()-1);
    TGraphErrors* gr_mean = new TGraphErrors(pT_bins.size()-1);

    for (uint i = 0; i < pT_bins.size()-1; ++i) {
        h_ptRes[i]->Fit("gaus", "Q");  // Quiet fit

        TF1* f = h_ptRes[i]->GetFunction("gaus");
        double sigma = f->GetParameter(2);
        double sigmaErr = f->GetParError(2);

        double ptCenter = 0.5 * (pT_bins[i] + pT_bins[i+1]);
        double ptWidth  = 0.5 * (pT_bins[i+1] - pT_bins[i]);

        gr_res->SetPoint(i, ptCenter, sigma);
        gr_res->SetPointError(i, ptWidth, sigmaErr);

        gr_mean->SetPoint(i, ptCenter, h_ptResAbs[i]->GetMean());
        gr_mean->SetPointError(i, ptWidth, h_ptResAbs[i]->GetRMS()/sqrt(h_ptResAbs[i]->GetEntries()));
    }
    gr_res->SetName("gr_ptResolution");
    gr_mean->SetName("gr_ptMean");

    for (Long64_t i = 0; i < muonsGen->GetEntries(); ++i) {
        std::cout << "Processing gen muons entry " << i+1 << " of " << muonsGen->GetEntries() << "\r" << std::flush;
        muonsGen->GetEntry(i);
        pTMuonGenHist->Fill(pTMuonGen);
    }
    std::cout << std::endl;
    std::cout << "Entries in gen pT histogram: " << pTMuonGenHist->GetEntries() << std::endl;

    for (Long64_t i = 0; i < JPsiReco->GetEntries(); ++i) {
        std::cout << "Processing reco J/Psi entry " << i+1 << " of " << JPsiReco->GetEntries() << "\r" << std::flush;
        JPsiReco->GetEntry(i);
        if (massJPsiReco > signal_range_min && massJPsiReco < signal_range_max) {
            pTJPsiRecoHist->Fill(pTJPsiReco);
            pTJPsiRecoTrueHist->Fill(pTJPsiReco_true);
        } else if (massJPsiReco > background_range_min && massJPsiReco < background_range_max) {
            pTJPsiRecoHist->Fill(pTJPsiReco, -1.0);
            pTJPsiRecoTrueHist->Fill(pTJPsiReco_true, -1.0);
        }
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

    fillEfficiencHist(pTMuonRecoHist, pTMuonGenHist, muonEffHist);
    fillEfficiencHist(pTMuonRecoTrueHist, pTMuonGenHist, muonEffTrueHist);
    fillEfficiencHist(pTJPsiRecoHist, pTJPsiGenHist, JPsiEffHist);
    fillEfficiencHist(pTJPsiRecoTrueHist, pTJPsiGenHist, JPsiEffTrueHist);

    JPsiEffHist->Scale(1.0 / JPsi_branching_ratio);
    JPsiEffTrueHist->Scale(1.0 / JPsi_branching_ratio);

    // TEfficiency *muonEff = createEfficiencyGraph(pTMuonRecoHist, pTMuonGenHist);
    // TEfficiency *muonEffTrue = createEfficiencyGraph(pTMuonRecoTrueHist, pTMuonGenHist);
    // TEfficiency *JPsiEff = createEfficiencyGraph(pTJPsiRecoHist, pTJPsiGenHist);
    // TEfficiency *JPsiEffTrue = createEfficiencyGraph(pTJPsiRecoTrueHist, pTJPsiGenHist);

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

    // muonEff->Write();
    // muonEffTrue->Write();
    // JPsiEff->Write();
    // JPsiEffTrue->Write();

    pTMuonRecoTrueRecoHist->Write();
    gr_res->Write();
    gr_mean->Write();
    outFile->Write();
    outFile->Close();
}