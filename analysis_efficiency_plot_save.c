#include "setALICEStyle.c"
#include <nlohmann/json.hpp>
using json = nlohmann::json;

std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

void drawHist(TH1F* hist, TString title, int line_color, float scale_factor, bool same = false) {
    // Draw any histogram with consistent styling
    hist->SetLineColor(line_color);
    hist->SetLineWidth(2);
    hist->SetTitle(title);
    // hist->Scale(scale_factor);
    hist->SetMinimum(0);
    if (same) {
        hist->Draw("SAME");
    } else {
        hist->Draw();
    }
}

void analysis_efficiency_plot_save() {

    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    // Data
    std::string data_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    TString MC_name = TString::Format("%s_%s", data_name.c_str(), muon_type.c_str());

    TString in_file = TString::Format("results/%s/efficiency.root", MC_name.Data());
    TFile* file = TFile::Open(in_file, "READ");

    // Load metadata
    TTree* metaData = nullptr;
    file->GetObject("MetaData", metaData);
    std::vector<float>* pTCuts = nullptr;
    std::vector<float>* etaCuts = nullptr;
    metaData->SetBranchAddress("pTCuts", &pTCuts);
    metaData->SetBranchAddress("etaCuts", &etaCuts);
    metaData->GetEntry(0);

    // Efficiency plots
    TH1F *muonEffHist = (TH1F*)file->Get("muonEffHist");
    TH1F *muonEffTrueHist = (TH1F*)file->Get("muonEffTrueHist");
    TH1F *JPsiEffHist = (TH1F*)file->Get("JPsiEffHist");
    TH1F *JPsiEffTrueHist = (TH1F*)file->Get("JPsiEffTrueHist");

    TCanvas *c1 = new TCanvas("c1", "Muon and J/Psi eff", 1300, 600);
    c1->Divide(2,1);

    gStyle->SetOptTitle(1);
    c1->cd(1);
    
    drawHist(muonEffHist, "#mu Efficiency;p_{T} (GeV/c);Eff_{#mu}", kBlue, 1.2);
    drawLabel_cuts(MC_name, "", pTCuts, etaCuts, 0.85, 0.59);

    c1->cd(2);
    drawHist(JPsiEffHist, "J/#Psi Efficiency;p_{T} (GeV/c);Eff_{J/#Psi}", kBlue, 1.2);
    drawLabel_cuts(MC_name, "", pTCuts, etaCuts, 0.85, 0.59);

    c1->SaveAs(TString::Format("results/%s/efficiency_plots.png", MC_name.Data()));

    // Efficiency plots (true pT)
    TCanvas *c2 = new TCanvas("c2", "Muon and J/Psi eff", 1300, 600);
    c2->Divide(2,1);

    c2->cd(1);
    drawHist(muonEffTrueHist, "#mu Efficiency;p_{T} (GeV/c);Eff_{#mu}", kBlue, 1.2);
    drawLabel_cuts(MC_name, "", pTCuts, etaCuts, 0.85, 0.59);

    c2->cd(2);
    drawHist(JPsiEffTrueHist, "J/#Psi Efficiency;p_{T} (GeV/c);Eff_{J/#Psi}", kBlue, 1.2);
    drawLabel_cuts(MC_name, "", pTCuts, etaCuts, 0.85, 0.59);

    c2->SaveAs(TString::Format("results/%s/efficiency_plots_true.png", MC_name.Data()));

    // pT spectrum plots
    TH1F *pTMuonRecoHist = (TH1F*)file->Get("pTMuonRecoHist");
    TH1F *pTMuonRecoTrueHist = (TH1F*)file->Get("pTMuonRecoTrueHist");
    TH1F *pTJPsiRecoHist = (TH1F*)file->Get("pTJPsiRecoHist");
    TH1F *pTJPsiRecoTrueHist = (TH1F*)file->Get("pTJPsiRecoTrueHist");
    TH1F *pTMuonGenHist = (TH1F*)file->Get("pTMuonGenHist");
    TH1F *pTJPsiGenHist = (TH1F*)file->Get("pTJPsiGenHist");
    
    pTMuonRecoTrueHist->Scale(1, "width");
    pTMuonRecoHist->Scale(1, "width");
    pTJPsiRecoTrueHist->Scale(1, "width");
    pTJPsiRecoHist->Scale(1, "width");
    pTMuonGenHist->Scale(1, "width");
    pTJPsiGenHist->Scale(1, "width");
    
    TCanvas *c3 = new TCanvas("c3", "Muon and J/Psi p_{T}", 1300, 600);
    c3->Divide(2,1);
    
    c3->cd(1);
    drawHist(pTMuonRecoTrueHist, "#mu p_{T};p_{T} (GeV/c);Counts", kRed, 1.2);
    drawHist(pTMuonRecoHist, "#mu p_{T};p_{T} (GeV/c);Counts", kBlue, 1.2, true);
    TLegend *legend1 = new TLegend(0.7, 0.7, 0.85, 0.85);
    legend1->SetBorderSize(0);
    legend1->SetFillStyle(0);
    legend1->AddEntry(pTMuonRecoTrueHist, "true p_{T}", "l");
    legend1->AddEntry(pTMuonRecoHist, "reco p_{T}", "l");
    legend1->Draw();
    drawLabel_cuts(MC_name, "", pTCuts, etaCuts, 0.85, 0.59);
    
    c3->cd(2);
    drawHist(pTJPsiRecoTrueHist, "J/#Psi p_{T};p_{T} (GeV/c);Counts", kRed, 1.2);
    drawHist(pTJPsiRecoHist, "J/#Psi p_{T};p_{T} (GeV/c);Counts", kBlue, 1.2, true);
    TLegend *legend2 = new TLegend(0.7, 0.7, 0.85, 0.85);
    legend2->SetBorderSize(0);
    legend2->SetFillStyle(0);
    legend2->AddEntry(pTJPsiRecoTrueHist, "true p_{T}", "l");
    legend2->AddEntry(pTJPsiRecoHist, "reco p_{T}", "l");
    legend2->Draw();
    drawLabel_cuts(MC_name, "", pTCuts, etaCuts, 0.85, 0.59);
    
    c3->SaveAs(TString::Format("results/%s/pT_plots_true_reco.png", MC_name.Data()));
    
    TCanvas *c4 = new TCanvas("c4", "Muon and J/Psi p_{T}", 1300, 600);
    c4->Divide(2,1);
    
    c4->cd(1);
    drawHist(pTMuonGenHist, "#mu p_{T};p_{T} (GeV/c);Counts", kRed, 1.2);
    drawHist(pTMuonRecoHist, "#mu p_{T};p_{T} (GeV/c);Counts", kBlue, 1.2, true);
    TLegend *legend3 = new TLegend(0.7, 0.7, 0.85, 0.85);
    legend3->SetBorderSize(0);
    legend3->SetFillStyle(0);
    legend3->AddEntry(pTMuonGenHist, "gen p_{T}", "l");
    legend3->AddEntry(pTMuonRecoHist, "reco p_{T}", "l");
    legend3->Draw();
    drawLabel_cuts(MC_name, "", pTCuts, etaCuts, 0.85, 0.59);

    c4->cd(2);
    drawHist(pTJPsiGenHist, "J/#Psi p_{T};p_{T} (GeV/c);Counts", kRed, 1.2);
    drawHist(pTJPsiRecoHist, "J/#Psi p_{T};p_{T} (GeV/c);Counts", kBlue, 1.2, true);
    TLegend *legend4 = new TLegend(0.7, 0.7, 0.85, 0.85);
    legend4->SetBorderSize(0);
    legend4->SetFillStyle(0);
    legend4->AddEntry(pTJPsiGenHist, "gen p_{T}", "l");
    legend4->AddEntry(pTJPsiRecoHist, "reco p_{T}", "l");
    legend4->Draw();
    drawLabel_cuts(MC_name, "", pTCuts, etaCuts, 0.85, 0.59);

    c4->SaveAs(TString::Format("results/%s/pT_plots_gen_reco.png", MC_name.Data()));
}