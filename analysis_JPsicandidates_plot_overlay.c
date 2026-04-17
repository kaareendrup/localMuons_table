
#include <TFile.h>
#include <TTree.h>
#include <TKey.h>
#include <TDirectory.h>
#include <TString.h>
#include <TMath.h>
#include <Math/Vector4D.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TLegend.h>
#include <TStyle.h>

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <iomanip>

#include <nlohmann/json.hpp>
using json = nlohmann::json;

#include "effUtils.c"

void scale_histogram(TH1F* hist, const std::vector<double>& efficiency, double cuts_eta_JPsi_min, double cuts_eta_JPsi_max, int nEvents) {

    for (int i = 1; i <= hist->GetNbinsX(); i++) {

        // double bin_center = pT_reco_scale->GetBinCenter(i);
        // double deltaY = getDeltaY(bin_center, cuts_eta_JPsi_min, cuts_eta_JPsi_max);
        // std::cout << "Bin " << i << ": pT = " << bin_center << " GeV/c, Δy = " << deltaY << ", efficiency = " << efficiency[i-1] << std::endl;

        double w   = efficiency[i-1];//*deltaY;              // since ROOT bins start at 1
        double c   = hist->GetBinContent(i);
        double e   = hist->GetBinError(i);

        hist->SetBinContent(i, c / w);
        hist->SetBinError(i, e / w);              // scale uncertainties too
    }
    // Scale by number of events to get absolute yields
    hist->Scale(1./nEvents);
}

void analysis_JPsicandidates_plot_overlay() {
    
    std::cout << std::fixed << std::setprecision(1);
    SetALICEStyle();

    // Load config
    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    std::string data_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    std::string eff_source = config["data_eff_source"];
    std::vector<double> pT_bins = config["hists"]["pT_bins"].get<std::vector<double>>();
    double cuts_pT_JPsi_min = config["cuts_JPsi"]["pT_JPsi_min"];
    double cuts_pT_JPsi_max = config["cuts_JPsi"]["pT_JPsi_max"];
    double cuts_eta_JPsi_min = config["cuts_JPsi"]["eta_JPsi_min"];
    double cuts_eta_JPsi_max = config["cuts_JPsi"]["eta_JPsi_max"];

    TString data = TString::Format("%s_%s", data_name.c_str(), muon_type.c_str());
    TString MC_name = TString::Format("%s_%s", eff_source.c_str(), muon_type.c_str());

    // Import metadata
    TString in_file = TString::Format("results/%s/reco/invMassSpektra.root", data.Data());
    TFile* file = TFile::Open(in_file, "READ");
    TTree *metaTree = nullptr;
    file->GetObject("MetaData", metaTree);

    int nEvents;
    metaTree->SetBranchAddress("nEvents", &nEvents);
    metaTree->GetEntry(0);
    file->Close();

    TString in_file_MC = TString::Format("results/%s/reco/invMassSpektra.root", MC_name.Data());
    TFile* file_MC = TFile::Open(in_file_MC, "READ");
    TTree *metaTree_MC = nullptr;
    file_MC->GetObject("MetaData", metaTree_MC);

    int nEvents_MC;
    metaTree_MC->SetBranchAddress("nEvents", &nEvents_MC);
    metaTree_MC->GetEntry(0);
    file_MC->Close();

    // Import efficiency
    std::vector<double> efficiency = get_efficiency(config, MC_name);

    // Draw histograms overlayed
    TCanvas *c3 = new TCanvas("c3", "pT bin counts", 900, 400);
    c3->Divide(2,1);
    
    // Create histograms
    TH1F* pT_reco_data = createInvMassHist("reco", config, data);
    TH1F* pT_reco_MC = createInvMassHist("reco", config, MC_name);
    
    // Create uncorrected histogram
    c3->cd(1);
    pT_reco_data->Draw();
    pT_reco_MC->SetLineColor(kBlue);
    pT_reco_MC->Draw("same");
    setMax({pT_reco_data, pT_reco_MC});
    TLegend *leg3s = new TLegend(0.4,0.6,0.9,0.9);
    leg3s->AddEntry(pT_reco_data, "Data\n (not corrected)", "l");
    leg3s->AddEntry(pT_reco_MC, "MC\n (not corrected)", "l");
    leg3s->Draw();

    // Create scaled histogram
    c3->cd(2);
    
    TH1F *pT_reco_scale = (TH1F*)pT_reco_data->Clone("pTscale");
    TH1F *pT_reco_scale_MC = (TH1F*)pT_reco_MC->Clone("pTscaleMC");
    scale_histogram(pT_reco_scale, efficiency, cuts_eta_JPsi_min, cuts_eta_JPsi_max, nEvents);
    scale_histogram(pT_reco_scale_MC, efficiency, cuts_eta_JPsi_min, cuts_eta_JPsi_max, nEvents_MC);
    pT_reco_scale->Draw();
    pT_reco_scale_MC->Draw("same");

    setMax({pT_reco_scale, pT_reco_scale_MC});
    TLegend *leg3 = new TLegend(0.4,0.7,0.9,0.9);
    leg3->AddEntry(pT_reco_scale, "Data\n (corrected)", "l");
    leg3->AddEntry(pT_reco_scale_MC, "MC\n (corrected)", "l");
    leg3->Draw();
    c3->SaveAs(TString::Format("results/%s/pTspectrascaled_overlay.png", data.Data()));
}