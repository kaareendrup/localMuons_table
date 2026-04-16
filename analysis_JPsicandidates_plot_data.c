
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

#include "utils/plots.c"

void analysis_JPsicandidates_plot_data() {

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
    float cuts_pT_JPsi_min = config["cuts_JPsi"]["pT_JPsi_min"];
    float cuts_pT_JPsi_max = config["cuts_JPsi"]["pT_JPsi_max"];
    float cuts_eta_JPsi_min = config["cuts_JPsi"]["eta_JPsi_min"];
    float cuts_eta_JPsi_max = config["cuts_JPsi"]["eta_JPsi_max"];
    std::vector<float> pT_JPsi_cuts = {cuts_pT_JPsi_min, cuts_pT_JPsi_max, 0, 20};
    std::vector<float> eta_JPsi_cuts = {cuts_eta_JPsi_min, cuts_eta_JPsi_max, -4.0, 4.0};

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
    
    // Import hepdata points
    std::string hepdata_name = config["hepdata_name"];
    double crossSection = config["hepdata_crosssection"];

    TString hepdata_file = TString::Format("results/%s.json", hepdata_name.c_str());
    std::ifstream hepdata_json(hepdata_file);
    json hepdata;
    hepdata_json >> hepdata;

    std::vector<double> hepdata_values;
    std::vector<double> hepdata_errors;
    for (const auto& point : hepdata["values"]) {
        double pt_low  = std::stod(point["x"][0]["low"].get<std::string>());
        double pt_high = std::stod(point["x"][0]["high"].get<std::string>());

        const auto& y = point["y"][0];
        double value = std::stod(y["value"].get<std::string>());

        double err_stat = 0.0;
        double err_sys  = 0.0;

        for (const auto& err : y["errors"]) {
            std::string label = err["label"];
            double e = std::stod(err["symerror"].get<std::string>());

            if (label == "stat") err_stat = e;
            if (label == "sys")  err_sys  = e;
        }

        if (pt_low >= cuts_pT_JPsi_min && pt_high <= cuts_pT_JPsi_max) { // Only include points within the plotted range
            hepdata_values.push_back(value);
            hepdata_errors.push_back(sqrt(pow(err_stat, 2) + pow(err_sys, 2)));
        }
    }

    // Create histograms
    TH1F* pT_reco = createInvMassHist("reco", config, data);
    
    // Import efficiency
    std::vector<double> efficiency = get_efficiency(config, MC_name);

    TCanvas *c3 = new TCanvas("c3", "pT bin counts", 900, 400);
    c3->Divide(2,1);
    
    // Create uncorrected histogram
    c3->cd(1);
    pT_reco->Draw();
    TLegend *leg3s = new TLegend(0.4,0.6,0.9,0.9);
    leg3s->AddEntry(pT_reco, "Reconstructed (not corrected)", "l");
    leg3s->SetBorderSize(0);
    leg3s->SetFillStyle(0);
    leg3s->Draw();

    // Create scaled histogram
    c3->cd(2);

    TH1F *pT_reco_scale = (TH1F*)pT_reco->Clone("pTscale");
    for (int i = 1; i <= pT_reco_scale->GetNbinsX(); i++) {
        // double bin_center = pT_reco_scale->GetBinCenter(i);
        // double deltaY = getDeltaY(bin_center, cuts_eta_JPsi_min, cuts_eta_JPsi_max);
        // std::cout << "Bin " << i << ": pT = " << bin_center << " GeV/c, Δy = " << deltaY << ", efficiency = " << efficiency[i-1] << std::endl;

        double w   = efficiency[i-1];//*deltaY;              // since ROOT bins start at 1
        double c   = pT_reco_scale->GetBinContent(i);
        double e   = pT_reco_scale->GetBinError(i);

        pT_reco_scale->SetBinContent(i, c / w);
        pT_reco_scale->SetBinError(i, e / w);              // scale uncertainties too
    }
    // Scale by number of events to get absolute yields
    pT_reco_scale->Scale(1./nEvents);
    pT_reco_scale->Draw("same");
    
    pT_reco_scale->GetYaxis()->SetTitle("d^{2}N/(dp_{T} dy) (GeV/c)^{-1}");

    setMax({pT_reco_scale});
    TLegend *leg3 = new TLegend(0.4,0.7,0.9,0.9);
    leg3->AddEntry(pT_reco_scale, "Reconstructed\n (corrected)", "l");
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->Draw();
    c3->SaveAs(TString::Format("results/%s/pTspectrascaled.png", data.Data()));

    // Create a canvas with two pads: top for the histogram, bottom for the ratio
    TCanvas *c4 = new TCanvas("c4", "pT bin counts", 700, 600);

    pT_reco_scale->Draw();
    
    // Add hepdata points with error bars
    TH1F *pTHepData = new TH1F("pTHepData", "pTHepData", pT_bins.size()-1, pT_bins.data());
    for (size_t i = 0; i < hepdata_values.size(); ++i) {
        pTHepData->SetBinContent(i+1, hepdata_values[i]);
        pTHepData->SetBinError(i+1, hepdata_errors[i]);
    }
    pTHepData->Scale(crossSection*1e-6); // Scale by cross section if needed
    // pTHepData->SetMarkerStyle(20);
    pTHepData->SetMarkerColor(kBlack);
    pTHepData->SetLineColor(kBlack);
    pTHepData->Draw("E1 SAME");
    
    pT_reco_scale->SetMinimum(1e-5);
    
    TLegend *leg4 = new TLegend(0.4,0.75,0.9,0.9);
    leg4->AddEntry(pT_reco_scale, "Reconstructed\n (corrected)", "l");
    leg4->AddEntry(pTHepData, "ALICE 2017", "lep");
    leg4->SetBorderSize(0);
    leg4->SetFillStyle(0);
    leg4->Draw();

    setMax({pT_reco_scale, pTHepData});
    gPad->SetLogy();
    increaseMargins(c4);
    drawLabel_cuts(data_name, "", &pT_JPsi_cuts, &eta_JPsi_cuts, 0.45, 0.55);
    c4->SaveAs(TString::Format("results/%s/pTspectracompare.png", data.Data()));
}