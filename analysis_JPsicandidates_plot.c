
#include "setALICEStyle.c"
#include <nlohmann/json.hpp>
using json = nlohmann::json;

std::cout << std::fixed << std::setprecision(1);
SetALICEStyle();

void setMax(std::vector<TH1F*> hists) {
    // Adjust y-axis maximum to be 1.5 times the largest maximum among the provided histograms
    double max_val = 0;
    for (auto hist : hists) {
        if (hist->GetMaximum() > max_val) {
            max_val = hist->GetMaximum();
        }
    }
    for (auto hist : hists) {
        hist->SetMaximum(1.5 * max_val);
    }
}

std::vector<double> get_efficiency(json config, TString MC_name) {

    // Data
    TString in_file = TString::Format("results/%s/efficiency.root", MC_name.Data());
    TFile* file = TFile::Open(in_file, "READ");

    // Efficiency plots
    TH1F *JPsiEffHist = (TH1F*)file->Get("JPsiEffHist");

    int nbins = config["hists"]["pT_bins"].size();
    std::vector<double> efficiency(nbins);
    for (int i = 0; i < nbins; ++i) {
        efficiency[i] = JPsiEffHist->GetBinContent(i+1);
    }
    return efficiency;
}

TH1F *createInvMassHist(TString type, json config, TString MC_name) {

    // Load json config
    std::vector<double> pT_bins = config["hists"]["pT_bins"].get<std::vector<double>>();
    
    float signal_range_min = config["signal_range"]["min"];
    float signal_range_max = config["signal_range"]["max"];
    float background_range_min = config["background_range"]["min"];
    float background_range_max = config["background_range"]["max"];
    
    // Load data
    TString in_file = TString::Format("results/%s/%s/invMassSpektra.root", MC_name.Data(), type.Data());
    TFile* file = TFile::Open(in_file, "READ");
    
    TH1F *pT_sig = (TH1F*)file->Get("pT_sig");
    TH1F *pT_bkg_low = (TH1F*)file->Get("pT_bkg_low");
    TH1F *pT_bkg_high = (TH1F*)file->Get("pT_bkg_high");

    int n_bins = pT_bins.size() - 1;
    int n_cols = std::ceil(std::sqrt(n_bins));
    int n_rows = std::ceil((double)n_bins / n_cols);
    int canvas_width = n_cols * 400; // 400 pixels per plot
    int canvas_height = n_rows * 400; // 400 pixels per plot

    TCanvas *c1 = new TCanvas("c1", "Invariant Mass Spectra by pT", canvas_width, canvas_height);
    c1->Divide(n_cols, n_rows);
    
    gStyle->SetOptTitle(1);
    for (int i = 0; i < n_bins; ++i) {
        
        c1->cd(i+1);
        TString category_str = TString::Format("pT_%.1f_%.1f_invMass", pT_bins[i], pT_bins[i+1]);
        TH1F *invMass_pt = (TH1F*)file->Get(category_str);

        invMass_pt->SetLineWidth(2);
        invMass_pt->SetLineColor(kRed);
        invMass_pt->SetTitle(category_str);
        invMass_pt->GetXaxis()->SetTitle("Invariant Mass (GeV/c^{2})");
        invMass_pt->GetYaxis()->SetTitle("Counts");
        // invMass_pt->Sumw2();
        invMass_pt->Draw();

        double yMinSingle = invMass_pt->GetMinimum();
        double yMaxSingle = invMass_pt->GetMaximum();

        // Draw shaded boxes for background regions
        TBox* box1s = new TBox(background_range_min, yMinSingle, signal_range_min, yMaxSingle);
        box1s->SetFillColorAlpha(kGray, 1); // semi-transparent
        box1s->SetFillStyle(3004);
        box1s->Draw("same");

        TBox* box2s = new TBox(signal_range_max, yMinSingle, background_range_max, yMaxSingle);
        box2s->SetFillColorAlpha(kGray, 1); // semi-transparent
        box2s->SetFillStyle(3004);
        box2s->Draw("same");    
    }
    gStyle->SetOptTitle(0);

    c1->SaveAs(TString::Format("results/%s/invMasspTBins_%s.png", MC_name.Data(), type.Data()));

    TCanvas *c2 = new TCanvas("c2", "pT bin counts", 450, 400);
    // pT_sig->Sumw2();
    pT_sig->Scale(1, "width");
    pT_sig->SetLineWidth(2);
    pT_sig->SetLineColor(kRed);
    pT_sig->SetMinimum(1e1); // Set minimum to 1 for log scale
    pT_sig->Draw();
    // pT_bkg_low->Sumw2();
    pT_bkg_low->Scale(1, "width");
    pT_bkg_low->SetLineWidth(2);
    pT_bkg_low->SetLineColor(kBlue);
    pT_bkg_low->Draw("same");
    // pT_bkg_high->Sumw2();
    pT_bkg_high->Scale(1, "width");
    pT_bkg_high->SetLineWidth(2);
    pT_bkg_high->SetLineColor(kBlue+1);
    pT_bkg_high->Draw("same");
    
    TLegend *leg2 = new TLegend(0.6,0.7,0.9,0.9);
    leg2->AddEntry(pT_sig, "Signal region", "l");
    leg2->AddEntry(pT_bkg_low, "Lower background", "l");
    leg2->AddEntry(pT_bkg_high, "Higher background", "l");
    leg2->SetBorderSize(0);
    leg2->SetFillStyle(0);
    leg2->Draw();
    
    c2->SetLogy();
    c2->SaveAs(TString::Format("results/%s/pTspectra_%s.png", MC_name.Data(), type.Data()));
    
    TH1F *pT_sig_sub = (TH1F*)pT_sig->Clone("pTsub");
    pT_sig_sub->Add(pT_bkg_low, -1);
    pT_sig_sub->Add(pT_bkg_high, -1);
    return pT_sig_sub;

}
    
void analysis_JPsicandidates_plot() {
        
    std::ifstream jsonFile("localMuons_table/config/config_analysis.json");
    json config;
    jsonFile >> config;

    std::string data_name = config["data_name"];
    std::string muon_type = config["muon_type"];
    TString MC_name = TString::Format("%s_%s", data_name.c_str(), muon_type.c_str());

    TH1F* pT_reco = createInvMassHist("reco", config, MC_name);
    TH1F* pT_gen = createInvMassHist("gen", config, MC_name);

    std::vector<double> efficiency = get_efficiency(config, MC_name);

    TCanvas *c3 = new TCanvas("c3", "pT bin counts", 900, 400);
    c3->Divide(2,1);
    
    c3->cd(1);
    pT_reco->Draw();
    TLegend *leg3s = new TLegend(0.4,0.6,0.9,0.9);
    leg3s->AddEntry(pT_reco, "Reconstructed (not corrected)", "l");
    leg3s->SetBorderSize(0);
    leg3s->SetFillStyle(0);
    leg3s->Draw();

    c3->cd(2);
    pT_gen->SetLineColor(kBlue);
    pT_gen->Draw();

    TH1F *pT_reco_scale = (TH1F*)pT_reco->Clone("pTscale");
    for (int i = 1; i <= pT_reco_scale->GetNbinsX(); i++) {
        double w   = efficiency[i-1];              // since ROOT bins start at 1
        double c   = pT_reco_scale->GetBinContent(i);
        double e   = pT_reco_scale->GetBinError(i);

        pT_reco_scale->SetBinContent(i, c / w);
        pT_reco_scale->SetBinError(i, e / w);              // scale uncertainties too
    }
    pT_reco_scale->Draw("same");
    
    setMax({pT_gen, pT_reco_scale});
    TLegend *leg3 = new TLegend(0.4,0.5,0.9,0.9);
    leg3->AddEntry(pT_gen, "Generated", "l");
    leg3->AddEntry(pT_reco_scale, "Reconstructed\n (corrected)", "l");
    leg3->SetBorderSize(0);
    leg3->SetFillStyle(0);
    leg3->Draw();
    c3->SaveAs(TString::Format("results/%s/pTspectrascaled.png", MC_name.Data()));

    TCanvas *c4 = new TCanvas("c4", "pT bin counts", 600, 600);
    c4->Divide(1,2);
  
    c4->cd(1);
    gPad->SetPad(0,0.3,1,1);    // Top 70%
    pT_gen->Draw();
    pT_reco_scale->Draw("same");
    gPad->SetBottomMargin(0); // Remove bottom margin for top pad
    
    c4->cd(2);
    gPad->SetPad(0,0,1,0.3);    // Bottom 30%
    gPad->SetTopMargin(0); // Remove top margin for bottom pad
    TH1F *ratio_hist = (TH1F*)pT_gen->Clone("ratio_hist");
    ratio_hist->Divide(pT_reco_scale);
    ratio_hist->Draw();
    gPad->SetGridy();
    ratio_hist->SetMinimum(.6);
    ratio_hist->SetMaximum(2.3);
    ratio_hist->GetYaxis()->SetTitleOffset(0.5);
    ratio_hist->GetYaxis()->SetTitleSize(0.09);
    ratio_hist->GetYaxis()->SetLabelSize(0.08);
    ratio_hist->GetYaxis()->SetNdivisions(505);
    ratio_hist->GetXaxis()->SetTitleSize(0.09);
    ratio_hist->GetXaxis()->SetLabelSize(0.08);

    gPad->SetBottomMargin(0.3);
    c4->SetBottomMargin(0.8);
    c4->SaveAs(TString::Format("results/%s/pTspectraratio.png", MC_name.Data()));
}