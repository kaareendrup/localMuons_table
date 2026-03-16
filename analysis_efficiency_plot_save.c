#include "setALICEStyle.c"

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

    TString MC_name = "c3_global";
    // TString MC_name = "c3_standalone";
    // TString MC_name = "c3_global_temp";

    TString in_file = TString::Format("results/%s/efficiency.root", MC_name.Data());
    TFile* file = TFile::Open(in_file, "READ");

    TH1F *muonEffHist = (TH1F*)file->Get("muonEffHist");
    TH1F *JPsiEffHist = (TH1F*)file->Get("JPsiEffHist");

    TCanvas *c1 = new TCanvas("c1", "Muon and J/Psi eff", 1300, 600);
    c1->Divide(2,1);

    std::vector<float> pTCuts = {0, 20,.7, 20};
    std::vector<float> etaCuts = {-4, 4, -3.6, -2.5};
    gStyle->SetOptTitle(1);
    c1->cd(1);
    
    drawHist(muonEffHist, "#mu Efficiency;p_{T} (GeV/c);Eff_{#mu}", kBlue, 1.2);
    drawLabel_cuts(MC_name, "", &pTCuts, &etaCuts, 0.85, 0.59);

    c1->cd(2);
    drawHist(JPsiEffHist, "J/#Psi Efficiency;p_{T} (GeV/c);Eff_{J/#Psi}", kBlue, 1.2);
    drawLabel_cuts(MC_name, "", &pTCuts, &etaCuts, 0.85, 0.59);


    c1->SaveAs(TString::Format("results/%s/efficiency_plots.png", MC_name.Data()));
}