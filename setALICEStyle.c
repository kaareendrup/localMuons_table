
#include "TLatex.h"

void SetALICEStyle() {
    gStyle->SetOptStat(0);          // No stat box
    gStyle->SetOptTitle(0);         // No default title

    gStyle->SetPadTickX(1);         // Ticks on both sides
    gStyle->SetPadTickY(1);

    gStyle->SetFrameLineWidth(2);
    gStyle->SetLineWidth(2);

    gStyle->SetLabelSize(0.05, "XY");
    gStyle->SetTitleSize(0.06, "XY");
    gStyle->SetTitleOffset(1.2, "X");
    gStyle->SetTitleOffset(1.3, "Y");

    gStyle->SetLegendBorderSize(0);
    gStyle->SetLegendFont(42);

    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XY");
    gStyle->SetTitleFont(42, "XY");
}

void increaseMargins(TCanvas* c) {
    c->SetLeftMargin(0.15);
    c->SetBottomMargin(0.15);
    c->SetRightMargin(0.05);
    c->SetTopMargin(0.05);
}

void increasePadMargins(TCanvas* c, int n) {
    for (int i = 1; i <= n; ++i) {
        c->cd(i);
        gPad->SetLeftMargin(0.1); 
        gPad->SetRightMargin(0.05); 
        gPad->SetBottomMargin(0.18); 
        gPad->SetTopMargin(0.05);
    }
}

void drawLabel_cuts(
    TString MC_name, 
    TString type, 
    std::vector<float> *pTCuts = nullptr, 
    std::vector<float> *etaCuts = nullptr,
    Double_t x = 0.50, 
    Double_t y = 0.85, 
    bool ralign = true
) {

    TString label;
    std::vector<TString> details;

    // Add reco/gen to detailed info
    if (type == "reco") {
        details.push_back("Reconstructed");
    } else if (type == "gen") {
        details.push_back("Generator level");
    } else if (type == "") {
        // Do nothing
    } else {
        details.push_back(type);
    }

    // Add dataset info
    if (MC_name == "DQ") {
        label = "Pythia DQ prompt J/#Psi MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        // details.push_back("matchedQualityCuts");
        details.push_back("matchedMchMid");
        details.push_back("muonQualityCuts");
    } else if (MC_name == "k4h_baseline") {
        label = "Pythia General purpose MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("matchedQualityCuts");
    } else if (MC_name == "k4h_standalone") {
        label = "Pythia General purpose MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("matchedMchMid");
        details.push_back("muonQualityCuts");
    } else if (MC_name == "DQ_data") {
        label = "2024 pp data";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("matchedMchMid");
        details.push_back("muonQualityCuts");
    } else if (MC_name == "f4d_global") {
        label = "Pythia General purpose MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("muonQualityCutsGlobal");
    } else if (MC_name == "f4d_standalone") {
        label = "Pythia General purpose MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("muonQualityCutsStandalone");
    } else if (MC_name == "c3_global") {
        label = "Pythia DQ prompt J/#Psi MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("muonQualityCutsGlobal");
    } else if (MC_name == "c3_standalone") {
        label = "Pythia DQ prompt J/#Psi MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("muonQualityCutsStandalone");
    } else {
        label = MC_name;
    }

    int cuts_added = 0;
    if (pTCuts && etaCuts) {
        // Add cut info
        // Cuts are specified as:
        // {pT_trigger_min, pT_trigger_max, pT_assoc_min, pT_assoc_max}
        // {eta_trigger_min, eta_trigger_max, eta_assoc_min, eta_assoc_max}

        // Loop over mu and J/Psi
        for (int i = 0; i < 2; ++i) {
            // TString particle = (i == 0) ? "trig" : "assoc";
            TString particle = (i == 0) ? "trig" : "#mu";
            TString pT_str = Form("p_{T,%s}", particle.Data());
            TString eta_str = Form("#eta_{%s}", particle.Data());

            // pT cuts
            if (pTCuts->at(2*i+0) > 0) {
                pT_str = Form("%.1f < %s", pTCuts->at(2*i+0), pT_str.Data());
            }
            if (pTCuts->at(2*i+1) < 20) {
                pT_str = Form("%s < %.1f", pT_str.Data(), pTCuts->at(2*i+1));
            }

            // Eta cuts
            if (etaCuts->at(2*i+0) > -4.0) {
                eta_str = Form("%.1f < %s", etaCuts->at(2*i+0), eta_str.Data());
            }
            if (etaCuts->at(2*i+1) < 4.0) {
                eta_str = Form("%s < %.1f", eta_str.Data(), etaCuts->at(2*i+1));
            }

            // Check if cuts are non-trivial before adding to details
            if (pT_str != Form("p_{T,%s}", particle.Data())) {
                details.push_back(pT_str + " GeV/c");
                cuts_added++;
            }
            if (eta_str != Form("#eta_{%s}", particle.Data())) {
                details.push_back(eta_str);
                cuts_added++;
            }
        }
    }

    TLatex latex;
    latex.SetNDC();
    if (ralign) {
        latex.SetTextAlign(31);
    } else {
        latex.SetTextAlign(11);
    }
    latex.SetTextSize(0.035);
    latex.SetTextFont(42);
    latex.DrawLatex(x, y, label);
    for (size_t i = 0; i < details.size(); ++i) {
        float y_offset = 0.05 * (i + 1);
        if (i > details.size() - cuts_added) {
            // Add extra spacing between cut details
            y_offset = 0.05 * (details.size() - cuts_added + 1) + 0.06 * (i - (details.size() - cuts_added)); 
        }
        latex.DrawLatex(x, y - y_offset, details[i]);
    }
}