
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

void drawLabel(TString MC_name, TString type, Double_t x = 0.50, Double_t y = 0.85, TString extra_label = "") {

    TString label;
    std::vector<TString> details;

    if (type == "reco") {
        details.push_back("Reconstructed");
    } else if (type == "gen") {
        details.push_back("Generator level");
    } else if (type == "") {
        // Do nothing
    } else {
        details.push_back(type);
    }

    if (MC_name == "genpurp") {
        label = "Pythia General purpose MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("matchedQualityCuts");
    } else if (MC_name == "DQ") {
        label = "Pythia DQ prompt J/#Psi MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("matchedQualityCuts");
    } else if (MC_name == "HF") {
        label = "Pythia HF C#rightarrow#mu MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("matchedQualityCuts");
    } else if (MC_name == "k4h_baseline") {
        label = "Pythia General purpose MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("matchedQualityCuts");
    } else if (MC_name == "k4h_cuts") {
        label = "Pythia General purpose MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("matchedQualityCuts");
        details.push_back("#chi2 MFTMCH < 40");
        details.push_back("Ncluster MFC > 5");
        details.push_back("Ncluster MCH > 5");
        details.push_back("-3.6 < #eta < -2.5");
    } else if (MC_name == "k4h_standalone") {
        label = "Pythia General purpose MC";
        details.push_back("#sqrt{#it{s}} = 13.6 TeV");
        details.push_back("matchedMchMid");
        details.push_back("muonQualityCuts");
    } else {
        label = MC_name;
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAlign(31);
    latex.SetTextSize(0.04);
    latex.SetTextFont(42);
    latex.DrawLatex(x, y, label);
    for (size_t i = 0; i < details.size(); ++i) {
        latex.DrawLatex(x, y - 0.05 * (i + 1), details[i]);
    }
    if (extra_label != "") {
        latex.DrawLatex(x, y - 0.05 * (details.size() + 1), extra_label);
    }
}