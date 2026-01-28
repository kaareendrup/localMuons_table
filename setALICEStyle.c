
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

void drawLabel(TString MC_name, Double_t x = 0.50, Double_t y = 0.85, TString extra_label = "") {
    TString label, s1, s2;
    if (MC_name == "genpurp") {
        label = "Pythia General purpose MC";
        s1 = "#sqrt{#it{s}} = 5.36 TeV";
        s2 = "matchedQualityCuts";
    } else if (MC_name == "DQ") {
        label = "Pythia DQ prompt J/#Psi MC";
        s1 = "#sqrt{#it{s}} = 13.6 TeV";
        s2 = "matchedQualityCuts";
    } else if (MC_name == "HF") {
        label = "Pythia HF C#rightarrow#mu MC";
        s1 = "#sqrt{#it{s}} = 13.6 TeV";
        s2 = "matchedQualityCuts";
    } else {
        label = MC_name;
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextAlign(31);
    latex.SetTextSize(0.04);
    latex.SetTextFont(42);
    latex.DrawLatex(x, y, label);
    latex.DrawLatex(x, y - 0.05, s1);
    latex.DrawLatex(x, y - 0.10, s2);
    if (extra_label != "") {
        latex.DrawLatex(x, y - 0.15, extra_label);
    }
}