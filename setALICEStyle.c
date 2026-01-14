
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

void drawLabel(TString MC_name, Double_t x = 0.18, Double_t y = 0.85, TString extra_label = "") {
    TString label;
    if (MC_name == "genpurp") {
        label = "General purpose MC";
    } else if (MC_name == "DQ") {
        label = "DQ prompt J/#Psi MC";
    } else if (MC_name == "HF") {
        label = "HF C#rightarrow#mu MC";
    } else {
        label = MC_name;
    }

    TLatex latex;
    latex.SetNDC();
    latex.SetTextSize(0.05);
    latex.SetTextFont(42);
    latex.DrawLatex(x, y, label + " " + extra_label);
}