
#include "setALICEStyle.c"

SetALICEStyle();

TTree* get_tree(TKey *key, TFile *file) {

    const char* dirName = key->GetName();
    TDirectoryFile *dir = (TDirectoryFile*)file->Get(dirName);
    TTree *tree = (TTree*)dir->Get("O2dqmuontable");

    return tree;
}

void plot_source_matrix(TString data_file, TString MC_name) {

    TFile *file = TFile::Open(data_file);
    auto *keys = file->GetListOfKeys();
    const TString sources[6] = {"fIsPrimary", "fIsProducedInTransport", "fIsProducedByGenerator", "fIsFromBackgroundEvent", "fIsHEPMCFinalState", "fIsPowhegDY"};
    bool flags[6];

    TH2I *sourceMatrix = new TH2I("sourceMatrix", "Source Matrix; ; ", 6, 0, 6, 6, 0, 6);

    // Loop over dataframes
    for (int i = 0; i < keys->GetEntries()-1; ++i) {

        std::cout << "Processing " << i << "/" << keys->GetEntries()-1 << std::endl;

        TTree *tree = get_tree((TKey*)keys->At(i), file);

        // Set branch addresses for source flags
        for (int i =0; i < 6; ++i){
            tree->SetBranchAddress(sources[i], &flags[i]);
        }

        // Loop over entries
        Long64_t n = tree->GetEntries();
        for (Long64_t j = 0; j < n; ++j) {
            tree->GetEntry(j);

            // Fill source matrix
            for (int m = 0; m < 6; ++m){
                if (!flags[m]) continue;

                for (int n = m; n < 6; ++n){
                    if (flags[m] && flags[n]){
                        sourceMatrix->Fill(m, n);
                    }
                }
            }
        }
    }

    for (int i = 0; i < 6; ++i){
        sourceMatrix->GetXaxis()->SetBinLabel(i+1, sources[i]);
        sourceMatrix->GetYaxis()->SetBinLabel(i+1, sources[i]);
    }

    TCanvas *c1 = new TCanvas("c1", "Source Matrix", 900, 600);
    gStyle->SetOptStat(0);
    sourceMatrix->Draw("COLZ TEXT");

    c1->SetLeftMargin(0.25);
    c1->SetBottomMargin(0.15);
    c1->SetRightMargin(0.15);
    c1->SetTopMargin(0.15);

    c1->SaveAs("results/" + MC_name + "/source_matrix.png");
    delete c1;
}

void plot_source_pT(TString data_file, TString MC_name) {

    TFile *file = TFile::Open(data_file);
    auto *keys = file->GetListOfKeys();

    const TString sources[2] = {"fIsPrimary", "fIsProducedInTransport"};
    bool flags[2];
    std::vector<Double_t> pT_Primary, pT_Transport, pT_Other;
    float pT;

    // Loop over dataframes
    for (int i = 0; i < keys->GetEntries()-1; ++i) {

        std::cout << "Processing " << i << "/" << keys->GetEntries()-1 << std::endl;

        TTree *tree = get_tree((TKey*)keys->At(i), file);

        // Set branch addresses for source flags and pT
        tree->SetBranchAddress("fPtassoc", &pT);
        tree->SetBranchAddress("fIsPrimary", &flags[0]);
        tree->SetBranchAddress("fIsProducedInTransport", &flags[1]);

        // Loop over entries
        for (Long64_t j = 0; j < tree->GetEntries(); ++j) {
            tree->GetEntry(j); 

            if (flags[0]) {
                pT_Primary.push_back(pT);
            }
            else if (flags[1]) {
                pT_Transport.push_back(pT);
            }
            else {
                pT_Other.push_back(pT);
            }
        }
    }

    TCanvas *c1 = new TCanvas("c1", "Source pT Distributions", 1200, 600);
    TH1F *pTPrimaryHist = new TH1F("h1","Primary Muon pT Distribution;p_{T} (GeV/c);Counts",100,0,10);
    pTPrimaryHist->FillN(pT_Primary.size(), pT_Primary.data(), nullptr);
    pTPrimaryHist->SetLineColor(kBlue+1);
    pTPrimaryHist->SetLineWidth(2);
    pTPrimaryHist->Draw();

    TH1F *pTTransportHist = new TH1F("h2","Transport Muon pT Distribution;p_{T} (GeV/c);Counts",100,0,10);
    pTTransportHist->FillN(pT_Transport.size(), pT_Transport.data(), nullptr);
    pTTransportHist->SetLineColor(kRed+1);
    pTTransportHist->SetLineWidth(2);
    pTTransportHist->Draw("SAME");

    TH1F *pTOtherHist = new TH1F("h3","Other Muon pT Distribution;p_{T} (GeV/c);Counts",100,0,10);
    pTOtherHist->FillN(pT_Other.size(), pT_Other.data(), nullptr);
    pTOtherHist->SetLineColor(kBlack);
    pTOtherHist->SetLineWidth(2);
    pTOtherHist->Draw("SAME");

    increaseMargins(c1);
    drawLabel(MC_name);
    TLegend *legend = new TLegend(0.6,0.7,0.88,0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(pTPrimaryHist, "Primary Muons", "l");
    legend->AddEntry(pTTransportHist, "Transport Muons", "l");
    legend->AddEntry(pTOtherHist, "Other Muons", "l");
    legend->Draw();

    // Adjust ymax 
    double maxY = std::max(pTPrimaryHist->GetMaximum(), pTTransportHist->GetMaximum());
    pTPrimaryHist->SetMaximum(1.2 * maxY);
    pTTransportHist->SetMaximum(1.2 * maxY);

    // Set log Y axis
    gPad->SetLogy();

    c1->SaveAs("results/" + MC_name + "/source_pT_distributions.png");
    delete c1;
}

void plot_source_PID(TString data_file, TString MC_name, TString source_flag) {

    int n_bins = 20;
    float range_min = 0;
    float range_max = 10;

    // Load the dataframe keys
    TFile *file = TFile::Open(data_file);
    auto *keys = file->GetListOfKeys();

    // Initialize inv mass variables
    std::vector<Double_t> all_pT, JPsi_pT, Psi2S_pT, charm_pT, b_pT, K_pT, Pi_pT, LM_pT, noMC_pT, other_pT;

    // Loop over dataframes
    for (int i = 0; i < keys->GetEntries()-1; ++i) {

        std::cout << "Processing " << i << "/" << keys->GetEntries()-1 << std::endl;

        TTree *tree = get_tree((TKey*)keys->At(i), file);
                
        // Prepare to read muon info
        Long64_t fMotherPDG;
        float fPt;
        bool fSourceFlag;
        // tree->SetBranchAddress("fMotherPDG", &fMotherPDG);
        tree->SetBranchAddress("fGrandmotherPDG", &fMotherPDG);
        tree->SetBranchAddress("fPtassoc", &fPt);
        tree->SetBranchAddress(source_flag, &fSourceFlag);
        
        // First pass: build groups of muons from the same event
        Long64_t n = tree->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            tree->GetEntry(i);

            if (!fSourceFlag) continue;
            
            all_pT.push_back(fPt);

            if (std::abs(fMotherPDG) == 443) {
                JPsi_pT.push_back(fPt);
            } else if (std::abs(fMotherPDG) == 100443) {
                Psi2S_pT.push_back(fPt);
            } else if ((std::abs(fMotherPDG) >= 411 && std::abs(fMotherPDG) <= 445) || 
                       (std::abs(fMotherPDG) >= 4101 && std::abs(fMotherPDG) <= 4444)) {
                charm_pT.push_back(fPt);
            } else if ((std::abs(fMotherPDG) >= 511 && std::abs(fMotherPDG) <= 557) || 
                       (std::abs(fMotherPDG) >= 5101 && std::abs(fMotherPDG) <= 5554)) {
                b_pT.push_back(fPt);
            } else if (std::abs(fMotherPDG) == 321 || std::abs(fMotherPDG) == 311 || 
                       std::abs(fMotherPDG) == 130 || std::abs(fMotherPDG) == 323 ) {
                K_pT.push_back(fPt);
            } else if (std::abs(fMotherPDG) == 211 || std::abs(fMotherPDG) == 113 || 
                       std::abs(fMotherPDG) == 111 || std::abs(fMotherPDG) == 213 ) {
                Pi_pT.push_back(fPt);
            } else if (std::abs(fMotherPDG) == 221 || std::abs(fMotherPDG) == 331 || 
                       std::abs(fMotherPDG) == 223 || std::abs(fMotherPDG) == 333 ) {
                LM_pT.push_back(fPt);
            } else if (fMotherPDG == -9999) {
                noMC_pT.push_back(fPt);
            } else {
                std::cout << "Found unassigned PID:" << fMotherPDG << std::endl;
                other_pT.push_back(fPt);
            }
        }

        tree->ResetBranchAddresses();
    }

    // Plot histogram of pT distributions
    TCanvas *c1 = new TCanvas("c1", "pT of muons", 800, 600);
    TH1F *pTHist = new TH1F("h1","pT of Muons;pT (GeV/c);Counts",n_bins,range_min,range_max);
    pTHist->FillN(all_pT.size(), all_pT.data(), nullptr);
    pTHist->SetLineWidth(3); 
    pTHist->SetLineColor(kBlack); 
    pTHist->Draw();

    // Add secondary histograms
    TH1F *jpsiHist = new TH1F("h2","",n_bins,range_min,range_max);
    jpsiHist->FillN(JPsi_pT.size(), JPsi_pT.data(), nullptr);
    jpsiHist->SetLineColor(kGreen+1);
    jpsiHist->SetLineWidth(2);
    jpsiHist->Sumw2();
    jpsiHist->Draw("HIST SAME");
    jpsiHist->Draw("E SAME");

    TH1F *psi2sHist = new TH1F("h3","",n_bins,range_min,range_max);
    psi2sHist->FillN(Psi2S_pT.size(), Psi2S_pT.data(), nullptr);
    psi2sHist->SetLineColor(kGreen+2);
    psi2sHist->SetLineWidth(2);
    psi2sHist->Sumw2();
    psi2sHist->Draw("HIST SAME");
    psi2sHist->Draw("E SAME");

    TH1F *charmHist = new TH1F("h4","",n_bins,range_min,range_max);
    charmHist->FillN(charm_pT.size(), charm_pT.data(), nullptr);
    charmHist->SetLineColor(kMagenta);
    charmHist->SetLineWidth(2);
    charmHist->Sumw2();
    charmHist->Draw("HIST SAME");
    charmHist->Draw("E SAME");

    TH1F *bHist = new TH1F("h5","",n_bins,range_min,range_max);
    bHist->FillN(b_pT.size(), b_pT.data(), nullptr);
    bHist->SetLineColor(kRed);
    bHist->SetLineWidth(2);
    bHist->Sumw2();
    bHist->Draw("HIST SAME");
    bHist->Draw("E SAME");

    TH1F *KHist = new TH1F("h6","",n_bins,range_min,range_max);
    KHist->FillN(K_pT.size(), K_pT.data(), nullptr);
    KHist->SetLineColor(kBlue);
    KHist->SetLineWidth(2);
    KHist->Sumw2();
    KHist->Draw("HIST SAME");
    KHist->Draw("E SAME");

    TH1F *PiHist = new TH1F("h7","",n_bins,range_min,range_max);
    PiHist->FillN(Pi_pT.size(), Pi_pT.data(), nullptr);
    PiHist->SetLineColor(kOrange);
    PiHist->SetLineWidth(2);
    PiHist->Sumw2();
    PiHist->Draw("HIST SAME");
    PiHist->Draw("E SAME");

    TH1F *LMHist = new TH1F("h8","",n_bins,range_min,range_max);
    LMHist->FillN(LM_pT.size(), LM_pT.data(), nullptr);
    LMHist->SetLineColor(kOrange+2);
    LMHist->SetLineWidth(2);
    LMHist->Sumw2();
    LMHist->Draw("HIST SAME");
    LMHist->Draw("E SAME");

    TH1F *noMCHist = new TH1F("h9","",n_bins,range_min,range_max);
    noMCHist->FillN(noMC_pT.size(), noMC_pT.data(), nullptr);
    noMCHist->SetLineColor(kGray);
    noMCHist->SetLineWidth(2);
    noMCHist->Sumw2();
    noMCHist->Draw("HIST SAME");
    noMCHist->Draw("E SAME");

    TH1F *otherHist = new TH1F("h10","",n_bins,range_min,range_max);
    otherHist->FillN(other_pT.size(), other_pT.data(), nullptr);
    otherHist->SetLineColor(kCyan);
    otherHist->SetLineWidth(2);
    otherHist->Sumw2();
    otherHist->Draw("HIST SAME");
    otherHist->Draw("E SAME");
    
    // Style
    gPad->SetLogy();
    increaseMargins(c1);

    // Legend
    TLegend *legend = new TLegend(0.58,0.6,0.92,0.88);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->AddEntry(pTHist, "All #mu", "l");
    legend->AddEntry(jpsiHist, "#mu from J/#psi", "l");
    legend->AddEntry(psi2sHist, "#mu from #psi(2S)", "l");
    legend->AddEntry(charmHist, "#mu from Charm", "l");
    legend->AddEntry(bHist, "#mu from Beauty", "l");
    legend->AddEntry(KHist, "#mu from Kaons", "l");
    legend->AddEntry(PiHist, "#mu from Pions", "l");
    legend->AddEntry(LMHist, "#mu from other light mesons (#eta,#omega,#phi)", "l");
    legend->AddEntry(otherHist, "Other #mu", "l");
    legend->AddEntry(noMCHist, "No mother", "l");
    legend->Draw();

    // Adjust ymax
    double maxY = pTHist->GetMaximum();
    pTHist->SetMaximum(4 * maxY);

    drawLabel(MC_name, 0.57, 0.89, source_flag);
    TString out_name = TString::Format("results/%s/muon_pT_split_%s", MC_name.Data(), source_flag.Data());
    out_name.ReplaceAll(".", "_");
    c1->SaveAs(out_name + ".png");

    delete c1;
}

void inspect_sources() {

    // TString MC_name = "DQ";
    // TString MC_name = "HF";
    TString MC_name = "genpurp";
    TString data_file = "results/" + MC_name + "/muonAOD.root";

    plot_source_matrix(data_file, MC_name);
    plot_source_pT(data_file, MC_name);
    plot_source_PID(data_file, MC_name, "fIsPrimary");
    plot_source_PID(data_file, MC_name, "fIsProducedInTransport");
}