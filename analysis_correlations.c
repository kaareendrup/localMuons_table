
void fillHist(const std::string& name, std::map<std::string, std::unique_ptr<TH1F>>& hists, double value, int n_bins, double x_min, double x_max) {

    // If histogram doesn't exist yet, create it
    if (hists.find(name) == hists.end()) {
        hists[name] = std::make_unique<TH1F>(name.c_str(), name.c_str(), n_bins, x_min, x_max);
        hists[name]->SetDirectory(nullptr); 
    }

    hists[name]->Fill(value);
}

void analysis_correlations() {

    TString data_name = "DQ";
    // TString data_name = "DQ_data";

    TString type = "gen";
    // TString type = "reco";

    bool is_MC = !(data_name == "DQ_data");

    // pT cuts
    // float pT_trigger_min = 1.0;
    float pT_trigger_min = 0.0;
    // float pT_trigger_max = 20.0;
    float pT_trigger_max = 12.0;
    // float pT_assoc_min = 3.0;
    // float pT_assoc_min = 2.0;
    // float pT_assoc_min = 1.0;
    float pT_assoc_min = 0.0;
    // float pT_assoc_max = 20.0;
    float pT_assoc_max = 12.0;

    // Eta cuts    
    float eta_trigger_min = -3.6;
    float eta_trigger_max = -2.5;
    float eta_assoc_min = -3.6;
    float eta_assoc_max = -2.5;

    // Histogram parameters
    int n_bins_mass = 100;
    int n_bins = 20;
    float signal_range_min = 2.7;
    float signal_range_max = 3.4;
    float deltaEta_min = -2.0;
    float deltaEta_max = 2.0;

    float signal_width = signal_range_max - signal_range_min;
    float background_range_min = signal_range_min - (signal_width / 2.0);
    float background_range_max = signal_range_max + (signal_width / 2.0);

    std::vector<float> segments_pt = {0.0, 0.5, 1.0, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0, 20.0};

    // Set outfile
    TFile* outFile = TFile::Open(TString::Format("results/%s/%s/analysis.root", data_name.Data(), type.Data()), "RECREATE");
    TTree* triggerCounts = new TTree("Correlations", "JPsi correlations");
    TTree* metaData = new TTree("MetaData", "Event selection metadata");

    // Define branches for output tree
    std::string category_out; 
    int count_out;
    std::vector<float> signal_ranges;
    std::vector<float> pTCuts = {pT_trigger_min, pT_trigger_max, pT_assoc_min, pT_assoc_max};
    std::vector<float> etaCuts = {eta_trigger_min, eta_trigger_max, eta_assoc_min, eta_assoc_max};
    triggerCounts->Branch("category", &category_out);
    triggerCounts->Branch("count", &count_out, "count/I");
    metaData->Branch("signal_range", &signal_ranges);
    metaData->Branch("segments_pt", &segments_pt);
    metaData->Branch("pTCuts", &pTCuts);
    metaData->Branch("etaCuts", &etaCuts);

    // Store signal range info
    for (float range : {background_range_min, signal_range_min, signal_range_max, background_range_max}) {
        signal_ranges.push_back(range);
    }
    metaData->Fill();

    // Load events
    TString data_file = TString::Format("results/%s/%s/eventmuons.root", data_name.Data(), type.Data());
    TFile *file = TFile::Open(data_file);
    if (!file || file->IsZombie()) {
        std::cerr << "Cannot open file\n";
        return;
    }

    TTree *tree = nullptr;
    file->GetObject("Triggers", tree);
    if (!tree) {
        std::cerr << "Tree 'Triggers' not found\n";
        return;
    }

    std::string *category = nullptr;
    double pT, eta, phi, mass;
    std::vector<double> *pT_assocs = nullptr;
    std::vector<double> *eta_assocs = nullptr;
    std::vector<double> *phi_assocs = nullptr;
    std::vector<int> *MotherPID = nullptr;

    tree->SetBranchAddress("category", &category);
    tree->SetBranchAddress("pT",  &pT);
    tree->SetBranchAddress("eta",  &eta);
    tree->SetBranchAddress("phi",  &phi);
    tree->SetBranchAddress("mass",  &mass);
    tree->SetBranchAddress("pT_assocs",  &pT_assocs);
    tree->SetBranchAddress("eta_assocs",  &eta_assocs);
    tree->SetBranchAddress("phi_assocs",  &phi_assocs);
    tree->SetBranchAddress("MotherPID", &MotherPID);

    std::map<std::string, int> trigger_counts;
    std::map<std::string, std::unique_ptr<TH1F>> invMassHists;
    std::map<std::string, std::unique_ptr<TH1F>> deltaEtaHists;
    std::map<std::string, std::unique_ptr<TH1F>> deltaPhiHists;
    std::map<std::string, std::unique_ptr<TH1F>> deltaPhiHistspT;

    // Loop over entries and create the necessary histograms
    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {

        tree->GetEntry(i);
        std::cout << "Processing entry " << i+1 << " of " << tree->GetEntries() << "\r" << std::flush;

        std::string category_str = *category;
        trigger_counts[category_str]++;
        fillHist(category_str + "_invMass", invMassHists, mass, n_bins_mass, 1.0, 5.0);

        // Calculate delta eta and phi for associated particles
        for (size_t j = 0; j < eta_assocs->size(); ++j) {

            double deltaEta = eta - eta_assocs->at(j);
            double deltaPhi = phi - phi_assocs->at(j);
            if (deltaPhi > 3.0 / 2.0 * M_PI) {
                deltaPhi -= 2.0 * M_PI;
            }
            if (deltaPhi < -0.5 * M_PI) {
                deltaPhi += 2.0 * M_PI;
            }

            std::string full_category = category_str;

            // pT cuts
            if (pT < pT_trigger_min || pT > pT_trigger_max) continue;
            if (pT_assocs->at(j) < pT_assoc_min || pT_assocs->at(j) > pT_assoc_max) continue;

            // Eta cuts
            if (eta < eta_trigger_min || eta > eta_trigger_max) continue;
            if (eta_assocs->at(j) < eta_assoc_min || eta_assocs->at(j) > eta_assoc_max) continue;

            if (is_MC) {
                if ((std::abs(MotherPID->at(j)) >= 411 && std::abs(MotherPID->at(j)) <= 445) || 
                    (std::abs(MotherPID->at(j)) >= 4101 && std::abs(MotherPID->at(j)) <= 4444)) { // J/psi from charm
                    full_category += "_charm";
                } else { // J/psi not from charm
                    full_category += "_noncharm";
                }
            } else {
                if (mass > signal_range_min && mass < signal_range_max) {
                    full_category += "_signal";
                } else if (mass > background_range_min && mass < signal_range_min || 
                           mass > signal_range_max && mass < background_range_max) {
                    full_category += "_background";
                } else {
                    continue; // Skip if not in signal or background range
                }
            }

            fillHist(full_category + "_deltaEta", deltaEtaHists, deltaEta, n_bins, deltaEta_min, deltaEta_max);
            fillHist(full_category + "_deltaPhi", deltaPhiHists, deltaPhi, n_bins, -0.5 * M_PI, 3.0 / 2.0 * M_PI);

            for (int p = 0; p < 10; ++p) {
                float ptmin = segments_pt[p];
                float ptmax = segments_pt[p + 1];
                if (pT_assocs->at(j) >= ptmin && pT_assocs->at(j) < ptmax) {
                    fillHist(std::format("{}_deltaPhi_pT_{:.1f}_{:.1f}", full_category, ptmin, ptmax), deltaPhiHistspT, deltaPhi, n_bins, -0.5 * M_PI, 3.0 / 2.0 * M_PI);
                }
            }
        }
    }
    std::cout << std::endl;

    outFile->cd();
    for (auto& mHist : invMassHists) {
        mHist.second->Write();
    }
    for (auto& dEtaHist : deltaEtaHists) {
        dEtaHist.second->Write();
    }
    for (auto& dPhiHist : deltaPhiHists) {
        dPhiHist.second->Write();
    }
    for (auto& dPhiHist : deltaPhiHistspT) {
        dPhiHist.second->Write();
    }

    for (auto& count : trigger_counts) {
        category_out = count.first;
        count_out = count.second;
        triggerCounts->Fill();
    }
    triggerCounts->Write();
    metaData->Write();
    outFile->Close();
}