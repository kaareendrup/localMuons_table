
#include <vector>

double getRapidity(double pT, double eta) {
    double mJPsi = 3.096916; // J/Psi mass in GeV/c^2
    return log((sqrt(pow(mJPsi, 2) + (pow(pT, 2) * pow(cosh(eta), 2))) + pT * sinh(eta)) / (sqrt(pow(mJPsi, 2) + pow(pT, 2))));
}

double getDeltaY(double pT, double eta_min, double eta_max) {

    double y_min = getRapidity(pT, eta_min);
    double y_max = getRapidity(pT, eta_max);
    // std::cout << "pT = " << pT << " GeV/c, y_min = " << y_min << ", y_max = " << y_max << std::endl;
    return y_max - y_min;
}

void fillHist(const TString& name, std::map<TString, std::unique_ptr<TH1F>>& hists, double value, int n_bins, double x_min, double x_max, double weight = 1.0) {

    // If histogram doesn't exist yet, create it
    if (hists.find(name) == hists.end()) {
        hists[name] = std::make_unique<TH1F>(name, name, n_bins, x_min, x_max);
        hists[name]->SetDirectory(nullptr); 
    }

    hists[name]->Fill(value, weight);
}

double get_weight(double pT, const std::vector<double>& pT_bins, const std::vector<double>& efficiency, double eta_trigger_min, double eta_trigger_max) {
    
    int eff_bin = -1;
    for (size_t b = 0; b < pT_bins.size() - 1; ++b) {
        if (pT >= pT_bins[b] && pT < pT_bins[b + 1]) {
            eff_bin = b;
            break;
        }
    }
    if (eff_bin == -1) {
        std::cerr << "Warning: pT " << pT << " is out of efficiency histogram range!\n";
        return 1.0; // No correction if out of range
    }
    return 1.0 / (efficiency[eff_bin] * getDeltaY(pT, eta_trigger_min, eta_trigger_max));
}