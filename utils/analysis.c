
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

std::map<std::string, int> get_trigger_counts(TTree* tree) {
    // Load trigger counts from the "Triggers" tree into a map for easy access when normalizing histograms
    std::map<std::string, int> trigger_counts;
    std::string *category = nullptr;
    int count;

    tree->SetBranchAddress("category", &category);
    tree->SetBranchAddress("count", &count);

    for (Long64_t i = 0; i < tree->GetEntries(); ++i) {

        tree->GetEntry(i);
        std::string category_str = *category;
        trigger_counts[category_str] = count;
    }

    return trigger_counts;
}
