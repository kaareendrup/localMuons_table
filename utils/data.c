
bool charm_beauty_cut(Long64_t motherPDG) {
    return ((std::abs(motherPDG) == 443) || (std::abs(motherPDG) == 100443) ||
            (std::abs(motherPDG) >= 411 && std::abs(motherPDG) <= 445) || 
            (std::abs(motherPDG) >= 4101 && std::abs(motherPDG) <= 4444) || 
            (std::abs(motherPDG) >= 511 && std::abs(motherPDG) <= 557) || 
            (std::abs(motherPDG) >= 5101 && std::abs(motherPDG) <= 5554));
}
