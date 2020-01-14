#include <iostream>

#include "Math/VectorUtil.h"
#include "MetSelections.h"
#include "Config.h"

using namespace tas;


bool passesMETfilters(bool is_data) {

  if (gconf.year == 2016) {
    if (!Flag_goodVertices()) return false; // primary vertex filter?
    if (!Flag_globalSuperTightHalo2016Filter()) return false; // beam halo filter
    if (!Flag_HBHENoiseFilter()) return false;
    if (!Flag_HBHENoiseIsoFilter()) return false;
    if (!Flag_EcalDeadCellTriggerPrimitiveFilter()) return false; // ECAL TP filter
    if (!Flag_BadPFMuonFilter()) return false; // bad PF muon filter
    if (!Flag_BadChargedCandidateFilter()) return false;  // bad charged hadron filter
    if (is_data) {
      if (!Flag_eeBadScFilter()) return false;
    }
  } else if (gconf.year == 2017) {
    if (!Flag_goodVertices()) return false; // primary vertex filter
    if (!Flag_globalSuperTightHalo2016Filter()) return false; // beam halo filter
    if (!Flag_HBHENoiseFilter()) return false; // HBHE noise filter
    if (!Flag_HBHENoiseIsoFilter()) return false; // HBHEiso noise filter
    if (!Flag_EcalDeadCellTriggerPrimitiveFilter()) return false; // ECAL TP filter
    if (!Flag_BadPFMuonFilter()) return false; // bad PF muon filter
    if (!Flag_BadChargedCandidateFilter()) return false;  // bad charged hadron filter
    if (!Flag_ecalBadCalibFilterV2()) return false;  // ECAL bad calibration filter update
    if (is_data) {
      if (!Flag_eeBadScFilter()) return false;
    }
  } else if (gconf.year == 2018) {
    if (!Flag_goodVertices()) return false; // primary vertex filter
    if (!Flag_globalSuperTightHalo2016Filter()) return false; // beam halo filter
    if (!Flag_HBHENoiseFilter()) return false; // HBHE noise filter
    if (!Flag_HBHENoiseIsoFilter()) return false; // HBHEiso noise filter
    if (!Flag_EcalDeadCellTriggerPrimitiveFilter()) return false; // ECAL TP filter
    if (!Flag_BadPFMuonFilter()) return false; // bad PF muon filter
    if (!Flag_BadChargedCandidateFilter()) return false;  // bad charged hadron filter
    if (!Flag_ecalBadCalibFilterV2()) return false;  // ECAL bad calibration filter update
    if (is_data) {
      if (!Flag_eeBadScFilter()) return false;
    }
  } else {
    std::cout << "[MetSelections] >> uh, I don't know what year this is" << std::endl;
  }

  // otherwise, good
  return true;
}

