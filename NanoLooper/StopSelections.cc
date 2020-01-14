#include "StopSelections.h"
// #include "Utilities.h"
#include <bitset>

namespace {

inline bool isCloseObject(const float eta1, const float phi1, const float eta2, const float phi2, const float conesize, float* deltaR = nullptr) {
  const float PI = TMath::Pi();
  float deltaEta = fabs(eta1 - eta2);
  if (deltaEta > conesize) return false;
  float deltaPhi = fabs(phi1 - phi2);
  if (deltaPhi > PI) deltaPhi = 2*PI - deltaPhi;
  if (deltaPhi > conesize) return false;
  float deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
  if (deltaR2 > conesize*conesize) return false;
  if (deltaR) *deltaR = sqrt(deltaR2);

  return true;
}

}

using namespace std;
using namespace tas;

// RoccoR* muoncorr = nullptr;
TRandom3* randomGenerator = nullptr;

bool passTriggerSelections(int trigtype) {
  if (!gconf.is_data) return true;  // not using the trigger emulatino in MC

  // if (trigtype == 1) {
  //   // Muon triggers
  //   switch (gconf.year) {
  //     case 2016:
  //       return ( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL() || HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ() ||
  //                HLT_IsoMu24() || HLT_IsoTkMu24());
  //     case 2017:
  //       return ( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ() ||
  //                ((run() >= 299337)? HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8() : false) ||
  //                HLT_IsoMu24() || HLT_IsoMu27());
  //     case 2018:
  //     default:
  //       return ( HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8() ||
  //                HLT_IsoMu24());
  //   }
  // } else if (trigtype == 2) {
  //   // Electron triggers
  //   switch (gconf.year) {
  //     case 2016:
  //       return ( HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() ||
  //                HLT_Ele27_WPTight_Gsf() );
  //     case 2017:
  //     case 2018:
  //     default:
  //       return ( HLT_Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() ||
  //                ((run() >= 302026)? HLT_Ele32_WPTight_Gsf() : HLT_Ele35_WPTight_Gsf()) ||
  //                HLT_Ele32_WPTight_Gsf_L1DoubleEG() );
  //   }
  // } else if (trigtype == 3) {
  //   // emu triggers
  //   switch (gconf.year) {
  //     case 2016:
  //     case 2017:
  //     case 2018:
  //     default:
  //       return (
  //           HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL() || HLT_Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL_DZ() ||
  //           HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL() || HLT_Mu12_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ() ||
  //           HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL()  || HLT_Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL_DZ()  );
  //   }
  // } else if (trigtype == 4) {
  //   // Photon triggers
  //   bool prescaled = false;
  //   switch (gconf.year) {
  //     case 2016:
  //       // prescaled photon triggers
  //       if (false) {
  //         prescaled = (
  //             HLT_Photon50_R9Id90_HE10_IsoM() ||  // eff lumi: 0.
  //             HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50() ||
  //             HLT_Photon75_R9Id90_HE10_IsoM() ||
  //             HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3() ||
  //             HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3() ||
  //             HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3() ||
  //             HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3() ||
  //             HLT_Photon90_R9Id90_HE10_IsoM() ||
  //             HLT_Photon120_R9Id90_HE10_IsoM() );
  //       }

  //       // unprescaled in 2016, prescaled in 2017 & 2018
  //       return ( prescaled || HLT_Photon165_R9Id90_HE10_IsoM() || HLT_Photon250_NoHE() || HLT_Photon300_NoHE());
  //       // HLT_Photon175() || HLT_Photon200()  // unprescaled, but not use
  //       // HLT_Photon250_NoHE(); // in 2016 only that doesn't exist
  //     case 2017:
  //     case 2018:
  //     default:
  //       if (false) {
  //         prescaled = (
  //             HLT_Photon50_R9Id90_HE10_IsoM() ||  // eff lumi: 0.
  //             HLT_Photon50_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3_PFMET50() ||
  //             HLT_Photon75_R9Id90_HE10_IsoM() ||
  //             HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ300_PFJetsMJJ400DEta3() ||
  //             HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_CaloMJJ400_PFJetsMJJ600DEta3() ||
  //             HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ300DEta3() ||
  //             HLT_Photon75_R9Id90_HE10_IsoM_EBOnly_PFJetsMJJ600DEta3() ||
  //             HLT_Photon90_R9Id90_HE10_IsoM() ||
  //             HLT_Photon120_R9Id90_HE10_IsoM() );
  //       }

  //       // unprescaled in 2016, prescaled in 2017 & 2018
  //       return ( prescaled || HLT_Photon300_NoHE() );
  //       // HLT_Photon175() || HLT_Photon200()  // unprescaled, but not use
  //   }
  //   // Lepton triggers that doesn't exist
  //   // HLT_IsoTkMu24() ||
  //   // HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL() ||
  //   // HLT_Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL_DZ() ||
  // }

  // return false;
  return true;
}

std::tuple<vector<int>, vector<int>> getElectronIdxs() {

  vector<int> vetoElectrons;
  vector<int> goodElectrons;

  for (unsigned i = 0; i < Electron_pt().size(); ++i) {
    // double const absEtaSc = fabs(Electron_deltaEtaSC()[i] + Electron_eta()[i]);
    double const absEta = fabs(Electron_eta()[i]);
    bool const passVetoId = (Electron_cutBased()[i] >= 1);

    if (Electron_pt()[i] < k_minPt_el_veto || absEta > k_maxEta_el_veto || !passVetoId)
      continue;

    vetoElectrons.push_back(i);

    // if (absEtaSc > 1.4442 and absEtaSc < 1.5660) continue; // EB-EE gap

    bool const passGoodId = (Electron_cutBased()[i] >= 3);  // medium ID
    if (Electron_pt()[i] < k_minPt_el_good or absEta > k_maxEta_el_good or not passGoodId)
      continue;

    goodElectrons.emplace_back(i);
  }

  std::sort(vetoElectrons.begin(), vetoElectrons.end(), [](int a, int b){ return Electron_pt()[a] > Electron_pt()[b]; });
  std::sort(goodElectrons.begin(), goodElectrons.end(), [](int a, int b){ return Electron_pt()[a] > Electron_pt()[b]; });

  return {goodElectrons, vetoElectrons};
}

std::tuple<vector<int>, vector<int>> getMuonIdxs(bool applyRocCorr, float* shiftx, float* shifty) {

  vector<int> vetoMuons;
  vector<int> goodMuons;

  for (unsigned i = 0; i < Muon_pt().size(); ++i) {
    // Veto ID as per https://twiki.cern.ch/twiki/bin/view/CMS/SWGuideMuonIdRun2#Veto_Muon
    // bool const passVetoId = Muon_isPFcand()[i] && Muon_isGlobal()[i] && Muon_isTracker()[i];
    bool const passVetoId = (Muon_isPFcand()[i] && Muon_isGlobal()[i] && Muon_isTracker()[i]);
    bool const passVetoIso = Muon_miniIsoId()[i] >= 1;  // FIXME: 
    if (fabs(Muon_eta()[i]) > k_maxEta_mu_veto or not passVetoId or passVetoIso)
      continue;

    if (Muon_pt()[i] < k_minPt_mu_veto)  // minPtVeto
      continue;

    vetoMuons.emplace_back(i);

    bool const passGoodId = (Muon_mediumPromptId()[i] && (Muon_miniIsoId()[i] >= 3));

    if (Muon_pt()[i] < k_minPt_mu_good or Muon_eta()[i] > k_maxEta_mu_good or not passGoodId or passVetoIso) // minPtGood
      continue;

    goodMuons.emplace_back(i);
  }

  std::sort(vetoMuons.begin(), vetoMuons.end(), [](int a, int b){ return Muon_pt()[a] > Muon_pt()[b]; });
  std::sort(goodMuons.begin(), goodMuons.end(), [](int a, int b){ return Muon_pt()[a] > Muon_pt()[b]; });

  return {goodMuons, vetoMuons};
}


float getJetCorrectionFactorFromFile(int jetidx, LorentzVector injet, bool applyJER) {
    double corrFactor = 1.;
    return corrFactor;
}


// return jets, bjets
tuple<vector<int>, vector<int>> getJetIdxs(const vector<int>& mus, const vector<int>& els, bool reapplyJEC, bool isSim) {
  vector<int> jets;
  vector<int> bjets;

  for (unsigned j = 0; j < Jet_pt().size(); ++j) {
    if (not (Jet_jetId()[j] & 0b0010))  // bit1 always false in 2017 since it does not exist
      continue;

    // Perform angular cleaning w.r.t. recognized leptons and photons
    for (int l : mus) {
      if (isCloseObject(Jet_eta()[j], Jet_phi()[j], Muon_eta()[l], Muon_phi()[l], 0.4))
        goto end_of_loop_jets;
    }
    for (int l : els) {
      if (isCloseObject(Jet_eta()[j], Jet_phi()[j], Electron_eta()[l], Electron_phi()[l], 0.4))
        goto end_of_loop_jets;
    }

    // if (reapplyJEC) jet *= getJetCorrectionFactorFromFile();

    // Kinematical cuts for jets to be stored in the collection
    if (Jet_pt()[j] < k_minPt_jet or fabs(Jet_eta()[j]) > k_maxEta_jet)
      continue;

    jets.emplace_back(j);

    if (Jet_btagDeepB()[j] > gconf.WP_DeepCSV_medium)
      bjets.emplace_back(j);

 end_of_loop_jets:;
  }

  // Make sure jets are sorted in pt
  std::sort(jets.begin(), jets.end(), [](int a, int b){ return Jet_pt()[a] > Jet_pt()[b]; });
  std::sort(bjets.begin(), bjets.end(), [](int a, int b){ return Jet_pt()[a] > Jet_pt()[b]; });

  return {jets, bjets};
}
