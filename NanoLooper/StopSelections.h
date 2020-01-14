#ifndef StopSelections_H
#define StopSelections_H

#include "TRandom3.h"
#include "../NanoCORE/Nano.h"
#include "../NanoCORE/Config.h"
// #include "PhysicsObjects.h"

const float k_minPt_el_veto = 5;
const float k_minPt_mu_veto = 5;
const float k_maxEta_el_veto = 2.4;
const float k_maxEta_mu_veto = 2.4;
const float k_minPt_el_good = 20;
const float k_minPt_mu_good = 20;
const float k_maxEta_el_good = 1.4442;
const float k_maxEta_mu_good = 2.4;
const float k_minPt_jet = 30;
const float k_maxEta_jet = 2.4;

const float mZ = 91.1876;

extern TRandom3* randomGenerator;

bool passTriggerSelections(int trigtype);

vector<int> getJetIdxs(const vector<int>& mus, const vector<int>& els, const vector<int>& phs, bool reapplyJEC = false, bool isSim = true);
std::tuple<vector<int>, vector<int>> getMuonIdxs(bool applyRocCorr = true, float* shiftx = nullptr, float* shifty = nullptr);
std::tuple<vector<int>, vector<int>> getElectronIdxs();

float getJetCorrectionFactorFromFile(int jetidx, int injet, bool applyJER = true);

#endif

