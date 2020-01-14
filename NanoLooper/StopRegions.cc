#include "StopRegions.h"

const float fInf = std::numeric_limits<float>::max();
const float wpResTop = 0.95;
const float wpDeepTop = 0.40;


std::vector<SR> getStopControlRegionsNoBTags(std::vector<SR>&& SRvec) {

  std::vector<SR> CRvec;
  for (SR cr : SRvec) {
    cr.SetName(cr.GetName().replace(0, 2, "cr0b"));
    cr.SetAllowDummyVars(1);
    if (cr.VarExists("nbjet")) cr.RemoveVar("nbjet");
    if (cr.VarExists("ntbtag")) cr.RemoveVar("ntbtag");
    if (cr.GetName().find("base") != std::string::npos || cr.GetName().find("cr0bI") == 0) {
      cr.SetVar("nbjet", 0, 1);
    } else {
      cr.SetVar("nbtag", 0, 1);
    }
    if (cr.GetLowerBound("tmod") < 10)
      cr.SetVar("nbjet", 0, 1);
    cr.SetVar("mlb_0b", cr.GetLowerBound("mlb"), cr.GetUpperBound("mlb"));
    cr.RemoveVar("mlb");
    CRvec.emplace_back(cr);
  }

  return CRvec;
}

std::vector<SR> getStopControlRegionsDilepton(std::vector<SR>&& SRvec) {

  std::vector<SR> CRvec;

  for (SR cr : SRvec) {
    cr.SetName(cr.GetName().replace(0, 2, "cr2l"));
    cr.SetAllowDummyVars(1);
    for (std::string var : {"met", "mt", "dphijmet", "nlep", "tmod", "dphilmet"}) {
      if (!cr.VarExists(var)) continue;
      cr.SetVar(var+"_rl", cr.GetLowerBound(var), cr.GetUpperBound(var));
      cr.RemoveVar(var);
    }
    cr.RemoveVar("passvetos");
    cr.SetVar("nvlep", 2, fInf);
    cr.SetVar("nlep_rl", 2, fInf);
    cr.SetVar("mt2_ll", 0, 100);  // to avoid overlap with the stop-2l SR
    CRvec.emplace_back(cr);
  }

  return CRvec;
}

// Old 2016 Analysis search bins
std::vector<SR> getStopSignalRegionsTopological() {

  SR srbase;
  srbase.SetAllowDummyVars(1);
  srbase.SetName("srbase");
  srbase.SetVar("mt", 150, fInf);
  srbase.SetVar("met", 250, fInf);
  srbase.SetVar("nlep", 1, 2);
  srbase.SetVar("nvlep", 1, 2);
  srbase.SetVar("passvetos", 1, 2);
  srbase.SetVar("njet", 2, fInf);
  srbase.SetVar("nbjet", 1, fInf);
  srbase.SetVar("ntbtag", 0, fInf);
  // srbase.SetVar("nbtag", 1, fInf);
  srbase.SetVar("mlb", 0, fInf);
  srbase.SetVar("tmod", -fInf, fInf);
  srbase.SetVar("dphijmet", 0.8, 3.14159);
  srbase.SetMETBins({0, 250, 350, 450, 550, 650, 800, 1500});

  SR sr;
  std::vector<SR> SRvec;

  SRvec.emplace_back(srbase);

  // SR 2-3j

  sr = srbase;
  sr.SetName("srA");
  sr.SetDetailName("2to3j_tmod10toInf_mlb0to175");
  sr.SetVar("njet", 2, 4);
  sr.SetVar("tmod", 10, fInf);
  sr.SetVar("mlb", 0, 175);
  sr.SetVar("ntbtag", 0, fInf);
  sr.SetMETBins({250, 350, 450, 600, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srB");
  sr.SetDetailName("2to3j_tmod10toInf_mlb175toInf");
  sr.SetVar("njet", 2, 4);
  sr.SetVar("tmod", 10, fInf);
  sr.SetVar("mlb", 175, fInf);
  sr.SetVar("ntbtag", 1, fInf);
  sr.SetMETBins({250, 450, 600, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srC");
  sr.SetDetailName("geq4j_tmodlt0_mlb0to175");
  sr.SetVar("njet", 4, fInf);
  sr.SetVar("tmod", -fInf, 0);
  sr.SetVar("mlb", 0, 175);
  sr.SetVar("ntbtag", 0, fInf);
  sr.SetMETBins({250, 350, 450, 550, 650, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srD");
  sr.SetDetailName("geq4j_tmodlt0_mlb175toInf");
  sr.SetVar("njet", 4, fInf);
  sr.SetVar("tmod", -fInf, 0);
  sr.SetVar("mlb", 175, fInf);
  sr.SetVar("ntbtag", 1, fInf);
  sr.SetMETBins({250, 350, 450, 550, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srE");
  sr.SetDetailName("geq4j_tmod0to10_mlb0to175");
  sr.SetVar("njet", 4, fInf);
  sr.SetVar("tmod", 0, 10);
  sr.SetVar("mlb", 0, 175);
  sr.SetVar("ntbtag", 0, fInf);
  sr.SetMETBins({250, 350, 550, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srF");
  sr.SetDetailName("geq4j_tmod0to10_mlb175toInf");
  sr.SetVar("njet", 4, fInf);
  sr.SetVar("tmod", 0, 10);
  sr.SetVar("mlb", 175, fInf);
  sr.SetVar("ntbtag", 1, fInf);
  sr.SetMETBins({250, 450, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srG");
  sr.SetDetailName("geq4j_tmod10toInf_mlb0to175");
  sr.SetVar("njet", 4, fInf);
  sr.SetVar("tmod", 10, fInf);
  sr.SetVar("mlb", 0, 175);
  sr.SetVar("ntbtag", 0, fInf);
  sr.SetMETBins({250, 350, 450, 600, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srH");
  sr.SetDetailName("geq4j_tmod10toInf_mlb175toInf");
  sr.SetVar("njet", 4, fInf);
  sr.SetVar("tmod", 10, fInf);
  sr.SetVar("mlb", 175, fInf);
  sr.SetVar("ntbtag", 1, fInf);
  sr.SetMETBins({250, 450, 1500});
  SRvec.emplace_back(sr);

  // Compressed regions
  sr = srbase;
  sr.SetName("srI");
  sr.SetDetailName("geq5j_lpt0to150_j1notb");
  sr.SetVar("mt", 150, fInf);
  sr.SetVar("njet", 5, fInf);
  sr.SetVar("nbjet", 1, fInf);
  sr.SetVar("lep1pt", 0, 150);
  sr.SetVar("dphilmet", 0, 2.0);
  sr.SetVar("dphijmet", 0.5, 3.1416);
  sr.SetVar("j1passbtag", 0, 1);  // Require j1 not b-tagged
  sr.SetMETBins({250, 350, 450, 550, 1500});
  SRvec.emplace_back(sr);

  return SRvec;
}


std::vector<SR> getStopCrosscheckRegionsEMuRun2() {
  std::vector<SR> CRvec;

  SR crbase;

  crbase.SetName("cremu_base");
  crbase.SetVar(nvlep, 1, 3);
  crbase.SetVar(met, 50, fInf);
  crbase.SetVar(njet, 2, fInf);
  crbase.SetMETBins({50, 100, 150, 200, 250, 350, 450, 600, 800, 1500});
  crbase.SetAllowDummyVars(1);

  // The strict lepton selections are in the looper

  SR cr = crbase;

  cr.SetName("cremuA0");
  cr.SetDetailName("emu_ge2j_ge0b_met50toInf");
  cr.SetVar(nbjet, 0, fInf);
  CRvec.emplace_back(cr);

  cr.SetName("cremuA1");
  cr.SetDetailName("emu_ge2j_ge1b_met50toInf");
  cr.SetVar(nbjet, 1, fInf);
  CRvec.emplace_back(cr);

  cr.SetName("cremuA2");
  cr.SetDetailName("emu_ge2j_ge2b_met50toInf");
  cr.SetVar(nbjet, 2, fInf);
  CRvec.emplace_back(cr);

  cr.SetName("cremuA3");
  cr.SetDetailName("emu_ge2j_2b_met50toInf");
  cr.SetVar(nbjet, 2, 3);
  CRvec.emplace_back(cr);

  cr.SetName("cremuB1");
  cr.SetDetailName("emu_ge1j_1b_met50toInf");
  cr.SetVar(njet, 1, fInf);
  cr.SetVar(nbjet, 1, 2);
  CRvec.emplace_back(cr);

  cr.SetName("cremuB2");
  cr.SetDetailName("emu_1j_1b_met50toInf");
  cr.SetVar(njet, 1, 2);
  cr.SetVar(nbjet, 1, 2);
  CRvec.emplace_back(cr);

  cr.SetName("cremuC2");
  cr.SetDetailName("emu_2to3j_ge1b_met50toInf");
  cr.SetVar(njet, 2, 4);
  cr.SetVar(nbjet, 1, fInf);
  CRvec.emplace_back(cr);

  cr.SetName("cremuC3");
  cr.SetDetailName("emu_ge3j_ge1b_met50toInf");
  cr.SetVar(njet, 3, fInf);
  cr.SetVar(nbjet, 1, fInf);
  CRvec.emplace_back(cr);

  cr.SetName("cremuC4");
  cr.SetDetailName("emu_ge4j_ge1b_met50toInf");
  cr.SetVar(njet, 4, fInf);
  cr.SetVar(nbjet, 1, fInf);
  CRvec.emplace_back(cr);

  cr.SetName("cremuC5");
  cr.SetDetailName("emu_ge5j_ge1b_met50toInf");
  cr.SetVar(njet, 5, fInf);
  cr.SetVar(nbjet, 1, fInf);
  CRvec.emplace_back(cr);

  // cr = crbase;
  // cr.SetName("cremuM1");
  // cr.SetDetailName("ge2j_ge1b_met50to100");
  // cr.SetVar(met, 50, 100);
  // cr.SetVar(njet, 2, fInf);
  // cr.SetVar(nbjet, 1, fInf);
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuM2");
  // cr.SetDetailName("ge2j_ge1b_met150toInf");
  // cr.SetVar(met, 150, 250);
  // cr.SetVar(njet, 2, fInf);
  // cr.SetVar(nbjet, 1, fInf);
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuM3");
  // cr.SetDetailName("ge2j_ge1b_met150toInf");
  // cr.SetVar(met, 250, fInf);
  // cr.SetVar(njet, 2, fInf);
  // cr.SetVar(nbjet, 1, fInf);
  // CRvec.emplace_back(cr);


  // // SRs
  // cr.SetName("cremuA");
  // cr.SetDetailName("2to3j_ge1b_met50toInf_tmod10toInf_mlb0to175");
  // cr.SetVar(njet, 2, 4);
  // cr.SetVar(nbjet, 1, fInf);
  // cr.SetVar(mlb, 0, 175);
  // cr.SetVar(tmod, 10, fInf);
  // cr.SetMETBins({50, 100, 150, 200, 250, 350, 450, 600, 750, 1500});
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuB");
  // cr.SetDetailName("2to3j_ge1b_met50toInf_tmod10toInf_mlb175toInf");
  // cr.SetVar(njet, 2, 4);
  // cr.SetVar(nbjet, 1, fInf);
  // cr.SetVar(mlb, 175, fInf);
  // cr.SetVar(tmod, 10, fInf);
  // cr.SetMETBins({50, 100, 150, 200, 250, 450, 700, 1500});
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuC");
  // cr.SetDetailName("ge4j_ge1b_met50toInf_tmodInfto0_mlb0to175");
  // cr.SetVar(njet, 4, fInf);
  // cr.SetVar(nbjet, 1, fInf);
  // cr.SetVar(mlb, 0, 175);
  // cr.SetVar(tmod, -fInf, 0);
  // cr.SetMETBins({50, 100, 150, 200, 250, 350, 450, 550, 650, 800, 1500});
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuE");
  // cr.SetDetailName("ge4j_ge1b_met50toInf_tmod0to10_mlb0to175");
  // cr.SetVar(njet, 4, fInf);
  // cr.SetVar(nbjet, 1, fInf);
  // cr.SetVar(mlb, 0, 175);
  // cr.SetVar(tmod, 0, 10);
  // cr.SetMETBins({50, 100, 150, 200, 250, 350, 450, 600, 1500});
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuG");
  // cr.SetDetailName("ge4j_ge1b_met50toInf_tmod10toInf_mlb0to175");
  // cr.SetVar(njet, 4, fInf);
  // cr.SetVar(nbjet, 1, fInf);
  // cr.SetVar(mlb, 0, 175);
  // cr.SetVar(tmod, 0, 10);
  // cr.SetMETBins({50, 100, 150, 200, 250, 350, 450, 550, 750, 1500});
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuH");
  // cr.SetDetailName("ge4j_ge1b_met50toInf_tmod10toInf_mlb175toInf");
  // cr.SetVar(njet, 4, fInf);
  // cr.SetVar(nbjet, 1, fInf);
  // cr.SetVar(mlb, 0, 175);
  // cr.SetVar(tmod, 0, 10);
  // cr.SetMETBins({50, 100, 150, 200, 250, 500, 800, 1500});
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuI");
  // cr.SetDetailName("ge5j_ge1b_met50toInf");
  // cr.SetVar(njet, 5, fInf);
  // cr.SetVar(nbjet, 1, fInf);
  // cr.RemoveVar(mlb);
  // cr.RemoveVar(tmod);
  // cr.SetVar(passlmetcor, 1, 2);
  // cr.SetVar(j1passbtag, 0, 1);  // Require j1 not b-tagged
  // cr.SetMETBins({50, 100, 150, 200, 250, 350, 450, 550, 750, 1500});
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuJ");
  // cr.SetDetailName("ge5j_ge1b_met50toInf");
  // cr.SetVar(njet, 5, fInf);
  // cr.SetVar(nbjet, 0, fInf);
  // cr.SetVar(nsbtag, 1, fInf);
  // cr.SetMETBins({50, 100, 150, 200, 250, 350, 450, 550, 800, 1500});
  // CRvec.emplace_back(cr);

  // cr.SetVar(mt, 150, fInf);
  // cr.SetName("cremuC0");
  // cr.SetDetailName("ge2j_0b_met50toInf_mt150toInf");
  // cr.SetVar(nbjet, 0, 1);
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuC1");
  // cr.SetDetailName("ge2j_1b_met50toInf_mt150toInf");
  // cr.SetVar(nbjet, 1, 2);
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuC2");
  // cr.SetDetailName("ge2j_ge2b_met50toInf_mt150toInf");
  // cr.SetVar(nbjet, 2, fInf);
  // CRvec.emplace_back(cr);

  // cr.SetVar(mt, 0, fInf);
  // cr.SetVar(tmod, 0, 10);

  // cr.SetName("cremuD0");
  // cr.SetDetailName("ge2j_0b_met50toInf_tmod0to10");
  // cr.SetVar(nbjet, 0, 1);
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuD1");
  // cr.SetDetailName("ge2j_ge1b_met50toInf_tmod0to10");
  // cr.SetVar(nbjet, 1, fInf);
  // CRvec.emplace_back(cr);

  // cr.SetVar(tmod, 10, fInf);

  // cr.SetName("cremuE0");
  // cr.SetDetailName("ge2j_0b_met50toInf_tmod10toInf");
  // cr.SetVar(nbjet, 0, 1);
  // CRvec.emplace_back(cr);

  // cr.SetName("cremuE1");
  // cr.SetDetailName("ge2j_ge1b_met50toInf_tmod10toInf");
  // cr.SetVar(nbjet, 1, fInf);
  // CRvec.emplace_back(cr);


  return CRvec;
}


std::vector<SR> getStopSignalRegionsRun2() {

  SR srbase;
  srbase.SetAllowDummyVars(1);
  srbase.SetName("srbase");
  srbase.SetVar(mt, 150, fInf);
  srbase.SetVar(met, 250, fInf);
  srbase.SetVar(nlep, 1, 2);
  srbase.SetVar(nvlep, 1, 2);
  srbase.SetVar(passvetos, 1, 2);
  srbase.SetVar(njet, 2, fInf);
  srbase.SetVar(nbjet, 1, fInf);
  srbase.SetVar(mlb, 0, fInf);
  srbase.SetVar(tmod, -fInf, fInf);
  srbase.SetVar(dphijmet, 0.8, 3.1416);
  srbase.SetMETBins({0, 250, 350, 450, 550, 650, 800, 1500});

  SR sr;
  std::vector<SR> SRvec;

  SRvec.emplace_back(srbase);

  // Start of nominal signal regions
  sr = srbase;
  sr.ReplaceVar(nbjet, nbtag);

  // SR 2-3j

  sr.SetName("srA0");
  sr.SetDetailName("2to3j_tmod10toInf_mlb0to175_nottag");
  sr.SetVar(njet, 2, 4);
  sr.SetVar(tmod, 10, fInf);
  sr.SetVar(mlb, 0, 175);
  sr.SetVar(met, 600, fInf);
  sr.SetMETBins({600, 750, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srA1");
  sr.SetDetailName("2to3j_tmod10toInf_mlb0to175_nottag");
  sr.SetVar(njet, 2, 4);
  sr.SetVar(tmod, 10, fInf);
  sr.SetVar(mlb, 0, 175);
  sr.SetVar(deepttag, -1, wpDeepTop);
  sr.SetVar(met, 350, 600);
  sr.SetMETBins({350, 450, 600});
  SRvec.emplace_back(sr);

  sr.SetName("srA2");
  sr.SetDetailName("2to3j_tmod10toInf_mlb0to175_deepttag");
  sr.SetVar(deepttag, wpDeepTop, 1);
  sr.SetVar(met, 250, 600);
  sr.SetMETBins({250, 600});
  SRvec.emplace_back(sr);


  sr.SetName("srB");
  sr.SetDetailName("2to3j_tmod10toInf_mlb175toInf_nottag");
  sr.SetVar(njet, 2, 4);
  sr.SetVar(tmod, 10, fInf);
  sr.SetVar(mlb, 175, fInf);
  sr.SetVar(deepttag, -1, 1);
  sr.SetVar(met, 250, fInf);
  sr.SetMETBins({250, 450, 700, 1500});
  SRvec.emplace_back(sr);

  // SR ge4j

  sr.SetName("srC");
  sr.SetDetailName("ge4j_tmodlt0_mlb0to175_met550");
  sr.SetVar(njet, 4, fInf);
  sr.SetVar(tmod, -fInf, 0);
  sr.SetVar(mlb, 0, 175);
  sr.SetVar(deepttag, -2, 1);
  sr.SetVar(tfttag, -2, 1);
  sr.SetVar(met, 350, fInf);
  sr.SetMETBins({350, 450, 550, 650, 800, 1500});
  SRvec.emplace_back(sr);


  sr.SetName("srD");
  sr.SetDetailName("ge4j_tmodlt0_mlb175toInf");
  sr.SetVar(njet, 4, fInf);
  sr.SetVar(tmod, -fInf, 0);
  sr.SetVar(mlb, 175, fInf);
  sr.SetVar(deepttag, -2, 1);
  sr.SetVar(tfttag, -2, 1);
  sr.SetVar(met, 250, fInf);
  sr.SetMETBins({250, 350, 450, 600, 1500});
  SRvec.emplace_back(sr);


  sr.SetName("srE0");
  sr.SetDetailName("ge4j_tmod0to10_mlb0to175");
  sr.SetVar(njet, 4, fInf);
  sr.SetVar(tmod, 0, 10);
  sr.SetVar(mlb, 0, 175);
  sr.SetVar(deepttag, -2, 1);
  sr.SetVar(tfttag, -2, 1);
  sr.SetVar(met, 450, fInf);
  sr.SetMETBins({450, 600, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srE1");
  sr.SetDetailName("ge4j_tmod0to10_mlb0to175_nottag");
  sr.SetVar(njet, 4, fInf);
  sr.SetVar(tmod, 0, 10);
  sr.SetVar(mlb, 0, 175);
  sr.SetVar(deepttag, -2, wpDeepTop);
  sr.SetVar(tfttag, -2, wpResTop);
  sr.SetVar(met, 250, 450);
  sr.SetMETBins({250, 350, 450});
  SRvec.emplace_back(sr);

  sr.SetName("srE2");
  sr.SetDetailName("ge4j_tmod0to10_mlb0to175_deepttag");
  sr.SetVar(deepttag, wpDeepTop, 1);
  sr.SetVar(tfttag, -2, 1);
  sr.SetVar(met, 250, 450);
  sr.SetMETBins({250, 350, 450});
  SRvec.emplace_back(sr);

  sr.SetName("srE3");
  sr.SetDetailName("ge4j_tmod0to10_mlb0to175_tfttag");
  sr.SetVar(deepttag, -2, wpDeepTop);
  sr.SetVar(tfttag, wpResTop, 1);
  sr.SetVar(met, 250, 450);
  sr.SetMETBins({250, 350, 450});
  SRvec.emplace_back(sr);


  sr.SetName("srF");
  sr.SetDetailName("ge4j_tmod0to10_mlb175toInf");
  sr.SetVar(njet, 4, fInf);
  sr.SetVar(tmod, 0, 10);
  sr.SetVar(mlb, 175, fInf);
  sr.SetVar(deepttag, -2, 1);
  sr.SetVar(tfttag, -2, 1);
  sr.SetVar(met, 250, fInf);
  sr.SetMETBins({250, 350, 450, 1500});
  SRvec.emplace_back(sr);


  sr.SetName("srG0");
  sr.SetDetailName("ge4j_tmod10toInf_mlb0to175");
  sr.SetVar(njet, 4, fInf);
  sr.SetVar(tmod, 10, fInf);
  sr.SetVar(mlb, 0, 175);
  sr.SetVar(deepttag, -2, 1);
  sr.SetVar(tfttag, -2, 1);
  sr.SetVar(met, 450, fInf);
  sr.SetMETBins({450, 550, 750, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srG1");
  sr.SetDetailName("ge4j_tmod10toInf_mlb0to175_nottag");
  sr.SetVar(njet, 4, fInf);
  sr.SetVar(tmod, 10, fInf);
  sr.SetVar(mlb, 0, 175);
  sr.SetVar(deepttag, -2, wpDeepTop);
  sr.SetVar(tfttag, -2, wpResTop);
  sr.SetVar(met, 250, 450);
  sr.SetMETBins({250, 350, 450});
  SRvec.emplace_back(sr);

  sr.SetName("srG2");
  sr.SetDetailName("ge4j_tmod10toInf_mlb0to175_deepttag");
  sr.SetVar(deepttag, wpDeepTop, 1);
  sr.SetVar(tfttag, -2, 1);
  sr.SetVar(met, 250, 450);
  sr.SetMETBins({250, 350, 450});
  SRvec.emplace_back(sr);

  sr.SetName("srG3");
  sr.SetDetailName("ge4j_tmod10toInf_mlb0to175_tfttag");
  sr.SetVar(deepttag, -2, wpDeepTop);
  sr.SetVar(tfttag, wpResTop, 1);
  sr.SetVar(met, 250, 450);
  sr.SetMETBins({250, 350, 450});
  SRvec.emplace_back(sr);


  sr.SetName("srH");
  sr.SetDetailName("ge4j_tmod10toInf_mlb175toInf");
  sr.SetVar(njet, 4, fInf);
  sr.SetVar(tmod, 10, fInf);
  sr.SetVar(mlb, 175, fInf);
  sr.SetVar(deepttag, -2, 1);
  sr.SetVar(tfttag, -2, 1);
  sr.SetVar(met, 250, fInf);
  sr.SetMETBins({250, 500, 1500});
  SRvec.emplace_back(sr);

  // SRvec.clear();  // TEMPORARY!!

  // Compressed (corridor) region
  sr = srbase;
  sr.SetName("srI");
  sr.SetDetailName("ge5j_lpt0to150_j1notb");
  sr.SetVar(njet, 5, fInf);
  sr.SetVar(nbjet, 1, fInf);
  // sr.SetVar(lep1pt, 0, 150);
  // sr.SetVar(dphilmet, 0, 2.0);
  sr.SetVar(passlmetcor, 1, 2); // shape cut on lep1pt & dphilmet
  sr.SetVar(dphijmet, 0.5, 3.1416);
  sr.SetVar(j1passbtag, 0, 1);  // Require j1 not b-tagged
  sr.RemoveVar(tmod);
  sr.RemoveVar(mlb);
  sr.SetMETBins({250, 350, 450, 550, 750, 1500});
  SRvec.emplace_back(sr);

  // sr = srbase;
  sr.SetName("srJ");
  sr.SetDetailName("ge5j_ge1sb_passlmetcor_j1notb");
  sr.SetVar(njet, 3, fInf);
  sr.SetVar(nbjet, 0, fInf);
  sr.SetVar(nsbtag, 1, fInf);
  sr.SetVar(passlmetcor, 1, 2); // shape cut on lep1pt & dphilmet
  sr.SetVar(dphijmet, 0.5, 3.1416);
  sr.SetVar(j1passbtag, 0, 1);  // Require j1 not b-tagged
  sr.SetMETBins({250, 350, 450, 550, 750, 1500});
  SRvec.emplace_back(sr);

  return SRvec;
}

std::vector<SR> getStopControlRegionsNoBTagsRun2(const std::vector<SR>& SRvec = getStopSignalRegionsRun2()) {
  std::vector<SR> CRvec;

  for (SR cr : SRvec) {
    cr.SetName(cr.GetName().replace(0, 2, "cr0b"));
    cr.SetAllowDummyVars(1);
    if (cr.VarExists(ntbtag)) cr.RemoveVar(ntbtag);
    if (cr.VarExists(nbjet)) {
      cr.SetVar(nbjet, 0, 1);
    } else if (cr.VarExists(nbtag)) {
      cr.SetVar(nbtag, 0, 1);
    }
    if (cr.VarExists(tmod) && cr.GetLowerBound(tmod) < 10)
      cr.SetVar(nbjet, 0, 1);

    if (cr.VarExists(tfttag))   cr.RemoveVar(tfttag);
    if (cr.VarExists(deepttag)) cr.RemoveVar(deepttag);

    if (cr.VarExists(nsbtag)) cr.SetVar(nsbtag, 0, 1);
    if (cr.VarExists(mlb)) cr.ReplaceVar(mlb, mlb_0b);
    CRvec.emplace_back(cr);
  }

  return CRvec;
}

std::vector<SR> getStopControlRegionsDileptonRun2(const std::vector<SR>& SRvec = getStopSignalRegionsRun2()) {
  std::vector<SR> CRvec;

  for (SR cr : SRvec) {
    cr.SetName(cr.GetName().replace(0, 2, "cr2l"));
    cr.SetAllowDummyVars(1);
    cr.ReplaceVar(met, met_rl);
    cr.ReplaceVar(mt, mt_rl);
    cr.ReplaceVar(dphijmet, dphijmet_rl);
    cr.ReplaceVar(nlep, nlep_rl);
    cr.RemoveVar(passvetos);
    if (cr.VarExists(tmod)) cr.ReplaceVar(tmod, tmod_rl);
    if (cr.VarExists(dphilmet)) cr.ReplaceVar(dphilmet, dphilmet_rl);
    if (cr.VarExists(passlmetcor)) cr.ReplaceVar(passlmetcor, passlmet_rl);
    cr.SetVar(nvlep, 2, fInf);
    cr.SetVar(nlep_rl, 2, fInf);
    cr.SetVar(mt2_ll, 0, 100);  // to avoid overlap with the stop-2l SR
    CRvec.emplace_back(cr);
  }

  return CRvec;
}

std::vector<SR> getStopInclusiveRegionsRun2() {

  SR srbase;
  srbase.SetAllowDummyVars(1);
  srbase.SetName("srbase");
  srbase.SetVar(mt, 150, fInf);
  srbase.SetVar(met, 250, fInf);
  srbase.SetVar(nlep, 1, 2);
  srbase.SetVar(nvlep, 1, 2);
  srbase.SetVar(passvetos, 1, 2);
  srbase.SetVar(njet, 2, fInf);
  srbase.SetVar(nbjet, 1, fInf);
  srbase.SetVar(mlb, 0, fInf);
  srbase.SetVar(tmod, -fInf, fInf);
  srbase.SetVar(dphijmet, 0.8, 3.1416);
  srbase.SetMETBins({0, 250, 350, 450, 550, 650, 800, 1500});

  SR sr;
  std::vector<SR> SRvec;

  SRvec.emplace_back(srbase);

  // The actual base region that is the sum of all SRs
  sr = srbase;
  sr.SetName("srincl0");
  sr.SetVar(njettmod, 1, 2);
  sr.ReplaceVar(nbjet, nbtag);
  SRvec.emplace_back(sr);

  // Inclusive for A+B
  sr.SetName("srincl1");
  sr.SetVar(njet, 2, 4);
  SRvec.emplace_back(sr);

  // Inclusive for C+D
  sr.SetName("srincl2");
  sr.SetVar(njet, 4, fInf);
  sr.SetVar(tmod, -fInf, 0);
  SRvec.emplace_back(sr);

  // Inclusive for E+F
  sr.SetName("srincl3");
  sr.SetVar(tmod, 0, 10);
  SRvec.emplace_back(sr);

  // Inclusive for G+H
  sr.SetName("srincl4");
  sr.SetVar(tmod, 10, fInf);
  SRvec.emplace_back(sr);

  // Inclusive for C-H (geq 4j )
  sr.SetName("srincl4J");
  sr.SetVar(tmod, -fInf, fInf);
  sr.SetVar(njet, 4, fInf);
  SRvec.emplace_back(sr);

  // Inclusive for E+G
  sr.SetName("srincl5");
  sr.SetVar(tmod, 0, fInf);
  sr.SetVar(mlb, 0, 175);
  SRvec.emplace_back(sr);

  // Inclusive for E+G, low MET
  sr.SetName("srincl5L");
  sr.SetVar(met, 250, 450);
  SRvec.emplace_back(sr);
  sr.SetVar(met, 250, fInf);

  // Inclusive for E+G, low MET + not merge tagged
  sr.SetName("srincl5LR");
  sr.SetVar(met, 250, 450);
  sr.SetVar(deepttag, -2, wpDeepTop);
  SRvec.emplace_back(sr);
  sr.SetVar(met, 250, fInf);
  sr.RemoveVar(deepttag);

  // Inclusive for F+H
  sr.SetName("srincl6");
  sr.SetVar(mlb, 175, fInf);
  SRvec.emplace_back(sr);

  // Inclusive for A+E+G, low MET, (deepAK8 tags)
  sr.SetName("srincl7L");
  sr.SetVar(njettmod, 1, 2);
  sr.SetVar(njet, 2, fInf);
  sr.SetVar(mlb, 0, 175);
  sr.SetVar(tmod, 0, fInf);
  sr.SetVar(met, 250, 450);
  SRvec.emplace_back(sr);

  // Inclusive for A, low MET, (deepAK8 tags)
  sr.SetName("srinclA");
  sr.SetVar(njet, 2, 4);
  sr.SetVar(mlb, 0, 175);
  sr.SetVar(tmod, 10, fInf);
  sr.SetVar(met, 250, 600);
  SRvec.emplace_back(sr);

  sr.SetName("srinclM");
  sr.SetDetailName("ge2j_met250to450_mtag");
  sr.SetVar(met, 250, 450);
  sr.SetVar(njet, 2, fInf);
  sr.SetVar(mlb, 0, fInf);
  sr.SetVar(tmod, -fInf, fInf);
  sr.SetVar(deepttag, wpDeepTop, 1.1);
  SRvec.emplace_back(sr);
  sr.SetName("srinclMI");
  sr.SetDetailName("ge2j_met250to450_invmtag");
  sr.SetVar(deepttag, -2, wpDeepTop);
  SRvec.emplace_back(sr);
  sr.SetVar(met, 250, fInf);
  sr.RemoveVar(deepttag);

  sr.SetName("srinclR");
  sr.SetDetailName("ge2j_met250to450_rtag");
  sr.SetVar(met, 250, 450);
  sr.SetVar(njet, 4, fInf);
  sr.SetVar(tfttag, wpResTop, 1.1);
  SRvec.emplace_back(sr);
  sr.SetName("srinclRI");
  sr.SetDetailName("ge2j_met250to450_invrtag");
  sr.SetVar(deepttag, -2, wpResTop);
  SRvec.emplace_back(sr);
  sr.SetVar(met, 250, fInf);
  sr.RemoveVar(tfttag);

  // Test region at MET sideband for top taggers
  sr = srbase;
  sr.SetName("srsbmet");  // MET sideband
  sr.SetDetailName("ge2j_met150to250");
  sr.SetVar(met, 150, 250);
  sr.SetVar(passlep1pt, 1, 2);
  sr.SetMETBins({150, 200, 250});
  SRvec.emplace_back(sr);

  sr.SetName("srsbmet2");  // MET sideband
  sr.SetDetailName("ge4j_met150to250_mt0toInf");
  sr.SetVar(met, 150, 250);
  sr.SetVar(mt, 0, fInf);
  sr.SetVar(njet, 4, fInf);
  sr.SetVar(passlep1pt, 1, 2);
  sr.SetMETBins({150, 200, 250});
  SRvec.emplace_back(sr);

  sr.SetName("srsbmet3");  // MET sideband
  sr.SetDetailName("ge3j_met150to250_mt0toInf");
  sr.SetVar(njet, 3, fInf);
  sr.SetVar(mt, 100, 250);
  SRvec.emplace_back(sr);

  sr = srbase;
  sr.SetName("srsbfmet4");  // MET sideband
  sr.SetDetailName("ge2j_met125toInf_mt50toInf_dphijmet0toPi");
  sr.SetVar(met, 125, fInf);
  sr.SetVar(mt, 50, fInf);
  sr.SetVar(dphijmet, 0, 3.1416);
  sr.SetVar(passlep1pt, 1, 2);
  sr.SetMETBins({125, 200, 250});
  SRvec.emplace_back(sr);

  sr = srbase;
  sr.SetName("srsbmt");  // MT sideband
  sr.SetDetailName("ge2j_met250toInf_mt0to150");
  sr.SetVar(mt, 0, 150);
  sr.SetVar(njet, 2, fInf);
  sr.SetVar(passlep1pt, 1, 2);
  sr.SetMETBins({150, 200, 250, 300, 400, 500, 800, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srsbmt2");  // MT sideband
  sr.SetDetailName("ge2j_met150toInf_mt0to150");
  sr.SetVar(met, 150, fInf);
  sr.SetMETBins({250, 300, 400, 500, 800, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srsbmlb");  // MT sideband
  sr.SetDetailName("ge4j_met150toInf_mt0to150");
  sr.SetVar(mt, 0, fInf);
  sr.SetVar(mlb, 175, fInf);
  sr.SetMETBins({250, 300, 400, 500, 800, 1500});
  SRvec.emplace_back(sr);
  sr.SetVar(mlb, 0, fInf);

  sr = srbase;
  sr.SetName("srsbfull");  // fullband
  sr.SetDetailName("ge2j_met150toInf_mt0toInf");
  sr.SetVar(met, 150, fInf);
  sr.SetVar(mt, 0, fInf);
  sr.SetVar(njet, 2, fInf);
  sr.SetVar(passlep1pt, 1, 2);
  sr.SetMETBins({150, 200, 250, 300, 400, 500, 800, 1500});
  SRvec.emplace_back(sr);

  sr.SetName("srsbfmt");
  sr.SetDetailName("ge2j_met250toInf_mt0toInf");
  sr.SetVar(njet, 2, fInf);
  sr.SetVar(met, 250, fInf);
  sr.SetVar(mt, 0, fInf);
  SRvec.emplace_back(sr);

  sr.SetName("srsbfmt2");
  sr.SetDetailName("ge2j_met250toInf_mt0toInf");
  sr.RemoveVar(passvetos);
  SRvec.emplace_back(sr);

  sr.SetName("srsbfmt4");
  sr.SetDetailName("ge4j_met250toInf_mt0toInf");
  sr.SetVar(njet, 4, fInf);
  sr.SetVar(passvetos, 1, 2);
  SRvec.emplace_back(sr);

  sr = srbase;
  sr.SetName("srsbfdphi");
  sr.SetDetailName("ge2j_met250toInf_dphi0topi");
  sr.SetVar(dphijmet, 0, 3.1416);
  SRvec.emplace_back(sr);

  sr = srbase;
  sr.SetName("srsbjfsb");
  sr.SetDetailName("ge5j_ge1sb_passlmetcor_j1notb");
  sr.SetVar(njet, 3, fInf);
  sr.SetVar(nbjet, 0, fInf);
  // sr.tVar(nsbtag, 0, fInf);
  sr.SetVar(passlmetcor, 1, 2); // shape cut on lep1pt & dphilmet
  sr.SetVar(dphijmet, 0.5, 3.1416);
  sr.SetVar(j1passbtag, 0, 1);  // Require j1 not b-tagged
  sr.SetMETBins({250, 350, 450, 550, 750, 1500});
  SRvec.emplace_back(sr);

  return SRvec;
}

std::vector<SR> getStopSignalRegionsCorridorRun2() {

  SR srbase;
  srbase.SetAllowDummyVars(1);
  srbase.SetName("srbasecor");
  srbase.SetVar(mt, 150, fInf);
  srbase.SetVar(met, 250, fInf);
  srbase.SetVar(nlep, 1, 2);
  srbase.SetVar(nvlep, 1, 2);
  srbase.SetVar(passvetos, 1, 2);
  srbase.SetVar(njet, 5, fInf);
  // srbase.SetVar(nbjet, 0, fInf);
  srbase.SetVar(dphijmet, 0.5, 3.1416);
  srbase.SetMETBins({0, 250, 350, 450, 550, 650, 800, 1500});

  std::vector<SR> SRvec;
  // SRvec.emplace_back(srbase);

  SR sr;
  sr = srbase;
  sr.SetName("srI");
  sr.SetDetailName("ge5j_lpt0to150_j1notb");
  sr.SetVar(njet, 5, fInf);
  sr.SetVar(nbjet, 1, fInf);
  // sr.SetVar(nsbtag, 1, fInf);
  sr.SetVar(passlmetcor, 1, 2); // shape cut on lep1pt & dphilmet
  sr.SetVar(dphijmet, 0.5, 3.1416);
  sr.SetVar(j1passbtag, 0, 1);  // Require j1 not b-tagged
  sr.SetMETBins({250, 350, 450, 550, 750, 1500});
  SRvec.emplace_back(sr);

  sr = srbase;
  sr.SetName("srJ0");
  sr.SetDetailName("ge5j_ge1sb_passlmetcor_j1notb");
  sr.SetVar(njet, 5, fInf);
  sr.SetVar(nbjet, 0, fInf);
  sr.SetVar(nsbtag, 1, fInf);
  sr.SetVar(passlmetcor, 1, 2); // shape cut on lep1pt & dphilmet
  sr.SetVar(dphijmet, 0.5, 3.1416);
  sr.SetVar(j1passbtag, 0, 1);  // Require j1 not b-tagged
  sr.SetMETBins({250, 350, 450, 550, 1500});
  SRvec.emplace_back(sr);

  // sr = srbase;
  // sr.SetName("srJ1");
  // sr.SetDetailName("ge5j_ge1sb_passlmetcor_j1notb");
  // sr.SetVar(njet, 4, fInf);
  // sr.SetVar(nbjet, 0, fInf);
  // sr.SetVar(nsbtag, 1, fInf);
  // sr.SetVar(passlmetcor, 1, 2); // shape cut on lep1pt & dphilmet
  // sr.SetVar(dphijmet, 0.5, 3.1416);
  // sr.SetVar(j1passbtag, 0, 1);  // Require j1 not b-tagged
  // sr.SetMETBins({250, 350, 450, 550, 750, 1500});
  // SRvec.emplace_back(sr);

  // sr = srbase;
  // sr.SetName("srJ2");
  // sr.SetDetailName("ge5j_ge1sb_passlmetcor_j1notb");
  // sr.SetVar(njet, 3, 4);
  // sr.SetVar(nbjet, 0, fInf);
  // sr.SetVar(nsbtag, 1, fInf);
  // sr.SetVar(passlmetcor, 1, 2); // shape cut on lep1pt & dphilmet
  // sr.SetVar(dphijmet, 0.5, 3.1416);
  // sr.SetVar(j1passbtag, 0, 1);  // Require j1 not b-tagged
  // sr.SetMETBins({250, 350, 450, 550, 750, 1500});
  // SRvec.emplace_back(sr);

  sr = srbase;
  sr.SetName("srJ3");
  sr.SetDetailName("ge5j_ge1sb_passlmetcor_j1notb");
  sr.SetVar(njet, 3, fInf);
  sr.SetVar(nbjet, 0, fInf);
  sr.SetVar(nsbtag, 1, fInf);
  sr.SetVar(passlmetcor, 1, 2); // shape cut on lep1pt & dphilmet
  sr.SetVar(dphijmet, 0.5, 3.1416);
  sr.SetVar(j1passbtag, 0, 1);  // Require j1 not b-tagged
  sr.SetMETBins({250, 350, 450, 550, 750, 1500});
  SRvec.emplace_back(sr);

  // sr = srbase;
  // sr.SetName("srJ4");
  // sr.SetDetailName("ge5j_ge1sb_passlmetcor_j1notb");
  // sr.SetVar(njet, 3, fInf);
  // sr.SetVar(nbjet, 1, fInf);
  // sr.SetVar(nsbtag, 1, fInf);
  // sr.SetVar(passlmetcor, 1, 2); // shape cut on lep1pt & dphilmet
  // sr.SetVar(dphijmet, 0.5, 3.1416);
  // sr.SetVar(j1passbtag, 0, 1);  // Require j1 not b-tagged
  // sr.SetMETBins({250, 350, 450, 550, 750, 1500});
  // SRvec.emplace_back(sr);

  // sr = srbase;
  // sr.SetName("srJ5");
  // sr.SetDetailName("ge5j_ge1sb_passlmetcor_j1notb");
  // sr.SetVar(njet, 3, fInf);
  // sr.SetVar(nbjet, 0, fInf);
  // sr.SetVar(nsbtag, 1, fInf);
  // sr.SetVar(passlmetcor, 1, 2); // shape cut on lep1pt & dphilmet
  // sr.SetVar(dphijmet, 0.5, 3.1416);
  // sr.SetVar(j1passbtag, 0, 1);  // Require j1 not b-tagged
  // sr.SetMETBins({350, 450, 550, 750, 1500});
  // SRvec.emplace_back(sr);

  return SRvec;
}

std::vector<SR> getStopControlRegionsNoBTagsCorridorRun2() {
  std::vector<SR> CRvec;
  std::vector<SR> SRvec = getStopSignalRegionsCorridorRun2();

  for (SR cr : SRvec) {
    cr.SetName(cr.GetName().replace(0, 2, "cr0b"));
    cr.SetAllowDummyVars(1);
    if (cr.VarExists(nbjet))  cr.RemoveVar(nbjet);
    if (cr.VarExists(ntbtag)) cr.RemoveVar(ntbtag);
    if (cr.VarExists(nsbtag)) cr.SetVar(nsbtag, 0, 1);
    cr.SetVar(nbjet, 0, 1);
    if (cr.VarExists(mlb)) cr.ReplaceVar(mlb, mlb_0b);
    CRvec.emplace_back(cr);
  }

  return CRvec;
}

std::vector<SR> getStopControlRegionsDileptonCorridorRun2() {
  std::vector<SR> CRvec;
  std::vector<SR> SRvec = getStopSignalRegionsCorridorRun2();

  for (SR cr : SRvec) {
    cr.SetName(cr.GetName().replace(0, 2, "cr2l"));
    cr.SetAllowDummyVars(1);
    cr.ReplaceVar(met, met_rl);
    cr.ReplaceVar(mt, mt_rl);
    cr.ReplaceVar(dphijmet, dphijmet_rl);
    cr.ReplaceVar(nlep, nlep_rl);
    // cr.ReplaceVar(tmod, tmod_rl);
    cr.RemoveVar(passvetos);
    if (cr.VarExists(dphilmet)) cr.ReplaceVar(dphilmet, dphilmet_rl);
    if (cr.VarExists(passlmetcor)) cr.ReplaceVar(passlmetcor, passlmet_rl);
    cr.SetVar(nvlep, 2, fInf);
    cr.SetVar(nlep_rl, 2, fInf);
    cr.SetVar(mt2_ll, 0, 100);  // to avoid overlap with the stop-2l SR
    CRvec.emplace_back(cr);
  }

  return CRvec;
}

