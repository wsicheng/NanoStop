#ifndef STOPREGIONS_H
#define STOPREGIONS_H

#include <string>
#include "TTree.h"
#include "SR.h"

// ivar:            0,          1,           2,           3,           4,          5,          6,          7,          8,          9,
enum Vars {       met,        mt,         nlep,       nvlep,   passvetos,       njet,      nbjet,      nbtag,     ntbtag,        mlb,
                 tmod,  dphijmet,     deepttag,     resttag,     bdtttag,     tfttag,   njettmod,
               met_rl,     mt_rl,      tmod_rl,     nlep_rl, dphijmet_rl,     mt2_ll,        mll,     mlb_0b,
           j1passbtag,  dphilmet,  dphilmet_rl, passlmetcor, passlmet_rl,     nsbtag,
               metphi,    lep1pt,       lep2pt,     lep1eta,     lep2eta,     jet1pt,     jet2pt,    jet1eta,    jet2eta,
             nvtx, ht, nak8jets, chi2, binttag, leadbpt, mllbbmet, mtttbar, ptttbar, ptll, ptbb, passlep1pt,
                  mht, mhtphi, met_rs, metphi_rs, mt_rs, dphijmet_rs, dphilmet_rs, metphi_rl, nsblep,
                  // rmet, rmetphi, cmet, cmetphi, mt_cmet, dphilcmet,
                nvars };


// RunII analysis selections
std::vector<SR> getStopSignalRegionsRun2();
std::vector<SR> getStopInclusiveRegionsRun2();
std::vector<SR> getStopCrosscheckRegionsEMuRun2();
std::vector<SR> getStopControlRegionsNoBTagsRun2(const std::vector<SR>& SRvec);
std::vector<SR> getStopControlRegionsDileptonRun2(const std::vector<SR>& SRvec);

std::vector<SR> getStopSignalRegionsCorridorRun2();
std::vector<SR> getStopControlRegionsNoBTagsCorridorRun2();
std::vector<SR> getStopControlRegionsDileptonCorridorRun2();

// Old helper functions
std::vector<SR> getStopControlRegionsNoBTags(std::vector<SR>&& SRvec);
std::vector<SR> getStopControlRegionsDilepton(std::vector<SR>&& SRvec);

// Old 2016 Analysis search bins
std::vector<SR> getStopSignalRegionsTopological();
std::vector<SR> getStopControlRegionsNoBTagsTopological();
std::vector<SR> getStopControlRegionsDileptonTopological();

#endif // STOPREGIONS_H
