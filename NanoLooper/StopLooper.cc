// C++
#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <set>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include <numeric>

// ROOT
#include "TDirectory.h"
#include "TTreeCache.h"
#include "Math/VectorUtil.h"
#include "TVector2.h"
#include "TBenchmark.h"
#include "TMath.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TF1.h"

// NanoCORE
// #include "../NanoCORE/Tools/utils.h"
// #include "../NanoCORE/Tools/badEventFilter.h"
#include "../NanoCORE/Nano.h"
#include "../NanoCORE/Config.h"
#include "../NanoCORE/MetSelections.h"
#include "../NanoCORE/goodrun.h"
#include "../NanoCORE/dorky.h"

// #include "../StopCORE/StopTree.h"
// #include "../StopCORE/TopTagger/ResolvedTopMVA.h"
// #include "../StopCORE/METCorr/METCorrectionHandler.h"

// Stop class
#include "StopSelections.h"
#include "StopRegions.h"
#include "StopLooper.h"
#include "Utilities.h"
#include "SR.h"

using namespace std;
using namespace tas;
using namespace utils;

class SR;

const bool verbose = false;
// turn on to apply btag sf - take from files defined in eventWeight_btagSF.cc
const bool applyBtagSFfromFiles = false; // default false
// turn on to apply lepton sf to central value - reread from files
const bool applyLeptonSFfromFiles = false; // default false
// turn on to apply json file to data
const bool applyGoodRunList = true;
// apply the MET resolution correction to the MET variables, required running
const bool applyMETResCorrection = true;
// turn on to enable plots of metbins with systematic variations applied. will only do variations for applied weights
const bool doSystVariations = false;
// turn on to enable plots of metbins with different gen classifications
const bool doGenClassification = true;
// turn on to apply Nvtx reweighting to MC / data2016
const bool doNvtxReweight = false;
// turn on top tagging studies, off for baby ver < 25
const bool doTopTagging = true;
// veto 2018 events with an electron land in the HEM region
const bool doHEMElectronVeto = true;
// veto 2018 events with an AK4 jet land in the HEM region
const bool doHEMJetVeto = true;
// re-run resolved top MVA locally
const bool runResTopMVA = false;
// run the MET resolution correction (and store in separate branches)
const bool runMETResCorrection = false;
// only produce yield histos
const bool runYieldsOnly = false;
// only running selected signal points to speed up
const bool runFullSignalScan = false;
// debug symbol, for printing exact event kinematics that passes
const bool printPassedEvents = false;
// switch to use the separate fine scan points in the corridor or combine them
const bool combineCorridorScans = false;
// fill the distribution of event weights
const bool fillWeights = false;
// apply the HT5/HT cut
const bool doHT5cut = true;

// some global helper variables to be used in member functions
int datayear = -1;
string samplever;
string sigtype;
float babyver;

const float fInf = std::numeric_limits<float>::max();

const float kSMSMassStep = 25;
const vector<float> mStopBins = []() { vector<float> bins; for (float m = 137.5; m < 2137.5; m += kSMSMassStep) bins.push_back(m); return bins; } ();
const vector<float> mLSPBins  = []() { vector<float> bins; for (float m = -12.5; m < 1437.5; m += kSMSMassStep) bins.push_back(m); return bins; } ();

// Helper function to pick out signal events
// auto checkMassPt = [&](double mstop, double mlsp) { return (mass_stop() == mstop) && (mass_lsp() == mlsp); };
auto checkMassPt = [&](double mstop, double mlsp) { return (false); };

std::ofstream ofile;

void StopLooper::SetSignalRegions() {


  SRVec = getStopSignalRegionsRun2();

  // Adding the inclusive regions
  if (!runYieldsOnly) {
    auto inclSRvec = getStopInclusiveRegionsRun2();
    SRVec.insert(SRVec.begin()+1, inclSRvec.begin()+1, inclSRvec.end());
    CRemuVec = getStopCrosscheckRegionsEMuRun2();
  }

  CR0bVec = getStopControlRegionsNoBTagsRun2(SRVec);
  CR2lVec = getStopControlRegionsDileptonRun2(SRVec);

  values_.resize(nvars, NAN);

  if (verbose) {
    cout << "SRVec.size = " << SRVec.size() << ", including the following:" << endl;
    for (auto it = SRVec.begin(); it != SRVec.end(); ++it) {
      cout << it-SRVec.begin() << "  " << it->GetName() << ":  " << it->GetDetailName() << endl;
    }
  }

  auto createRangesHists = [&] (vector<SR>& srvec) {
    for (auto& sr : srvec) {
      vector<string> vars = sr.GetListOfVariables();
      TDirectory * dir = (TDirectory*) outfile_->Get((sr.GetName() + "/ranges").c_str());
      if (dir == 0) dir = outfile_->mkdir((sr.GetName() + "/ranges").c_str());
      dir->cd("ranges");
      for (auto& var : vars) {
        plot1d("h_"+var+"_"+"LOW",  1, sr.GetLowerBound(var), sr.histMap, "", 1, 0, 2);
        plot1d("h_"+var+"_"+"HI",   1, sr.GetUpperBound(var), sr.histMap, "", 1, 0, 2);
      }
      if (sr.GetNMETBins() > 0) {
        plot1d("h_metbins", -1, 0, sr.histMap, (sr.GetName()+":"+sr.GetDetailName()+";E^{miss}_{T} [GeV]").c_str(), sr.GetNMETBins(), sr.GetMETBinsPtr());
      }
      if (false) {
        // temporary for the tests of fullsim corridor
        plot3d("hSMS_metbins", -1, 0, 0, 0, sr.histMap, (sr.GetName()+":"+sr.GetDetailName()+";E^{miss}_{T} [GeV];M_{stop};M_{LSP}").c_str(),
               sr.GetNMETBins(), sr.GetMETBinsPtr(), mStopBins.size()-1, mStopBins.data(), mLSPBins.size()-1, mLSPBins.data());
      }
    }
  };

  createRangesHists(SRVec);
  createRangesHists(CR0bVec);
  createRangesHists(CR2lVec);
  createRangesHists(CRemuVec);

  testVec.emplace_back("testGeneral");
  // testVec.emplace_back("testTopTagging");
  testVec.emplace_back("testCutflow");
}

bool StopLooper::PassingHLTriggers(const int type) {
  // if (type == 1) {
  //   switch (year_) {
  //     case 2016:
  //       return ( (HLT_MET110_MHT110() || HLT_MET120_MHT120() || HLT_MET()) ||
  //                (abs(lep1_pdgid()) == 11 && HLT_SingleEl()) || (abs(lep1_pdgid()) == 13 && HLT_SingleMu()) );
  //     case 2017:
  //     case 2018:
  //     default:
  //       return ( HLT_MET_MHT() || (abs(lep1_pdgid()) == 11 && HLT_SingleEl()) || (abs(lep1_pdgid()) == 13 && HLT_SingleMu()) );
  //   }
  // } else if (type == 2) {
  //   switch (year_) {
  //     case 2016:
  //       return ( (HLT_MET() || HLT_MET110_MHT110() || HLT_MET120_MHT120()) || (HLT_SingleEl() && (abs(lep1_pdgid())==11 || abs(lep2_pdgid())==11)) ||
  //                (HLT_SingleMu() && (abs(lep1_pdgid())==13 || abs(lep2_pdgid())==13)) );
  //     case 2017:
  //     case 2018:
  //     default:
  //       return ( HLT_MET_MHT() || (HLT_SingleEl() && (abs(lep1_pdgid())==11 || abs(lep2_pdgid())==11)) ||
  //                (HLT_SingleMu() && (abs(lep1_pdgid())==13 || abs(lep2_pdgid())==13)) );
  //   }
  // } else if (type == 3) {
  //   int dilepid = abs(lep1_pdgid()*lep2_pdgid());
  //   switch (year_) {
  //     case 2016:
  //       return ( (HLT_MET() || HLT_MET110_MHT110() || HLT_MET120_MHT120()) ||
  //                (dilepid == 121 && HLT_DiEl()) || (dilepid == 169 && HLT_DiMu()) || (dilepid == 143 && HLT_MuE()) ||
  //                (HLT_SingleEl() && (abs(lep1_pdgid())==11 || abs(lep2_pdgid())==11)) ||
  //                (HLT_SingleMu() && (abs(lep1_pdgid())==13 || abs(lep2_pdgid())==13)) );
  //     case 2017:
  //     case 2018:
  //     default:
  //       return ( HLT_MET_MHT() || (dilepid == 121 && HLT_DiEl()) || (dilepid == 169 && HLT_DiMu()) || (dilepid == 143 && HLT_MuE()) ||
  //                (HLT_SingleEl() && (abs(lep1_pdgid())==11 || abs(lep2_pdgid())==11)) ||
  //                (HLT_SingleMu() && (abs(lep1_pdgid())==13 || abs(lep2_pdgid())==13)) );
  //   }
  // }
  // return false;
  return true;
}

void StopLooper::looper(TChain* chain, string samplestr, string output_dir, int jes_type) {

  // Benchmark
  TBenchmark *bmark = new TBenchmark();
  bmark->Start("benchmark");

  TString output_name = Form("%s/%s.root",output_dir.c_str(),samplestr.c_str());
  cout << "[StopLooper::looper] creating output file: " << output_name << endl;  outfile_ = new TFile(output_name.Data(),"RECREATE") ;
  cout << "Complied with C++ standard: " << __cplusplus << endl;

  if (printPassedEvents) ofile.open(output_dir+"/passEventList_"+samplestr+".txt");

  // if (runResTopMVA)
  //   resTopMVA = new ResolvedTopMVA("../StopCORE/TopTagger/resTop_xGBoost_v2.weights.xml", "BDT");

  outfile_ = new TFile(output_name.Data(), "RECREATE") ;

  // Combined 2016 (35.922/fb), 2017 (41.529/fb) and 2018 (59.740/fb) json,
  const char* json_file = "../NanoCORE/goodrun_files/Cert_271036-325175_13TeV_Combined161718_JSON_snt.txt";
  if (applyGoodRunList) {
    cout << "Loading json file: " << json_file << endl;
    set_goodrun_file(json_file);
  }

  // Determine datayear from the sample name
  if (samplestr == "data_single_lepton_met") datayear = 16;
  else if (samplestr.find("data_2016") == 0) datayear = 2016;
  else if (samplestr.find("data_2017") == 0) datayear = 2017;
  else if (samplestr.find("data_2018") == 0) datayear = 2018;
  else datayear = -1;

  // TFile dummy( (output_dir+"/dummy.root").c_str(), "RECREATE" );
  SetSignalRegions();

  // MET resolution stuff
  // METCorrectionHandler metCorrector;

  int nDuplicates = 0;
  int nEvents = chain->GetEntries();
  unsigned int nEventsChain = nEvents;
  cout << "[StopLooper::looper] running on " << nEventsChain << " events" << endl;
  unsigned int nEventsTotal = 0;
  unsigned int nPassedTotal = 0;

  TObjArray *listOfFiles = chain->GetListOfFiles();
  TIter fileIter(listOfFiles);
  TFile *currentFile = 0;
  while ( (currentFile = (TFile*)fileIter.Next()) ) {

    TString fname = currentFile->GetTitle();
    TFile file( fname, "READ" );
    TTree *tree = (TTree*) file.Get("Events");
    TTreeCache::SetLearnEntries(10);
    tree->SetCacheSize(128*1024*1024);
    nt.Init(tree);

    // Use the first event to get dsname
    tree->LoadTree(0);
    nt.GetEntry(0);

    // Setup configs for sample dependent processes
    gconf.GetConfigsFromDatasetName(fname.Data());
    year_ = gconf.year;

    TString dsname = "";  // FIXME: where to find dataset info?

    // Find the stopbaby versions automatically from file path <-- FIXME: is this needed anymore?
    // if (int i = fname.Index("_v"); i >= 0) samplever = fname(i+1, min(fname.Index("/",i),5)); // include subversions
    // else if (fname.Contains("v24")) samplever = "v24";
    // else cout << "[looper] >> Cannot find the sample version!" << endl;
    // babyver = TString(samplever).ReplaceAll("v","").ReplaceAll("_",".").Atof();

    cout << "[StopLooper] >> Running on sample: " << dsname << endl;
    cout << "[StopLooper] >> Sample detected with year = " << year_ << " and version = " << samplever << " (" << babyver << ")" << endl;

    is_fastsim_ = fname.Contains("SMS") || fname.Contains("Signal");
    is_bkg_ = (!gconf.is_data && !is_fastsim_);
    // if (dsname.Contains("Fast")) is_fastsim_ = true;

    if (int i = dsname.Index("mStop-"); i >= 0) mstop_ = TString(dsname(i+6, dsname.Index("_",i))).Atof();
    if (int i = dsname.Index("mLSP-"); i >= 0) mlsp_ = TString(dsname(i+5, dsname.Index("_",i))).Atof();
    if (!is_fastsim_ && mlsp_ > 0 && mstop_ > 0)
      cout << "[StopLooper] >> Detect fullsim signal!! With mstop = " << mstop_ << ", mlsp = " << mlsp_ << endl;

    // Get event weight histogram from baby
    TH2D* h_sig_counter_nEvents = nullptr;
    if (is_fastsim_) h_sig_counter_nEvents = (TH2D*) file.Get("histNEvts");

    // TODO: // Setup the event weight calculator

    // evtWgt.verbose = true;
    // evtWgt.setDefaultSystematics(evtWgtInfo::stop_Run2, is_fastsim_);

    // evtWgt.apply_HEMveto_el_sf = doHEMElectronVeto;
    // evtWgt.apply_HEMveto_jet_sf = doHEMJetVeto;
    // if (babyver > 31.7 && dsname.Contains("T2tt"))
    //   evtWgt.combineCorridorScans = false;

    // evtWgt.Setup(samplestr, year_, doSystVariations, applyBtagSFfromFiles, applyLeptonSFfromFiles);

    // evtWgt.getCounterHistogramFromBaby(&file);
    // // Extra file weight for extension dataset, should move these code to other places
    // if (year_ == 2016 && samplever.find("v22") == 0)
    //   evtWgt.getExtSampleWeightSummer16v2(fname);
    // else if (year_ == 2016 && samplever.find("Summer16v3") != string::npos)
    //   evtWgt.getExtSampleWeightSummer16v3(fname);
    // else if (year_ == 2017)
    //   evtWgt.getExtSampleWeightFall17v2(fname);

    // evtWgt.getZSampleWeightFromCR3l(fname);

    // if (!gconf.is_data && runMETResCorrection) // setup MET resolution stuff
    //   metCorrector.setup(year_, to_string(year_), "../StopCORE/METCorr/METSFs");

    // dummy.cd();
    // Loop over Events in current file
    if (nEventsTotal >= nEventsChain) continue;
    unsigned int nEventsTree = tree->GetEntriesFast();
    for (unsigned int evt = 0; evt < nEventsTree; ++evt) {
      if (evt > 100) break;  // FIXME: debug
      // Read Tree
      if (nEventsTotal >= nEventsChain) continue;
      tree->LoadTree(evt);
      nt.GetEntry(evt);
      ++nEventsTotal;

      if ( gconf.is_data ) {
        if ( applyGoodRunList && !goodrun(run(), luminosityBlock()) ) continue;
        duplicate_removal::DorkyEventIdentifier id(run(), event(), luminosityBlock());
        if ( is_duplicate(id) ) {
          ++nDuplicates;
          continue;
        }
      }

      // fillEfficiencyHistos(testVec[0], "filters");

      // Apply met filters
      bool passMETfilt = passesMETfilters(gconf.is_data);
      if (!passMETfilt) continue;

      // Require at least 1 good vertex
      if (PV_npvsGood() < 1) continue;

      // Fill tirgger efficiency histos after the MET filters are applied
      // fillEfficiencyHistos(testVec[0], "triggers");

      // TODO: // Only consider events with nupt < 200 for the inclusive WNJetsToLNu samples
      // if (dsname.BeginsWith("/W") && dsname.Contains("JetsToLNu") && !dsname.Contains("NuPt-200") && nupt() > 200) continue;

      // mstop_ = 1200; mlsp_ = 1; 
      if (is_fastsim_) {
        // plot2d("h2d_signal_masspts", mass_stop(), mass_lsp(), 1, SRVec.at(0).histMap, ";M_{stop} [GeV]; M_{lsp} [GeV]", 320, 100, 2100, 120, 0, 1500);
        plot2d("h2d_sigpts_ylds", mstop_, mlsp_, 1, SRVec.at(0).histMap, ";M_{stop} [GeV]; M_{lsp} [GeV]", mStopBins.size()-1, mStopBins.data(), mLSPBins.size()-1, mLSPBins.data());
        // int nEventsPoint = h_sig_counter_nEvents->GetBinContent(h_sig_counter_nEvents->FindBin(mstop_, mlsp_));
        // evtweight_ = kLumi * xsec() * 1000 / nEventsPoint;
        plot2d("h2d_sigpts_xsecwgtd", mstop_, mlsp_, evtweight_, SRVec.at(0).histMap, ";M_{stop} [GeV]; M_{lsp} [GeV]", 320, 100, 2100, 120, 0, 1500);
      } else if (mlsp_ > 0 && mstop_ > 0) {
        plot2d("h2d_sigpts_ylds", mstop_, mlsp_, 1, SRVec.at(0).histMap, ";M_{stop} [GeV]; M_{lsp} [GeV]", mStopBins.size()-1, mStopBins.data(), mLSPBins.size()-1, mLSPBins.data());
        // plot2d("h2d_sigpts_xsecwgtd", mstop_, mlsp_, scale1fb(), SRVec.at(0).histMap, ";M_{stop} [GeV]; M_{lsp} [GeV]", 320, 100, 2100, 120, 0, 1500);
      }

      ++nPassedTotal;

      // TODO: // Calculate event weight
      // evtWgt.resetEvent(); // full event weights only get calculated if the event get selected for a SR

      // Simple weight with scale1fb only
      // if (is_bkg_) evtweight_ = kLumi * scale1fb();

      evtweight_ = 1;  // FIXME: find out how to get correct event weight

      // Plot nvtxs on the base selection of stopbaby for reweighting purpose
      plot1d("h_nvtxs", PV_npvs(), 1, testVec[0].histMap, ";Number of primary vertices", 100, 1, 101);

      // Calculate leading top tagger variables for the event
      float lead_restopdisc = -1.1;
      float lead_tftopdisc = -0.1;
      float lead_deepdisc_top = -0.1;
      float lead_bindisc_top = -0.1;
      float lead_tftopdisc_jup = lead_tftopdisc;
      float lead_tftopdisc_jdown = lead_tftopdisc;
      if (doTopTagging) {
        for (size_t iak8 = 0; iak8 < nFatJet(); ++iak8) {
          float qcddisc = FatJet_deepTag_QCD()[iak8];
          float bindisc = FatJet_deepTag_TvsQCD()[iak8];
          if (bindisc > lead_bindisc_top) lead_bindisc_top = bindisc;
          float rawdisc = (bindisc < 1)? qcddisc*bindisc / (1 - bindisc) : 1;
          if (rawdisc > lead_deepdisc_top) lead_deepdisc_top = rawdisc;
        }

        // TODO: get the DeepResolved score <-- can be run online
        // lead_restopdisc = (topcands_disc().size())? topcands_disc()[0] : -1.1;
        // for (auto disc : tftops_disc()) {
        //   if (disc > lead_tftopdisc) lead_tftopdisc = disc;
        // }
        // for (auto disc : jup_tftops_disc()) {
        //   if (disc > lead_tftopdisc_jup) lead_tftopdisc_jup = disc;
        // }
        // for (auto disc : jdown_tftops_disc()) {
        //   if (disc > lead_tftopdisc_jdown) lead_tftopdisc_jdown = disc;
        // }
      }

      // TODO: // Number of soft b-tags exlcuded from the leptons
      // int nsboverlep = 0;
      // for (auto sb : softtags_p4()) {
      //   if (nvetoleps() > 0 && isCloseObject(sb, lep1_p4(), 0.4)) nsboverlep++;
      //   if (nvetoleps() > 1 && isCloseObject(sb, lep2_p4(), 0.4)) nsboverlep++;
      // }

      int nsoftbjets = 1;  // Temporary!

      // TODO: make the lepton objects first
      // HEM issue for 2018, that excess of electron happens at -4.7 < eta < -1.4, -1.6 < phi < -0.8
      if (doHEMElectronVeto && gconf.year == 2018 && gconf.is_data && run() >= 319077) {
        // if (abs(lep2_pdgid()) == 11 && lep2_p4().eta() < -1.4 && lep2_p4().phi() > -1.6 && lep2_p4().phi() < -0.8)
        //   continue; // veto the events
      }

      // TODO: contruct the jet objects
      if (doHEMJetVeto && year_ == 2018 && gconf.is_data && run() >= 319077) {
        bool hasHEMjet = false;
        // for (auto jet : ak4pfjets_p4()) {  // the jets are cleaned with leptons already
        //   if (jet.eta() < -1.4 && jet.phi() > -1.6 && jet.phi() < -0.8)
        //     hasHEMjet = true;
        // }
        // if (hasHEMjet) continue;
      }

      LorentzVector sumMHTp4(0,0,0,0);
      // if (babyver >= 31.2) {
      //   for (auto jet : ak4pfjets_p4()) {
      //     sumMHTp4 -= jet;
      //   }
      // } else {
      //   sumMHTp4 = LorentzVector(ak4_MHT_pt()*cos(ak4_MHT_phi()),ak4_MHT_pt()*sin(ak4_MHT_phi()),0,ak4_MHT_pt());
      // }
      // sumMHTp4 -= lep1_p4();
      // if (nvetoleps() > 1) sumMHTp4 -= lep2_p4();
      float MHT_pt = sumMHTp4.pt();
      float MHT_phi = sumMHTp4.phi();

      // float dphij1mht = (ngoodjets() > 0)? deltaPhi(ak4pfjets_p4().at(0).phi(), MHT_phi) : 9;
      // if (babyver >= 31.2 && ngoodjets() > 0)
      //   plot2d("h2d_ht5_dphijmht", dphij1mht, ak4_HTeta5()/ak4_HT(), 1, testVec[0].histMap, ";#Delta#phi(j1,MHT);HT5/HT", 16, 0, 3.2, 20, 1, 3);

      // if (doHT5cut && year_ >= 2017) {
      //   float HT5overHT = ak4_HTeta5() / ak4_HT();
      //   // if (dphij1mht < 5.3*HT5overHT - 4.78) continue;
      //   if (HT5overHT > 1.5) continue;
      // }

      // TODO: // MET resolution stuff: recalculate MET using METObject
      // METObject metobj;
      // if (!gonf.is_data && runMETResCorrection) {
      //   metobj.extras.met
      //       = metobj.extras.met_METup = metobj.extras.met_METdn
      //       = metobj.extras.met_JERup = metobj.extras.met_JERdn
      //       = metobj.extras.met_PUup = metobj.extras.met_PUdn
      //       = pfmet();
      //   metobj.extras.met_JECup = pfmet_jup();
      //   metobj.extras.met_JECdn = pfmet_jdown();
      //   metobj.extras.phi
      //       = metobj.extras.phi_METup = metobj.extras.phi_METdn
      //       = metobj.extras.phi_JERup = metobj.extras.phi_JERdn
      //       = metobj.extras.phi_PUup = metobj.extras.phi_PUdn
      //       = pfmet_phi();
      //   metobj.extras.phi_JECup = pfmet_phi_jup();
      //   metobj.extras.phi_JECdn = pfmet_phi_jdown();
      //   metCorrector.correctMET(genmet(), genmet_phi(), &metobj, is_fastsim_); // Last flag is for fast sim., but there are no separate factors, so it doesn't matter
      // }

      /*
      // Filling the variables for analysis
      values_.clear();
      values_.resize(nvars, NAN);

      values_[mht] = MHT_pt;
      values_[mhtphi] = MHT_phi;

      /// Common variables for all JES type
      values_[nlep] = ngoodleps();
      values_[nvlep] = nvetoleps();
      values_[lep1pt] = lep1_p4().pt();
      values_[passvetos] = PassTrackVeto() && PassTauVeto();
      values_[nvtx] = nvtxs();
      values_[nlep_rl] = (ngoodleps() == 1 && nvetoleps() >= 2 && lep2_p4().Pt() > 10)? 2 : ngoodleps();
      values_[mll] = (lep1_p4() + lep2_p4()).M();

      // For toptagging, add correct switch later
      values_[nak8jets] = ak8pfjets_deepdisc_top().size();
      // values_[resttag] = lead_restopdisc;
      values_[resttag] = lead_tftopdisc;
      values_[deepttag] = lead_deepdisc_top;
      values_[tfttag] = lead_tftopdisc;
      values_[bdtttag] = lead_restopdisc;
      values_[binttag] = lead_bindisc_top;
      values_[nsblep] = nsboverlep;

      // // Temporary to for no top-tagging result
      // values_[deepttag] = 0;
      // values_[tfttag] = 0;

      if (runResTopMVA) {
        // Prepare deep_cvsl vector
        vector<float> ak4pfjets_dcvsl;
        for (size_t j = 0; j < ak4pfjets_deepCSV().size(); ++j) {
          ak4pfjets_dcvsl.push_back(ak4pfjets_deepCSVc().at(j) / (ak4pfjets_deepCSVc().at(j) + ak4pfjets_deepCSVl().at(j)));
        }
        resTopMVA->setJetVecPtrs(&ak4pfjets_p4(), &ak4pfjets_deepCSV(), &ak4pfjets_dcvsl, &ak4pfjets_ptD(), &ak4pfjets_axis1(), &ak4pfjets_mult());
        std::vector<TopCand> topcands = resTopMVA->getTopCandidates(-1);
        values_[resttag] = (topcands.size() > 0)? topcands[0].disc : -1.1;
      }

      // Temporary test for top tagging efficiency
      // testTopTaggingEffficiency(testVec[1]);
      // testGenMatching(testVec[0]);

      /// Values only for hist filling or testing
      values_[chi2] = hadronic_top_chi2();
      values_[lep1eta] = lep1_p4().eta();
      values_[passlep1pt] = (abs(lep1_pdgid()) == 13 && lep1_p4().pt() > 40) || (abs(lep1_pdgid()) == 11 && lep1_p4().pt() > 45);

      for (int systype = 0; systype < ((doSystVariations && !gconf.is_data)? 5 : 1); ++systype) {
        string suffix = "";

        /// JES type dependent variables
        if (systype == 0) {
          values_[mt] = mt_met_lep();
          values_[met] = pfmet();
          values_[mlb] = Mlb_closestb();
          values_[tmod] = topnessMod();
          values_[njet] = ngoodjets();
          values_[nbjet] = ngoodbtags();
          values_[ntbtag] = ntightbtags();
          values_[tfttag] = lead_tftopdisc;
          values_[dphijmet] = mindphi_met_j1_j2();
          values_[dphilmet] = lep1_dphiMET();
          values_[j1passbtag] = (ngoodjets() > 0)? ak4pfjets_passMEDbtag().at(0) : 0;

          values_[jet1pt] = (ngoodjets() > 0)? ak4pfjets_p4().at(0).pt() : 0;
          values_[jet2pt] = (ngoodjets() > 1)? ak4pfjets_p4().at(1).pt() : 0;
          values_[jet1eta] = (ngoodjets() > 0)? ak4pfjets_p4().at(0).eta() : -9;
          values_[jet2eta] = (ngoodjets() > 1)? ak4pfjets_p4().at(1).eta() : -9;

          values_[metphi] = pfmet_phi();
          values_[nbtag] = nanalysisbtags();
          values_[nsbtag] = nsoftbjets;
          values_[leadbpt] = ak4pfjets_leadbtag_p4().pt();
          values_[mlb_0b] = (ak4pfjets_leadbtag_p4() + lep1_p4()).M();

          jestype_ = 0;
          // suffix = "_nominal";
        } else if (systype == 1) {
          if (!applyMETResCorrection) continue;
          if (samplever.find("v31") == 0) {  // available on and after v31_2
            values_[mt] = mt_met_lep_resup();
            values_[met] = pfmet_resup();
            values_[tmod] = topnessMod_resup();
            values_[dphijmet] = mindphi_met_j1_j2_resup();
            values_[dphilmet] = lep1_dphiMET_resup();
            values_[metphi] = pfmet_phi_resup();
          } else if (runMETResCorrection) {
            values_[met] = metobj.extras.met_METup;
            values_[metphi] = metobj.extras.phi_METup;
            values_[mt] = calculateMT(lep1_p4().pt(), lep1_p4().phi(), values_[met], values_[metphi]);
            values_[dphilmet] = deltaPhi(lep1_p4().phi(), values_[metphi]);
            values_[dphijmet] = (ngoodjets() > 1)? min2deltaPhi(values_[metphi], ak4pfjets_p4().at(0).phi(), ak4pfjets_p4().at(1).phi()) : -9;
          }

          suffix = "_metResUp";
        } else if (systype == 2) {
          if (!applyMETResCorrection) continue;
          if (samplever.find("v31") == 0) {  // available on and after v31_2
            values_[mt] = mt_met_lep_resdown();
            values_[met] = pfmet_resdown();
            values_[tmod] = topnessMod_resdown();
            values_[dphijmet] = mindphi_met_j1_j2_resdown();
            values_[dphilmet] = lep1_dphiMET_resdown();
            values_[metphi] = pfmet_phi_resdown();
          } else if (runMETResCorrection) {
            values_[met] = metobj.extras.met_METdn;
            values_[metphi] = metobj.extras.phi_METdn;
            values_[mt] = calculateMT(lep1_p4().pt(), lep1_p4().phi(), values_[met], values_[metphi]);
            values_[dphilmet] = deltaPhi(lep1_p4().phi(), values_[metphi]);
            values_[dphijmet] = (ngoodjets() > 1)? min2deltaPhi(values_[metphi], ak4pfjets_p4().at(0).phi(), ak4pfjets_p4().at(1).phi()) : -9;
          }

          suffix = "_metResDn";
        } else if (systype == 3) {
          values_[mt] = mt_met_lep_jup();
          values_[met] = pfmet_jup();
          values_[mlb] = Mlb_closestb_jup();
          values_[tmod] = topnessMod_jup();
          values_[njet] = jup_ngoodjets();
          values_[nbjet] = jup_ngoodbtags();
          values_[ntbtag] = jup_ntightbtags();
          values_[tfttag] = lead_tftopdisc_jup;
          values_[dphijmet] = mindphi_met_j1_j2_jup();
          values_[dphilmet] = deltaPhi(lep1_p4().phi(), pfmet_phi_jup());
          values_[j1passbtag] = (jup_ngoodjets() > 0)? jup_ak4pfjets_passMEDbtag().at(0) : 0;

          values_[jet1pt] = (jup_ngoodjets() > 0)? jup_ak4pfjets_p4().at(0).pt() : 0;
          values_[jet2pt] = (jup_ngoodjets() > 1)? jup_ak4pfjets_p4().at(1).pt() : 0;
          values_[jet1eta] = (jup_ngoodjets() > 0)? jup_ak4pfjets_p4().at(0).eta() : -9;
          values_[jet2eta] = (jup_ngoodjets() > 1)? jup_ak4pfjets_p4().at(1).eta() : -9;

          values_[metphi] = pfmet_phi_jup();
          values_[nbtag] = jup_nanalysisbtags();
          // values_[nsbtag] = (doTopTagging)? jup_nsoftbtags() - nsboverlep : 0;
          values_[leadbpt] = jup_ak4pfjets_leadbtag_p4().pt();
          values_[mlb_0b] = (jup_ak4pfjets_leadbtag_p4() + lep1_p4()).M();

          jestype_ = 1;
          suffix = "_jesUp";
        } else if (systype == 4) {
          values_[mt] = mt_met_lep_jdown();
          values_[met] = pfmet_jdown();
          values_[mlb] = Mlb_closestb_jdown();
          values_[tmod] = topnessMod_jdown();
          values_[njet] = jdown_ngoodjets();
          values_[nbjet] = jdown_ngoodbtags();
          values_[ntbtag] = jdown_ntightbtags();
          values_[tfttag] = lead_tftopdisc_jdown;
          values_[dphijmet] = mindphi_met_j1_j2_jdown();
          values_[dphilmet] = deltaPhi(lep1_p4().phi(), pfmet_phi_jdown());
          values_[j1passbtag] = (jdown_ngoodjets() > 0)? jdown_ak4pfjets_passMEDbtag().at(0) : 0;

          values_[jet1pt] = (jdown_ngoodjets() > 0)? jdown_ak4pfjets_p4().at(0).pt() : 0;
          values_[jet2pt] = (jdown_ngoodjets() > 1)? jdown_ak4pfjets_p4().at(1).pt() : 0;
          values_[jet1eta] = (jdown_ngoodjets() > 0)? jdown_ak4pfjets_p4().at(0).eta() : -9;
          values_[jet2eta] = (jdown_ngoodjets() > 1)? jdown_ak4pfjets_p4().at(1).eta() : -9;

          values_[metphi] = pfmet_phi_jdown();
          values_[nbtag] = jdown_nanalysisbtags();
          // values_[nsbtag] = (doTopTagging)? jdown_nsoftbtags() - nsboverlep : 0;
          values_[leadbpt] = jdown_ak4pfjets_leadbtag_p4().pt();
          values_[mlb_0b] = (jdown_ak4pfjets_leadbtag_p4() + lep1_p4()).M();

          jestype_ = 2;
          suffix = "_jesDn";
        }
        /// should do the same job as nanalysisbtags
        values_[nbtag] = (values_[mlb] > 175)? values_[ntbtag] : values_[nbjet];
        values_[njettmod] = (values_[njet] >= 4 || values_[tmod] >= 10);
        values_[passlmetcor] = (values_[lep1pt] < 50) || (values_[lep1pt] < (250 - 100*values_[dphilmet]));

        if (!gconf.is_data && runMETResCorrection) {
          values_[met_rs] = metobj.extras.met;
          values_[metphi_rs] = metobj.extras.phi;
          values_[mt_rs] = calculateMT(lep1_p4().pt(), lep1_p4().phi(), metobj.extras.met, metobj.extras.phi);
          values_[dphilmet_rs] = deltaPhi(lep1_p4().phi(), metobj.extras.phi);
          values_[dphijmet_rs] = (ngoodjets() > 1)? min2deltaPhi(metobj.extras.phi, ak4pfjets_p4().at(0).phi(), ak4pfjets_p4().at(1).phi()) : -9;
          if (applyMETResCorrection) {
            values_[met] = values_[met_rs];
            values_[metphi] = values_[metphi_rs];
            values_[mt] = values_[mt_rs];
            values_[dphilmet] = values_[dphilmet_rs];
            values_[dphijmet] = values_[dphijmet_rs];
          }
        } else {
          values_[met_rs] = values_[met];
          values_[metphi_rs] = values_[metphi];
          values_[mt_rs] = values_[mt];
          values_[dphilmet_rs] = values_[dphilmet];
          values_[dphijmet_rs] = values_[dphijmet];
        }

        // // Uncomment following lines if want to use CSV instead
        // values_[nbtag] = (values_[mlb] > 175)? ntbtagCSV : nbtagCSV;
        // values_[nbjet] = nbtagCSV;
        // values_[ntbtag] = ntbtagCSV;

        // Filling histograms for SR
        fillHistosForSR(suffix);

        fillHistosForCR0b(suffix);

        // Filling analysis variables with removed leptons, for CR2l
        if (systype == 0) {
          values_[mt_rl] = mt_met_lep_rl();
          values_[mt2_ll] = (doTopTagging)? MT2_ll() : 90;
          values_[met_rl] = pfmet_rl();
          values_[dphijmet_rl]= mindphi_met_j1_j2_rl();
          values_[dphilmet_rl] = lep1_dphiMET_rl();
          values_[tmod_rl] = topnessMod_rl();
        } else if (systype == 1 && values_[nlep_rl] > 1) {
          if (samplever.find("v31")) {
            values_[mt_rl] = mt_met_lep_rl_resup();
            values_[mt2_ll] = (doTopTagging)? MT2_ll_resup() : 90;
            values_[met_rl] = pfmet_rl_resup();
            values_[tmod_rl] = topnessMod_rl_resup();
            values_[dphijmet_rl] = mindphi_met_j1_j2_rl_resup();
            values_[dphilmet_rl] = lep1_dphiMET_rl_resup();
          } else if (runMETResCorrection) {
            LorentzVector newmet(values_[met]*cos(values_[metphi]), values_[met]*sin(values_[metphi]), 0, values_[met]);
            newmet += lep2_p4();
            values_[met_rl] = newmet.pt();
            values_[metphi_rl] = newmet.phi();
            values_[mt_rl] = calculateMT(lep1_p4().pt(), lep1_p4().phi(), values_[met_rl], values_[metphi_rl]);
            values_[dphilmet_rl] = deltaPhi(lep1_p4().phi(), values_[metphi_rl]);
            values_[dphijmet_rl] = (ngoodjets() > 1)? min2deltaPhi(values_[metphi_rl], ak4pfjets_p4().at(0).phi(), ak4pfjets_p4().at(1).phi()) : -9;
          }
        } else if (systype == 2 && values_[nlep_rl] > 1) {
          if (samplever.find("v31") == 0) {
            values_[mt_rl] = mt_met_lep_rl_resdown();
            values_[mt2_ll] = (doTopTagging)? MT2_ll_resdown() : 90;
            values_[met_rl] = pfmet_rl_resdown();
            values_[tmod_rl] = topnessMod_rl_resdown();
            values_[dphijmet_rl] = mindphi_met_j1_j2_rl_resdown();
            values_[dphilmet_rl] = lep1_dphiMET_rl_resdown();
          } else if (runMETResCorrection) {
            LorentzVector newmet(values_[met]*cos(values_[metphi]), values_[met]*sin(values_[metphi]), 0, values_[met]);
            newmet += lep2_p4();
            values_[met_rl] = newmet.pt();
            values_[metphi_rl] = newmet.phi();
            values_[mt_rl] = calculateMT(lep1_p4().pt(), lep1_p4().phi(), values_[met_rl], values_[metphi_rl]);
            values_[dphilmet_rl] = deltaPhi(lep1_p4().phi(), values_[metphi_rl]);
            values_[dphijmet_rl] = (ngoodjets() > 1)? min2deltaPhi(values_[metphi_rl], ak4pfjets_p4().at(0).phi(), ak4pfjets_p4().at(1).phi()) : -9;
          }
        } else if (jestype_ == 1) {
          values_[mt_rl] = mt_met_lep_rl_jup();
          values_[mt2_ll] = (doTopTagging)? MT2_ll_jup() : 90;
          values_[met_rl] = pfmet_rl_jup();
          values_[dphijmet_rl]= mindphi_met_j1_j2_rl_jup();
          values_[dphilmet_rl] = lep1_dphiMET_rl_jup();
          values_[tmod_rl] = topnessMod_rl_jup();
        } else if (jestype_ == 2) {
          values_[mt_rl] = mt_met_lep_rl_jdown();
          values_[mt2_ll] = (doTopTagging)? MT2_ll_jdown() : 90;
          values_[met_rl] = pfmet_rl_jdown();
          values_[dphijmet_rl]= mindphi_met_j1_j2_rl_jdown();
          values_[dphilmet_rl] = lep1_dphiMET_rl_jdown();
          values_[tmod_rl] = topnessMod_rl_jdown();
        }
        values_[passlmet_rl] = (values_[lep1pt] < 50) || (values_[lep1pt] < (250 - 100*values_[dphilmet_rl]));

        if (!gconf.is_data && runMETResCorrection && applyMETResCorrection && values_[nlep_rl] > 1) {
          LorentzVector newmet(metobj.extras.met*cos(metobj.extras.phi), metobj.extras.met*sin(metobj.extras.phi), 0, metobj.extras.met);
          newmet += lep2_p4();
          values_[met_rl] = newmet.pt();
          values_[metphi_rl] = newmet.phi();
          values_[mt_rl] = calculateMT(lep1_p4().pt(), lep1_p4().phi(), values_[met_rl], values_[metphi_rl]);
          values_[dphilmet_rl] = deltaPhi(lep1_p4().phi(), values_[metphi_rl]);
          values_[dphijmet_rl] = (ngoodjets() > 1)? min2deltaPhi(values_[metphi_rl], ak4pfjets_p4().at(0).phi(), ak4pfjets_p4().at(1).phi()) : -9;
        }

        fillHistosForCRemu(suffix);

        fillHistosForCR2l(suffix);

        if (systype == 0)
          testCutFlowHistos(testVec[2]);
        // fillTopTaggingHistos(suffix);

        // Also do yield using genmet for fastsim samples <-- under development
        if (is_fastsim_ && systype == 0 && babyver >= 31.4) {
          values_[met] = genmet();
          values_[mt] = mt_genmet_lep();
          values_[tmod] = topnessMod_genmet();
          values_[dphijmet] = mindphi_genmet_j1_j2();
          values_[dphilmet] = lep1_dphiGenMET();
          values_[passlmetcor] = (values_[lep1pt] < 50) || (values_[lep1pt] < (250 - 100*values_[dphilmet]));

          values_[met_rl] = genmet();
          values_[mt2_ll] = MT2_ll_genmet();
          values_[mt_rl] = mt_genmet_lep();
          values_[tmod_rl] = topnessMod_genmet();
          values_[dphijmet_rl] = mindphi_genmet_j1_j2();
          values_[dphilmet_rl] = lep1_dphiGenMET();
          values_[passlmet_rl] = (values_[lep1pt] < 50) || (values_[lep1pt] < (250 - 100*values_[dphilmet_rl]));

          fillHistosForSR("_genmet");
          fillHistosForCR0b("_genmet");
          fillHistosForCR2l("_genmet");
          // fillHistosForSR("");
          // fillHistosForCR0b("");
          // fillHistosForCR2l("");
        }

      }  // end of jes variation

      */
    } // end of event loop


    delete tree;
    file.Close();
  } // end of file loop

  cout << "[StopLooper::looper] processed  " << nEventsTotal << " events" << endl;
  if ( nEventsChain != nEventsTotal )
    cout << "WARNING: Number of events from files is not equal to total number of events" << endl;

  outfile_->cd();

  auto writeHistsToFile = [&] (vector<SR>& srvec) {
    for (auto& sr : srvec) {
      TDirectory* dir = (TDirectory*) outfile_->Get(sr.GetName().c_str());
      if (dir == 0) dir = outfile_->mkdir(sr.GetName().c_str()); // shouldn't happen
      dir->cd();
      for (auto& h : sr.histMap) {
        if (h.first.find("HI") != string::npos || h.first.find("LOW") != string::npos) continue;
        // Move overflows of the yield hist to the last bin of histograms
        if (h.first.find("h_metbins") != string::npos) {
          moveOverFlowToLastBin1D(h.second);
          zeroOutNegativeYields(h.second);
        } else if (h.first.find("hSMS_metbins") != string::npos) {
          moveXOverFlowToLastBin3D((TH3*) h.second);
        }
        h.second->Write();
      }
    }
  };

  writeHistsToFile(testVec);
  writeHistsToFile(SRVec);
  writeHistsToFile(CR0bVec);
  writeHistsToFile(CR2lVec);
  writeHistsToFile(CRemuVec);

  auto writeRatioHists = [&] (const SR& sr) {
    for (const auto& h : sr.histMap) {
      if (h.first.find("hnum") != 0) continue;
      string hname = h.first;
      hname.erase(0, 4);
      // dummy.cd();
      TH1F* h_ratio = (TH1F*) h.second->Clone(("ratio"+hname).c_str());
      h_ratio->SetDirectory(0);
      h_ratio->Divide(h_ratio, sr.histMap.at("hden"+hname), 1, 1, "B");

      TDirectory * dir = (TDirectory*) outfile_->Get(sr.GetName().c_str());
      if (dir == 0) dir = outfile_->mkdir(sr.GetName().c_str()); // shouldn't happen
      dir->cd();
      h_ratio->Write();
    }
  };
  writeRatioHists(testVec[0]);
  writeRatioHists(testVec[1]);

  outfile_->Write();
  outfile_->Close();
  if (printPassedEvents) ofile.close();

  bmark->Stop("benchmark");
  cout << endl;
  cout << nEventsTotal << " Events Processed, where " << nDuplicates << " duplicates were skipped, and ";
  cout << nPassedTotal << " Events passed all selections." << endl;
  cout << "------------------------------" << endl;
  cout << "CPU  Time:   " << Form( "%.01f s", bmark->GetCpuTime("benchmark")  ) << endl;
  cout << "Real Time:   " << Form( "%.01f s", bmark->GetRealTime("benchmark") ) << endl;
  cout << endl;
  delete bmark;

  return;
}

void StopLooper::fillYieldHistos(SR& sr, float met, string suf, bool is_cr2l) {

  string srname = sr.GetName();
  int cortype = 0;
  if (srname.back() == '2') cortype = 2;
  if (srname.back() == '3') cortype = 3;
  if (srname.find("I") != string::npos) cortype = 4;
  if (srname.find("J") != string::npos) cortype = 5;

  // FIXME: evtweight_ = evtWgt.getWeight(evtWgtInfo::systID(jestype_), is_cr2l, cortype);
  evtweight_ = 1;

  // if (doNvtxReweight) {
  //   if (PV_npvs() < 100) evtweight_ *= nvtxscale_[PV_npvs()];  // only scale for data
  // }

  if (evtweight_ == 0) cout << "[StopLooper::fillYieldHistos]: WARNING: the event weight is 0 at evt = " << event() << "!!\n";
  if (isnan(evtweight_)) {
    cout << "[StopLooper::fillYieldHistos]: WARNING: the event weight is NAN!! Skipping!" << ", SR=" << srname << endl;
    return;
  }

  auto fillhists = [&](string s, int filldim = 1) {
    if (filldim == 3)
      plot3d("hSMS_metbins"+s+suf, met, mstop_, mlsp_, evtweight_, sr.histMap, ";E^{miss}_{T} [GeV];M_{stop};M_{LSP}",
             sr.GetNMETBins(), sr.GetMETBinsPtr(), mStopBins.size()-1, mStopBins.data(), mLSPBins.size()-1, mLSPBins.data());
    else
      plot1d("h_metbins"+s+suf, met, evtweight_, sr.histMap, (sr.GetName()+":"+sr.GetDetailName()+";E^{miss}_{T} [GeV]").c_str(),
             sr.GetNMETBins(), sr.GetMETBinsPtr());

    // TODO: systematic variation 
    // if (doSystVariations && (is_bkg_ || (is_fastsim_ && filldim == 3)) && suf == "") {
    //   // Only run once when filling for nominal. JES or METRes variation dealt with above. No need for signals?
    //   for (int isyst = 3; isyst < evtWgtInfo::k_nSyst; ++isyst) {
    //     auto syst = (evtWgtInfo::systID) isyst;
    //     if (evtWgt.doingSystematic(syst)) {
    //       float systwgt = evtWgt.getWeight(syst, is_cr2l, cortype);
    //       if (evtweight_ != 0 && systwgt == 0 && syst != evtWgtInfo::k_L1prefireUp)
    //         cout << "[StopLooper::fillYieldHistos]: WARNING: the event weight is 0 for syst = " << evtWgt.getLabel(syst) << " at evt = " << event() << "!!" << ", SR=" << srname << endl;
    //       if (isnan(systwgt))
    //         cout << "[StopLooper::fillYieldHistos]: WARNING: the event weight is NAN for syst = " << evtWgt.getLabel(syst) << " at evt = " << event() << "!!" << ", SR=" << srname << endl;
    //       if (filldim == 3)
    //         plot3d("hSMS_metbins"+s+"_"+evtWgt.getLabel(syst), met, mstop_, mlsp_, systwgt, sr.histMap, ";E^{miss}_{T} [GeV];M_{stop};M_{LSP}",
    //                sr.GetNMETBins(), sr.GetMETBinsPtr(), mStopBins.size()-1, mStopBins.data(), mLSPBins.size()-1, mLSPBins.data());
    //       else
    //         plot1d("h_metbins"+s+"_"+evtWgt.getLabel(syst), met, systwgt, sr.histMap,
    //                (sr.GetName()+":"+evtWgt.getLabel(syst)+";E^{miss}_{T} [GeV]").c_str(), sr.GetNMETBins(), sr.GetMETBinsPtr());
    //       if (fillWeights && s == "") {
    //         plot1d("h_weights"+s+"_"+evtWgt.getLabel(syst), systwgt/evtweight_, 1, sr.histMap,
    //                (sr.GetName()+":"+evtWgt.getLabel(syst)+";event weights").c_str(), 200, 0, 2);
    //       }
    //     }
    //   }
    // }
  };
  fillhists("", (is_fastsim_)? 3 : 1);
  if (!is_fastsim_ && mlsp_ > 0 && mstop_ > 0) fillhists("", 3);  // for fullsim signal

  // TODO: gen category
  // if (doGenClassification && is_bkg_) {
  //   // Only fill gen classification for background events, used for background estimation
  //   if (isZtoNuNu()) fillhists("_Znunu");
  //   else if (is2lep()) fillhists("_2lep");
  //   else if (is1lepFromW()) fillhists("_1lepW");
  //   else if (is1lepFromTop()) fillhists("_1lepTop");
  //   else fillhists("_unclass");  // either unclassified 1lep or 0lep, or something else unknown, shouldn't have (m)any
  // }

  // Block for debugging, active when setting printPassedEvents = true
  if (printPassedEvents && met > 750 && suf == "" && sr.GetName().find("sr") == 0) {
    ofile << run() << ' ' << luminosityBlock() << ' ' << event() << ' ' << sr.GetName();
    ofile << ": values_[met]= " << values_[Vars::met] << ", values_[mt]= " << values_[mt] << ", values_[njet]= " << values_[njet] << ", values_[nbjet] = " << values_[nbjet]  << ", values_[tmod]= " << values_[tmod] << ", values_[mlb]= " << values_[mlb];
    ofile << endl;
  }
}

void StopLooper::fillHistosForSR(string suf) {

  /*
  // Trigger requirements
  if (gconf.is_data && !PassingHLTriggers()) return;

  // // For getting into full trigger efficiency in 2017 data
  // if ( (abs(lep1_pdgid()) == 11 && values_[lep1pt] < 40) || (abs(lep1_pdgid()) == 13 && values_[lep1pt] < 30) ) return;

  for (auto& sr : SRVec) {
    if (!sr.PassesSelection(values_)) continue;
    fillYieldHistos(sr, values_[met], suf);


    if (runYieldsOnly) continue;
    // Skip histogram ploting for individul signal regions
    if (sr.GetName().at(2) >= 'A' && sr.GetName().at(2) <= 'H') continue;

    // Plot kinematics histograms
    auto fillKineHists = [&](string s) {
      // Simple plot function plot1d to add extra plots anywhere in the code, is great for quick checks
      plot1d("h_mt"+s,       values_[mt]      , evtweight_, sr.histMap, ";M_{T} [GeV]"                   , 20,  150, 650);
      plot1d("h_met"+s,      values_[met]     , evtweight_, sr.histMap, ";E_{T}^{miss} [GeV]"            , 24,  250, 850);
      plot1d("h_metphi"+s,   values_[metphi]  , evtweight_, sr.histMap, ";#phi(E_{T}^{miss})"            , 34, -3.4, 3.4);
      plot1d("h_lep1pt"+s,   values_[lep1pt]  , evtweight_, sr.histMap, ";p_{T}(lepton) [GeV]"           , 24,    0, 600);
      plot1d("h_lep1eta"+s,  values_[lep1eta] , evtweight_, sr.histMap, ";#eta(lepton)"                  , 30,   -3,   3);
      plot1d("h_nleps"+s,    values_[nlep]    , evtweight_, sr.histMap, ";Number of leptons"             ,  5,    0,   5);
      plot1d("h_njets"+s,    values_[njet]    , evtweight_, sr.histMap, ";Number of jets"                ,  8,    2,  10);
      plot1d("h_nbjets"+s,   values_[nbjet]   , evtweight_, sr.histMap, ";Number of b-tagged jets"       ,  4,    1,   5);
      plot1d("h_ntbtags"+s,  values_[ntbtag]  , evtweight_, sr.histMap, ";Number of tight b-tagged jets" ,  5,    0,   5);
      plot1d("h_mlepb"+s,    values_[mlb]     , evtweight_, sr.histMap, ";M_{#it{l}b} [GeV]"             , 24,    0, 600);
      plot1d("h_dphijmet"+s, values_[dphijmet], evtweight_, sr.histMap, ";#Delta#phi(jet,E_{T}^{miss})"  , 25,  0.8, 3.3);
      plot1d("h_tmod"+s,     values_[tmod]    , evtweight_, sr.histMap, ";Modified topness"              , 25,  -10,  15);
      plot1d("h_nvtxs"+s,    values_[nvtx]    , evtweight_, sr.histMap, ";Number of vertices"            , 100,   1, 101);

      if (doTopTagging) {
        plot1d("h_nak8jets"+s, values_[nak8jets], evtweight_, sr.histMap, ";Number of AK8 jets"   , 7, 0, 7);
        plot1d("h_resttag"+s,  values_[resttag] , evtweight_, sr.histMap, ";resolved top tag"     , 110, -1.1f, 1.1f);
        plot1d("h_bdtttag"+s,  values_[bdtttag] , evtweight_, sr.histMap, ";BDT resolved top tag" , 110, -1.1f, 1.1f);
        plot1d("h_tfttag"+s,   values_[tfttag]  , evtweight_, sr.histMap, ";Leading DeepResolved top tag", 120, -0.1f, 1.1f);
        plot1d("h_deepttag"+s, values_[deepttag], evtweight_, sr.histMap, ";Leading DeepAK8 top tag"     , 120, -0.1f, 1.1f);
        plot1d("h_nsbtags"+s,  values_[nsbtag]  , evtweight_, sr.histMap, ";Number of soft b-tagged jets",  5,  0, 5);
        plot1d("h_nsblep"+s,   values_[nsblep]  , evtweight_, sr.histMap, ";Number of soft b overlap lep"  ,  5,  0, 5);
      }

      // if ( (HLT_SingleEl() && abs(lep1_pdgid()) == 11 && values_[lep1pt] < 45) || (HLT_SingleMu() && abs(lep1_pdgid()) == 13 && values_[lep1pt] < 40) || (HLT_MET_MHT() && pfmet() > 250) ) {
      if (true) {
        plot1d("h_mt_h"+s,       values_[mt]       , evtweight_, sr.histMap, ";M_{T} [GeV]"                   , 24,   0, 600);
        plot1d("h_met_h"+s,      values_[met]      , evtweight_, sr.histMap, ";E_{T}^{miss} [GeV]"            , 32,  50, 850);
        plot1d("h_nbtags"+s,     values_[nbjet]    , evtweight_, sr.histMap, ";Number of b-tagged jets"       ,  5,   0,   5);
        plot1d("h_dphijmet_h"+s, values_[dphijmet] , evtweight_, sr.histMap, ";#Delta#phi(jet,E_{T}^{miss})"  , 33, 0.0, 3.3);
        plot1d("h_dphilmet"+s,   values_[dphilmet] , evtweight_, sr.histMap, ";#Delta#phi(l,MET)"             , 32,   0, 3.2);
      }

      plot1d("h_jet1pt"+s,  values_[jet1pt],  evtweight_, sr.histMap, ";p_{T}(jet1) [GeV]"  , 32,  0, 800);
      plot1d("h_jet2pt"+s,  values_[jet2pt],  evtweight_, sr.histMap, ";p_{T}(jet2) [GeV]"  , 32,  0, 800);
      plot1d("h_jet1eta"+s, values_[jet1eta], evtweight_, sr.histMap, ";#eta(jet1) [GeV]"   , 30,  -3,  3);
      plot1d("h_jet2eta"+s, values_[jet2eta], evtweight_, sr.histMap, ";#eta(jet2) [GeV]"   , 60,  -3,  3);

      if (suf == "" && babyver >= 31.2 && ngoodjets() > 0) {
        plot2d("h2d_ht5_dphijmht", deltaPhi(ak4pfjets_p4().at(0).phi(), values_[mhtphi]), ak4_HTeta5()/ak4_HT(), evtweight_, sr.histMap, ";#Delta#phi(j1,MHT);HT5/HT", 64, 0, 3.2, 80, 1, 3);
        plot1d("h_ht5overht", ak4_HTeta5()/ak4_HT(), evtweight_, sr.histMap, ";HT5/HT", 80, 1, 3);
      }
      for (auto sb : softtags_p4()) {
        plot1d("h_softbs_pt"+s, sb.pt() , evtweight_, sr.histMap, ";p_{T} (soft b)",  20,  0, 20);
      }
    };

    if (suf == "") fillKineHists(suf);
    if (is_fastsim_ && suf == "" && (checkMassPt(1050, 100) || checkMassPt(900, 500) || checkMassPt(950, 100) || checkMassPt(425, 325) ||
                                     checkMassPt(1000, 100) || checkMassPt(800, 450) || checkMassPt(750, 400) ||
                                     checkMassPt(1200, 100) || checkMassPt(850, 100) || checkMassPt(650, 350)))
      fillKineHists(Form("%s_%d_%d%s", sigtype.c_str(), mstop_, mlsp_, suf.c_str()));

    // Separate contribution by gen classification for background events
    if (doGenClassification && is_bkg_ && suf == "") {
      if (isZtoNuNu()) fillKineHists("_Znunu");
      else if (is2lep()) fillKineHists("_2lep");
      else if (is1lepFromW()) fillKineHists("_1lepW");
      else if (is1lepFromTop()) fillKineHists("_1lepTop");
      else fillKineHists("_unclass");  // either unclassified 1lep or 0lep, or something else unknown, shouldn't have (m)any
    }

  }
  */

}

void StopLooper::fillHistosForCR2l(string suf) {

  /*
  // Trigger requirements
  if (gconf.is_data && !PassingHLTriggers(2)) return;

  // For getting into full trigger efficiency in 2017 & 2018 data
  // if (not( (HLT_SingleEl() && abs(lep1_pdgid()) == 11 && values_[lep1pt] < 45) ||
  //          (HLT_SingleMu() && abs(lep1_pdgid()) == 13 && values_[lep1pt] < 40) || (HLT_MET_MHT() && pfmet() > 250) )) return;

  for (auto& cr : CR2lVec) {
    if (!cr.PassesSelection(values_)) continue;
    fillYieldHistos(cr, values_[met_rl], suf, true);

    if (runYieldsOnly) continue;
    if (cr.GetName().at(4) >= 'A' && cr.GetName().at(4) <= 'H') continue;

    auto fillKineHists = [&] (string s) {
      plot1d("h_finemet"+s,    values_[met]         , evtweight_, cr.histMap, ";E_{T}^{miss} [GeV]"            , 80,  0, 800);
      plot1d("h_met"+s,        values_[met]         , evtweight_, cr.histMap, ";E_{T}^{miss} [GeV]"            , 24, 50, 850);
      plot1d("h_metphi"+s,     values_[metphi]      , evtweight_, cr.histMap, ";#phi(E_{T}^{miss})"            , 40,  -4, 4);
      plot1d("h_mt"+s,         values_[mt]          , evtweight_, cr.histMap, ";M_{T} [GeV]"                   , 24,  0, 600);
      plot1d("h_tmod"+s,       values_[tmod]        , evtweight_, cr.histMap, ";Modified topness"              , 25, -10, 15);
      plot1d("h_njets"+s,      values_[njet]        , evtweight_, cr.histMap, ";Number of jets"                ,  8,  2, 10);
      plot1d("h_nbjets"+s,     values_[nbjet]       , evtweight_, cr.histMap, ";Number of b-tagged jets"       ,  4,  1, 5);
      plot1d("h_ntbtags"+s,    values_[ntbtag]      , evtweight_, cr.histMap, ";Number of tight b-tagged jets" ,  5,  0, 5);
      plot1d("h_nleps"+s,      values_[nlep_rl]     , evtweight_, cr.histMap, ";nleps (dilep)"                 ,  5,  0, 5);
      plot1d("h_lep1pt"+s,     values_[lep1pt]      , evtweight_, cr.histMap, ";p_{T}(lepton) [GeV]"           , 24,  0, 600);
      plot1d("h_lep1eta"+s,    values_[lep1eta]     , evtweight_, cr.histMap, ";#eta(lepton)"                  , 30, -3, 3);
      plot1d("h_mlepb"+s,      values_[mlb]         , evtweight_, cr.histMap, ";M_{#it{l}b} [GeV]"             , 24,  0, 600);
      plot1d("h_dphijmet"+s,   values_[dphijmet]    , evtweight_, cr.histMap, ";#Delta#phi(jet,E_{T}^{miss})"  , 33,  0, 3.3);
      plot1d("h_nvtxs"+s,      values_[nvtx]        , evtweight_, cr.histMap, ";Number of vertices"            , 70,  1, 71);
      plot1d("h_rlmet"+s,      values_[met_rl]      , evtweight_, cr.histMap, ";E_{T}^{miss} (removed lepton) [GeV]" , 24, 250, 850);
      plot1d("h_rlmt"+s,       values_[mt_rl]       , evtweight_, cr.histMap, ";M_{T} (removed lepton) [GeV]"  , 20,  150, 650);
      plot1d("h_rltmod"+s,     values_[tmod_rl]     , evtweight_, cr.histMap, ";Modified topness"              , 25, -10, 15);
      plot1d("h_rldphijmet"+s, values_[dphijmet_rl] , evtweight_, cr.histMap, ";#Delta#phi(jet,(E_{T}^{miss}+l_{2}))" , 33, 0.8, 3.3);
      plot1d("h_rldphilmet"+s, values_[dphilmet_rl] , evtweight_, cr.histMap, ";#Delta#phi(l,MET(rl))"         , 32,   0, 3.2);

      if (doTopTagging) {
        plot1d("h_nak8jets"+s, values_[nak8jets] , evtweight_, cr.histMap, ";Number of AK8 jets"   , 7, 0, 7);
        plot1d("h_resttag"+s,  values_[resttag]  , evtweight_, cr.histMap, ";resolved top tag"     , 110, -1.1f, 1.1f);
        plot1d("h_bdtttag"+s,  values_[bdtttag]  , evtweight_, cr.histMap, ";BDT resolved top tag" , 110, -1.1f, 1.1f);
        plot1d("h_tfttag"+s,   values_[tfttag]   , evtweight_, cr.histMap, ";Leading DeepResolved top tag", 120, -0.1f, 1.1f);
        plot1d("h_deepttag"+s, values_[deepttag] , evtweight_, cr.histMap, ";Leading DeepAK8 top tag"     , 120, -0.1f, 1.1f);
        plot1d("h_nsbtags"+s,  values_[nsbtag]   , evtweight_, cr.histMap, ";Number of soft b-tagged jets"  ,  5,  0, 5);
        plot1d("h_nsblep"+s,   values_[nsblep]   , evtweight_, cr.histMap, ";Number of soft b overlap lep"  ,  5,  0, 5);
      }

      // if ( (HLT_SingleEl() && abs(lep1_pdgid()) == 11 && values_[lep1pt] < 45) || (HLT_SingleMu() && abs(lep1_pdgid()) == 13 && values_[lep1pt] < 40) || (HLT_MET_MHT() && pfmet() > 250) ) {
      if (true) {
        plot1d("h_rlmet_h"+s,  values_[met_rl]   , evtweight_, cr.histMap, ";E_{T}^{miss} (removed lepton) [GeV]"  , 32, 50, 850);
        plot1d("h_rlmt_h"+s,   values_[mt_rl]    , evtweight_, cr.histMap, ";M_{T} (removed lepton) [GeV]"  , 26, 0, 650);
        plot1d("h_nbtags"+s,   values_[nbjet]    , evtweight_, cr.histMap, ";Number of b-tagged jets"       ,  5,  0, 5);
        plot1d("h_dphilmet"+s, values_[dphilmet] , evtweight_, cr.histMap, ";#Delta#phi(l,MET)"             , 32,   0, 3.2);
      }

      plot1d("h_jet1pt"+s,  values_[jet1pt],  evtweight_, cr.histMap, ";p_{T}(jet1) [GeV]"  , 32,  0, 800);
      plot1d("h_jet2pt"+s,  values_[jet2pt],  evtweight_, cr.histMap, ";p_{T}(jet2) [GeV]"  , 32,  0, 800);
      plot1d("h_jet1eta"+s, values_[jet1eta], evtweight_, cr.histMap, ";#eta(jet1) [GeV]"   , 30,  -3,  3);
      plot1d("h_jet2eta"+s, values_[jet2eta], evtweight_, cr.histMap, ";#eta(jet2) [GeV]"   , 60,  -3,  3);
      plot1d("h_lep2pt"+s,   lep2_p4().pt()  , evtweight_, cr.histMap, ";p_{T}(lep2) [GeV]" , 24,     0,  600);
      plot1d("h_lep2eta"+s,  lep2_p4().eta() , evtweight_, cr.histMap, ";#eta(lep2)"        , 50,  -2.5,  2.5);

      plot1d("h_mht"+s,           values_[mht]       , evtweight_, cr.histMap, ";H_{T}^{miss} [GeV]"                , 40, 0, 800);
      plot1d("h_mhtovermet"+s, values_[mht]/values_[met], evtweight_, cr.histMap, ";H_{T}^{miss}/E_{T}^{miss}{miss}"  , 40, 0, 4);
      plot1d("h_diffmhtmet"+s, fabs(values_[mht]-values_[met])/values_[met], evtweight_, cr.histMap, ";#Delta(H_{T}^{miss}, E_{T}^{miss}{miss})/E_{T}^{miss}{miss}"   , 40, 0, 1);
      // Luminosity test at Z peak
      if (suf == "" && babyver >= 31.2 && ngoodjets() > 0) {
        plot2d("h2d_ht5_dphijmht", deltaPhi(ak4pfjets_p4().at(0).phi(), values_[mhtphi]), ak4_HTeta5()/ak4_HT(), evtweight_, cr.histMap, ";#Delta#phi(j1,MHT);HT5/HT", 64, 0, 3.2, 80, 1, 3);
        plot1d("h_ht5overht", ak4_HTeta5()/ak4_HT(), evtweight_, cr.histMap, ";HT5/HT", 80, 1, 3);
      }
    };
    if (suf == "") fillKineHists(suf);
    if (is_fastsim_ && suf == "" && (checkMassPt(1200, 50) || checkMassPt(800, 400) || checkMassPt(800, 675)))
      fillKineHists(Form("_%d_%d%s", mstop_, mlsp_, suf.c_str()));

    // Separate contribution by gen classification for background events
    if (doGenClassification && is_bkg_ && suf == "") {
      if (isZtoNuNu()) fillKineHists("_Znunu");
      else if (is2lep()) fillKineHists("_2lep");
      else if (is1lepFromW()) fillKineHists("_1lepW");
      else if (is1lepFromTop()) fillKineHists("_1lepTop");
      else fillKineHists("_unclass");  // either unclassified 1lep or 0lep, or something else unknown, shouldn't have (m)any
    }

    auto fillExtraZHists = [&] (string s) {
      plot1d("h_mll"+s,   values_[mll], evtweight_, cr.histMap, ";M_{#it{ll}} [GeV]" , 120, 0, 240 );
      if (82 < values_[mll] && values_[mll] < 100) {
        plot1d("h_zpt"+s, (lep1_p4() + lep2_p4()).pt(), evtweight_, cr.histMap, ";p_{T}(Z) [GeV]"          , 200, 0, 200);
        plot1d("h_njets_onZ"+s,   values_[njet]  , evtweight_, cr.histMap, ";Number of jets"          , 12,  0, 12);
        plot1d("h_nbjets_onZ"+s,  values_[nbjet] , evtweight_, cr.histMap, ";Number of b-tagged jets" ,  6,  0, 6);
      } else {
        plot1d("h_njets_offZ"+s,  values_[njet]  , evtweight_, cr.histMap, ";Number of jets"          , 12,  0, 12);
        plot1d("h_nbjets_offZ"+s, values_[nbjet] , evtweight_, cr.histMap, ";Number of b-tagged jets" ,  6,  0, 6);
      }
    };
    if (suf == "" && lep1_pdgid() == -lep2_pdgid())
      fillExtraZHists(suf);

    auto fillExtraLep2Hists = [&] (string s) {
      plot1d("h_lep2pt"+s,   lep2_p4().pt()  , evtweight_, cr.histMap, ";p_{T}(lep2) [GeV]" , 24,     0,  600);
      plot1d("h_lep2eta"+s,  lep2_p4().eta() , evtweight_, cr.histMap, ";#eta(lep2)"        , 50,  -2.5,  2.5);
      plot1d("h_lep2phi"+s,  lep2_p4().phi() , evtweight_, cr.histMap, ";#phi(lep2)"        , 63, -3.15, 3.15);
      // HEM happening at -4.7 < eta < -1.4, -1.6 < phi < -0.8
      plot2d("h2d_lep2phi_eta"+s, lep2_p4().eta(), lep2_p4().phi(), evtweight_, cr.histMap, ";#eta(lep2);#phi(lep2)" , 24, -2.4, 2.4, 32, -3.2, 3.2);
    };
    if (suf == "" && year_ == 2018 && abs(lep2_pdgid()) == 11) {
      fillExtraLep2Hists("_el");
      if (gconf.is_data &&  run() < 319077)
        fillExtraLep2Hists("_el_preHEM");
      else if (gconf.is_data && run() >= 319077)
        fillExtraLep2Hists("_el_posHEM");
    }

    // if ( HLT_SingleMu() ) fillKineHists(suf+"_hltmu");
    // if ( HLT_SingleEl() ) fillKineHists(suf+"_hltel");
    // if ( HLT_MET_MHT() )  fillKineHists(suf+"_hltmet");

    // if ( abs(lep1_pdgid()*lep2_pdgid()) == 121 )
    //   fillKineHists(suf+"_ee");
    // else if ( abs(lep1_pdgid()*lep2_pdgid()) == 143 )
    //   fillKineHists(suf+"_emu");
    // else if ( abs(lep1_pdgid()*lep2_pdgid()) == 169 )
    //   fillKineHists(suf+"_mumu");

    } */
}

void StopLooper::fillHistosForCR0b(string suf) {

  /*
  // Trigger requirements
  if (gconf.is_data && !PassingHLTriggers()) return;

  for (auto& cr : CR0bVec) {
    if (!cr.PassesSelection(values_)) continue;
    fillYieldHistos(cr, values_[met], suf);

    if (runYieldsOnly) continue;
    if (cr.GetName().at(4) >= 'A' && cr.GetName().at(4) <= 'H') continue;

    auto fillKineHists = [&] (string s) {
      plot1d("h_mt"+s,       values_[mt]      , evtweight_, cr.histMap, ";M_{T} [GeV]"          , 20, 150, 650);
      plot1d("h_met"+s,      values_[met]     , evtweight_, cr.histMap, ";E_{T}^{miss} [GeV]"   , 24, 250, 850);
      plot1d("h_metphi"+s,   values_[metphi]  , evtweight_, cr.histMap, ";#phi(E_{T}^{miss})"   , 34, -3.4, 3.4);
      plot1d("h_lep1pt"+s,   values_[lep1pt]  , evtweight_, cr.histMap, ";p_{T}(lepton) [GeV]"  , 24,  0, 600);
      plot1d("h_lep1eta"+s,  values_[lep1eta] , evtweight_, cr.histMap, ";#eta(lepton)"         , 30, -3, 3);
      plot1d("h_nleps"+s,    values_[nlep]    , evtweight_, cr.histMap, ";Number of leptons"    ,  5,  0, 5);
      plot1d("h_njets"+s,    values_[njet]    , evtweight_, cr.histMap, ";Number of jets"       ,  8,  2, 10);
      plot1d("h_nbjets"+s,   values_[nbjet]   , evtweight_, cr.histMap, ";Number of b-tagged jets" ,  5,  0, 5);
      plot1d("h_mlepb"+s,    values_[mlb_0b]  , evtweight_, cr.histMap, ";M_{#it{l}b} [GeV]"  , 24,  0, 600);
      plot1d("h_dphijmet"+s, values_[dphijmet], evtweight_, cr.histMap, ";#Delta#phi(jet,E_{T}^{miss})"  , 33,  0, 3.3);
      plot1d("h_tmod"+s,     values_[tmod]    , evtweight_, cr.histMap, ";Modified topness"     , 25, -10, 15);
      plot1d("h_nvtxs"+s,    values_[nvtx]    , evtweight_, cr.histMap, ";Number of vertices"   , 70,  1, 71);
      plot1d("h_ntbtags"+s,  values_[ntbtag]  , evtweight_, cr.histMap, ";Number of tight b-tagged jets",  5,  0, 5);

      if (doTopTagging) {
        plot1d("h_nak8jets"+s, values_[nak8jets], evtweight_, cr.histMap, ";Number of AK8 jets"   , 7, 0, 7);
        plot1d("h_resttag"+s,  values_[resttag] , evtweight_, cr.histMap, ";resolved top tag"     , 110, -1.1f, 1.1f);
        plot1d("h_bdtttag"+s,  values_[bdtttag] , evtweight_, cr.histMap, ";BDT resolved top tag" , 110, -1.1f, 1.1f);
        plot1d("h_tfttag"+s,   values_[tfttag]  , evtweight_, cr.histMap, ";Leading DeepResolved top tag", 120, -0.1f, 1.1f);
        plot1d("h_deepttag"+s, values_[deepttag], evtweight_, cr.histMap, ";Leading DeepAK8 top tag"     , 120, -0.1f, 1.1f);
        plot1d("h_nsbtags"+s,  values_[nsbtag]  , evtweight_, cr.histMap, ";Number of soft b-tagged jets",  5,  0, 5);
      }

      if (true) {
        plot1d("h_mt_h"+s,       values_[mt]       , evtweight_, cr.histMap, ";M_{T} [GeV]"                   , 24,   0, 600);
        plot1d("h_met_h"+s,      values_[met]      , evtweight_, cr.histMap, ";E_{T}^{miss} [GeV]"            , 32,  50, 850);
        plot1d("h_nbtags"+s,     values_[nbjet]    , evtweight_, cr.histMap, ";Number of b-tagged jets"       ,  5,   0,   5);
        plot1d("h_dphijmet_h"+s, values_[dphijmet] , evtweight_, cr.histMap, ";#Delta#phi(jet,E_{T}^{miss})"  , 33, 0.0, 3.3);
        plot1d("h_dphilmet"+s,   values_[dphilmet] , evtweight_, cr.histMap, ";#Delta#phi(l,MET)"             , 32,   0, 3.2);
      }

      if (runMETResCorrection) {
        plot1d("h_met_rs"+s,      values_[met_rs]      , evtweight_, cr.histMap, ";E_{T}^{miss} [GeV]"            , 32,  50, 850);
        plot1d("h_metphi_rs"+s,   values_[metphi_rs]   , evtweight_, cr.histMap, ";#phi(E_{T}^{miss})"            , 34, -3.4, 3.4);
        plot1d("h_mt_rs"+s,       values_[mt_rs]       , evtweight_, cr.histMap, ";M_{T} [GeV]"                   , 40, 0, 400);
        plot1d("h_dphijmet_rs"+s, values_[dphijmet_rs] , evtweight_, cr.histMap, ";#Delta#phi(jet,E_{T}^{miss})"  , 25,  0.8, 3.3);
        plot1d("h_dphilmet_rs"+s, values_[dphilmet_rs] , evtweight_, cr.histMap, ";#Delta#phi(l,MET)"             , 32,   0, 3.2);
      }

      plot1d("h_jet1pt"+s,  values_[jet1pt],  evtweight_, cr.histMap, ";p_{T}(jet1) [GeV]"  , 32,  0, 800);
      plot1d("h_jet2pt"+s,  values_[jet2pt],  evtweight_, cr.histMap, ";p_{T}(jet2) [GeV]"  , 32,  0, 800);
      plot1d("h_jet1eta"+s, values_[jet1eta], evtweight_, cr.histMap, ";#eta(jet1) [GeV]"   , 30,  -3,  3);
      plot1d("h_jet2eta"+s, values_[jet2eta], evtweight_, cr.histMap, ";#eta(jet2) [GeV]"   , 60,  -3,  3);
      // Temporary test for low dphijmet excess
      if (jestype_ == 0) {
        plot1d("h_dphij1j2"+s, fabs(ak4pfjets_p4().at(0).phi()-ak4pfjets_p4().at(1).phi()), evtweight_, cr.histMap, ";#Delta#phi(j1,j2)" , 33,  0, 3.3);
        plot1d("h_mht"+s,           values_[mht]       , evtweight_, cr.histMap, ";H_{T}^{miss} [GeV]"                , 40, 0, 800);
        plot1d("h_mhtovermet"+s, values_[mht]/values_[met], evtweight_, cr.histMap, ";H_{T}^{miss}/E_{T}^{miss}{miss}"  , 40, 0, 4);
        plot1d("h_diffmhtmet"+s, fabs(values_[mht]-values_[met])/values_[met], evtweight_, cr.histMap, ";#Delta(H_{T}^{miss}, E_{T}^{miss}{miss})/E_{T}^{miss}{miss}"   , 40, 0, 1);
      }
      if (suf == "" && babyver >= 31.2 && ngoodjets() > 0) {
        plot2d("h2d_ht5_dphijmht", deltaPhi(ak4pfjets_p4().at(0).phi(), values_[mhtphi]), ak4_HTeta5()/ak4_HT(), evtweight_, cr.histMap, ";#Delta#phi(j1,MHT);HT5/HT", 64, 0, 3.2, 80, 1, 3);
        plot1d("h_ht5overht", ak4_HTeta5()/ak4_HT(), evtweight_, cr.histMap, ";HT5/HT", 80, 1, 3);
      }
    };
    if (suf == "") fillKineHists(suf);
    if (is_fastsim_ && suf == "" && (checkMassPt(1200, 50) || checkMassPt(800, 400) || checkMassPt(800, 675)))
      fillKineHists(Form("_%d_%d%s", mstop_, mlsp_, suf.c_str()));

    // Separate contribution by gen classification for background events
    if (doGenClassification && is_bkg_ && suf == "") {
      if (isZtoNuNu()) fillKineHists("_Znunu");
      else if (is2lep()) fillKineHists("_2lep");
      else if (is1lepFromW()) fillKineHists("_1lepW");
      else if (is1lepFromTop()) fillKineHists("_1lepTop");
      else fillKineHists("_unclass");  // either unclassified 1lep or 0lep, or something else unknown, shouldn't have (m)any
    }

    if (suf == "" && babyver >= 31.2 && ngoodjets() > 0) {
      plot2d("h2d_ht5_dphijmht", deltaPhi(ak4pfjets_p4().at(0).phi(), values_[mhtphi]), ak4_HTeta5()/ak4_HT(), evtweight_, cr.histMap, ";#Delta#phi(j1,MHT);HT5/HT", 32, 0, 3.2, 40, 1, 3);
      plot1d("h_ht5overht", ak4_HTeta5()/ak4_HT(), evtweight_, cr.histMap, ";HT5/HT", 80, 1, 3);
    }
  }*/
}

void StopLooper::fillHistosForCRemu(string suf, int trigType) {

}

////////////////////////////////////////////////////////////////////
// Functions that are not indispensible part of the main analysis

void StopLooper::fillTopTaggingHistos(string suffix) {

}

void StopLooper::testGenMatching(SR& sr) {

}

void StopLooper::testTopTaggingEffficiency(SR& sr) {

}

// Generate the yields for cut-flow table
void StopLooper::testCutFlowHistos(SR& sr) {
  /*
  if (gconf.is_data) return;

  // Defining cuts
  const vector<pair<string,bool(*)()>> cuts = {
    {"base", [](){ return (ngoodleps() >= 1 && ngoodjets() >= 2 && pfmet() > 150); }},
    {"mt", [](){ return (mt_met_lep() > 150); }},
    {"btag", [](){ return (ngoodbtags() >= 1); }},
    {"lepveto", [](){ return (nvetoleps() == 1); }},
    {"tauveto", [](){ return (PassTrackVeto() && PassTauVeto()); }},
    {"dphijmet", [](){ return (mindphi_met_j1_j2() > 0.8); }},
    {"met", [](){ return (pfmet() > 250); }},
  };
  const int ncuts = cuts.size();

  const vector<pair<string,bool(*)()>> cats = {
    {"tmod_gt0", [](){ return (topnessMod() > 0); }},
    {"tmod_gt10", [](){ return (topnessMod() > 10); }},
    {"Mlb_gt175", [](){ return (Mlb_closestb() > 175); }},
    {"Mlb_lt175", [](){ return (Mlb_closestb() > 0 && Mlb_closestb() < 175); }},
    {"mtag_ge1", [](){
      for (size_t iak8 = 0; iak8 < ak8pfjets_deepdisc_top().size(); ++iak8) {
        if (ak8pfjets_deepdisc_top()[iak8] > 0.4) return true;
      }
      return false;
    }},
    {"rtag_ge1", [](){
      for (float disc : tftops_disc()) {
        if (disc > 0.95) return true;
      }
      return false;
    }},
  };
  const int ncats = cats.size();

  const vector<pair<string,bool(*)()>> corcuts = {
    {"base", [](){ return (ngoodleps() >= 1 && ngoodjets() >= 3 && pfmet() > 150); }},
    {"mt", [](){ return (mt_met_lep() > 150); }},
    {"lepveto", [](){ return (nvetoleps() == 1); }},
    {"tauveto", [](){ return (PassTrackVeto() && PassTauVeto()); }},
    {"j1bveto", [](){ return (!ak4pfjets_passMEDbtag().at(0)); }},
    {"dphijmet", [](){ return (mindphi_met_j1_j2() > 0.5); }},
    {"lmet2d", [](){ return (lep1_p4().pt() < 50 || (lep1_p4().pt() < (250 - 100*lep1_dphiMET()))); }},
    {"met", [](){ return (pfmet() > 250); }},
  };
  const int ncorcuts = corcuts.size();

  const vector<pair<string,bool(*)()>> corcats = {
    {"medb", [](){ return (ngoodbtags() >= 1); }},
    {"medbn5jet", [](){ return (ngoodbtags() >= 1 && ngoodjets() >= 5); }},
    {"softb", [](){
      int nsoftb = nsoftbtags();
      for (auto sb : softtags_p4()) {
        if (nvetoleps() > 0 && isCloseObject(sb, lep1_p4(), 0.4)) nsoftb--;
        if (nvetoleps() > 1 && isCloseObject(sb, lep2_p4(), 0.4)) nsoftb--;
      }
      return nsoftb > 0;
    }},
  };
  const int ncorcats = corcats.size();

  string suf;
  if (isZtoNuNu()) suf = "_Znunu";
  else if (is2lep()) suf = "_2lep";
  else if (is1lepFromW()) suf = "_1lepW";
  else if (is1lepFromTop()) suf = "_1lepTop";
  else suf = "_unclass";  // either unclassified 1lep or 0lep, or something else unknown, shouldn't have (m)any

  if (is_fastsim_) {
    if (!(checkMassPt(1200, 100) || checkMassPt(850, 100) || checkMassPt(650, 350) || // T2tt fullsim available
          checkMassPt(1050,  50) || checkMassPt(800, 400) || checkMassPt(600, 400) || // New points
          checkMassPt(1050, 100) || checkMassPt(950, 100) || checkMassPt(750, 400) || // Signal points in SR figures
          checkMassPt(1000,  50) || checkMassPt(800, 400) || checkMassPt(500, 325) || // T2tt points in SUS-16-051
          checkMassPt(900,  50)  || checkMassPt(800, 350) || checkMassPt(500, 300) || // T2bW points in SUS-16-051
          checkMassPt(850,  50)  || checkMassPt(750, 300) || checkMassPt(450, 250) || // T2tb points in SUS-16-051
          checkMassPt(175,   1)  || checkMassPt(250,  50) || checkMassPt(250,  75) || // top corridor fullsim points
          checkMassPt(500, 325)  || checkMassPt(425, 338) || checkMassPt(400, 300) || // Corridor region points
          checkMassPt(508, 325)  || checkMassPt(492, 325) || checkMassPt(250, 100)
          )) return;
    suf = "_"+to_string(mstop_)+"_"+to_string(mlsp_);
  }

  for (int icut = 0; icut < ncuts; ++icut) {
    if (!cuts[icut].second()) goto corridor_cuts;
    plot1d("h_cutflow", icut, evtweight_, sr.histMap, ";Cuts", ncuts, 0, ncuts);
    plot1d("h_cutflow"+suf, icut, evtweight_, sr.histMap, ";Cuts", ncuts, 0, ncuts);
  }

  for (int icat = 0; icat < ncats; ++icat) {
    if (!cats[icat].second()) continue;
    plot1d("h_selcat", icat, evtweight_, sr.histMap, ";Selection Cats", ncats, 0, ncats);
    plot1d("h_selcat"+suf, icat, evtweight_, sr.histMap, ";Selection Cats", ncats, 0, ncats);
  }

corridor_cuts:
  for (int icorcut = 0; icorcut < ncorcuts; ++icorcut) {
    if (!corcuts[icorcut].second()) return;
    plot1d("h_corcutflow", icorcut, evtweight_, sr.histMap, ";cuts for corridor", ncorcuts, 0, ncorcuts);
    plot1d("h_corcutflow"+suf, icorcut, evtweight_, sr.histMap, ";cuts for corridor", ncorcuts, 0, ncorcuts);
  }

  for (int icorcat = 0; icorcat < ncorcats; ++icorcat) {
    if (!corcats[icorcat].second()) continue;
    plot1d("h_corselcat", icorcat, evtweight_, sr.histMap, ";Selection cats for corridor", ncorcats, 0, ncorcats);
    plot1d("h_corselcat"+suf, icorcat, evtweight_, sr.histMap, ";Selection cats for corridor", ncorcats, 0, ncorcats);
  }

  */
}
