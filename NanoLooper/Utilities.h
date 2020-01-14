#ifndef NANOLOOPER_UTILITIES_H
#define NANOLOOPER_UTILITIES_H

#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"

#include <map>
#include <string>

namespace utils {

// Histogram manipulation
inline void moveOverFlowToLastBin1D(TH1* hist) {
  int nbin = hist->GetNbinsX();
  if (hist->GetBinLowEdge(nbin+1) < 1499.9) return;  // only when the last bin is the infinity bin
  if (hist->GetBinContent(nbin+1) > 0) {
    // cout << "Moving the overflow for hist: " << hist->GetTitle() << " to its last bin!" << endl;
    double err = 0;
    hist->SetBinContent(nbin, hist->IntegralAndError(nbin, -1, err));
    hist->SetBinError(nbin, err);
    hist->SetBinContent(nbin+1, 0);
    hist->SetBinError(nbin+1, 0);
  }
}

inline void moveXOverFlowToLastBin3D(TH3* hist) {
  int nbinx = hist->GetNbinsX();
  if (hist->GetXaxis()->GetBinLowEdge(nbinx+1) < 1499.9) return; // only when the last bin is infinity
  for (int ibiny = 0; ibiny < hist->GetNbinsY(); ++ibiny) {
    for (int ibinz = 0; ibinz < hist->GetNbinsZ(); ++ibinz) {
      if (hist->GetBinContent(nbinx+1, ibiny, ibinz) <= 0) continue;
      double err = 0;
      double yield = hist->IntegralAndError(nbinx, -1, ibiny, ibiny, ibinz, ibinz, err);
      hist->SetBinContent(nbinx, ibiny, ibinz, yield);
      hist->SetBinError(nbinx, ibiny, ibinz, err);
      hist->SetBinContent(nbinx+1, ibiny, ibinz, 0);
      hist->SetBinError(nbinx+1, ibiny, ibinz, 0);
    }
  }
}

inline void zeroOutNegativeYields(TH1* hist) {
  int nbin = hist->GetNbinsX();
  for (int ibin = 1; ibin <= nbin; ++ibin) {
    if (hist->GetBinContent(ibin) < 0) {
      if (string(hist->GetName()).find("h_metbins_") == string::npos)  // only print out for central hist
        cout << "Reverting negative yield " << hist->GetBinContent(ibin) << " in: " << hist->GetTitle() << " bin " << ibin << " to 0!" << endl;
      hist->SetBinContent(ibin, 0);
      // hist->SetBinError(ibin, 0); // should we set the error to 0 also?
    }
  }
}

// Quick helper functions
inline float deltaPhi(float phi1, float phi2) {
  float dphi = fabs(phi1 - phi2);
  if (dphi > TMath::Pi()) dphi = TMath::TwoPi() - dphi;
  return dphi;
}

inline float min2deltaPhi(float phi0, float phi1, float phi2) {
  return std::min(deltaPhi(phi0, phi1), deltaPhi(phi0, phi2));
}

inline float calculateMT(double Et1, double phi1, double Et2, double phi2) {
  return sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
}

// Old functions that enforce float for ranges to be consistent with xval for floating point errors
void plot1D(string name, float xval, double weight, std::map<string, TH1*> &allhistos, string title, int numbinsx, float xmin, float xmax)
{
  if (title=="") title=name;
  std::map<string, TH1*>::iterator iter= allhistos.find(name);
  if (iter == allhistos.end()) { //no histo for this yet, so make a new one
    TH1D* currentHisto= new TH1D(name.c_str(), title.c_str(), numbinsx, xmin, xmax);
    currentHisto->SetDirectory(0);
    currentHisto->Sumw2();
    currentHisto->Fill(xval, weight);
    allhistos.insert(std::pair<string, TH1*>(name, currentHisto) );
  } else {
    iter->second->Fill(xval, weight);
  }
}

void linkHist(string hnew, string hexist, std::map<std::string, TH1*> &allhistos)
{
  // Could be useful when mutiple ratio hists sharing a common denominator
  if (allhistos.count(hnew)) return;
  auto iter = allhistos.find(hexist);
  if (iter == allhistos.end()) throw std::logic_error("linkHist(): Histogram "+hexist+" need to be plotted first");
  allhistos.insert( std::pair<std::string, TH1*>(hnew, iter->second) );
}


// Templated function
template<class LorentzVectorType>
inline float calculateMt(const LorentzVectorType& l1p4, const LorentzVectorType& l2p4, double met, double met_phi){
  LorentzVectorType zp4 = l1p4 + l2p4;
  float phi1 = zp4.Phi();
  float phi2 = met_phi;
  float Et1  = zp4.Et();
  float Et2  = met;

  return sqrt(2*Et1*Et2*(1.0 - cos(phi1-phi2)));
}

template<class LorentzVectorType>
bool isCloseObject(const LorentzVectorType& p1, const LorentzVectorType& p2, const float conesize = 0.4, float* deltaR = nullptr)
{
  const float PI = TMath::Pi();
  float deltaEta = fabs(p1.eta() - p2.eta());
  if (deltaEta > conesize) return false;
  float deltaPhi = fabs(p1.phi() - p2.phi());
  if (deltaPhi > PI) deltaPhi = 2*PI - deltaPhi;
  if (deltaPhi > conesize) return false;
  float deltaR2 = deltaEta*deltaEta + deltaPhi*deltaPhi;
  if (deltaR2 > conesize*conesize) return false;
  if (deltaR) *deltaR = sqrt(deltaR2);

  return true;
}

template<typename... TArgs>
void plot1d(std::string name, double xval, double weight, std::map<std::string, TH1*> &allhistos, TArgs... args)
{
  auto iter = allhistos.find(name);
  if (iter == allhistos.end()) {
    TH1D* currentHisto= new TH1D(name.c_str(), args...);
    currentHisto->SetDirectory(0);
    currentHisto->Sumw2();
    currentHisto->Fill(xval, weight);
    allhistos.insert( std::pair<std::string, TH1*>(name, currentHisto) );
  } else {
    iter->second->Fill(xval, weight);
  }
}

template<typename... TArgs>
void plot2d(std::string name, double xval, double yval, double weight, std::map<std::string, TH1*> &allhistos, TArgs... args)
{
  auto iter = allhistos.find(name);
  if (iter == allhistos.end()) {
    TH2D* currentHisto= new TH2D(name.c_str(), args...);
    currentHisto->SetDirectory(0);
    currentHisto->Sumw2();
    currentHisto->Fill(xval, yval, weight);
    allhistos.insert( std::pair<std::string, TH1*>(name, currentHisto) );
  } else {
    ((TH2D*) iter->second)->Fill(xval, yval, weight);
  }
}

template<typename... TArgs>
void plot3d(std::string name, double xval, double yval, double zval, double weight, std::map<std::string, TH1*> &allhistos, TArgs... args)
{
  auto iter = allhistos.find(name);
  if (iter == allhistos.end()) {
    TH3D* currentHisto= new TH3D(name.c_str(), args...);
    currentHisto->SetDirectory(0);
    currentHisto->Sumw2();
    currentHisto->Fill(xval, yval, zval, weight);
    allhistos.insert( std::pair<std::string, TH1*>(name, currentHisto) );
  } else {
    ((TH3D*) iter->second)->Fill(xval, yval, zval, weight);
  }
}

}

#endif  // HZZLOOPER_UTILITIES_H
