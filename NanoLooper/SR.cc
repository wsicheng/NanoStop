#include <stdexcept>
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

#include "SR.h"

using namespace std;

const pair<float,float> kemptycut{0,0};

void SR::SetName(string sr_name) {
  srname_ = sr_name;
}

void SR::SetDetailName(string detail_name) {
  detailname_ = detail_name;
}

void SR::SetVar(string var_name, float lower_bound, float upper_bound) {
  cuts_[var_name] = pair<float,float>(lower_bound, upper_bound);
}

void SR::SetVar(int ivar, float lower_bound, float upper_bound) {
  if (ivar >= (int)vcuts_.size()) vcuts_.resize(ivar+1, kemptycut);
  vcuts_[ivar] = pair<float,float>(lower_bound, upper_bound);
}

void SR::SetMETBins(std::vector<float> met_bins) {
  metbins_ = met_bins;
}

void SR::SetAllowDummyVars(int val) {
  kAllowDummyVars_ = val;
}

string SR::GetName() const {
  return srname_;
}

string SR::GetDetailName() const {
  return detailname_;
}

unsigned int SR::GetYield() const {
  return yield_;
}

float SR::GetLowerBound(string var_name) const {
  if (!cuts_.count(var_name))
    throw invalid_argument("Variable " + var_name + " not defined in SR " + srname_ + " ! Cannot get lower bound!");
  return (cuts_.at(var_name)).first;
}

float SR::GetUpperBound(string var_name) const {
  if (!cuts_.count(var_name))
    throw invalid_argument("Variable " + var_name + " not defined in SR " + srname_ + " ! Cannot get upper bound!");
  return (cuts_.at(var_name)).second;
}

float SR::GetLowerBound(int ivar) const {
  if (!VarExists(ivar))
    throw invalid_argument("Variable " + to_string(ivar) + " not defined in SR " + srname_ + " ! Cannot get lower bound!");
  return vcuts_[ivar].first;
}

float SR::GetUpperBound(int ivar) const {
  if (!VarExists(ivar))
    throw invalid_argument("Variable " + to_string(ivar) + " not defined in SR " + srname_ + " ! Cannot get upper bound!");
  return vcuts_[ivar].second;
}

unsigned int SR::GetNumberOfVariables() const {
  return cuts_.size();
}

unsigned int SR::GetNumberOfCuts() const {
  unsigned int ncut = vcuts_.size();
  for (const auto& ipair : vcuts_) {
    if (ipair == kemptycut) ncut--;
  }
  return ncut;
}

vector<string> SR::GetListOfVariables() const {
  vector<string> vars;
  for (auto it = cuts_.begin(); it != cuts_.end(); ++it) {
    vars.push_back(it->first);
  }
  return vars;
}

float* SR::GetMETBinsPtr() {
  return metbins_.data();
}

int SR::GetNMETBins() {
  return metbins_.size() - 1;
}

bool SR::PassesSelection(const vector<float>& values) {
  if ((kAllowDummyVars_ == 0 && GetNumberOfCuts() != values.size()) ||
      (kAllowDummyVars_ == 1 && GetNumberOfCuts()  > values.size())) {
    cout << "Number of variables to cut on != number of variables in signal region. Passed " << values.size() << ", expected " << GetNumberOfVariables() << endl;
    throw invalid_argument(srname_ + ": Number of variables to cut on != number of variables in signal region");
  }
  for (size_t icut = 0; icut < vcuts_.size(); icut++) {
    if (vcuts_[icut] == kemptycut) continue;
    float value = values[icut];
    if (!kAllowDummyVars_ && isnan(value)) {
      throw invalid_argument("The " + to_string(icut) + "th cut variable is not found in values");
    }
    float cut_lower = vcuts_[icut].first;
    float cut_upper = vcuts_[icut].second;
    if (value <  cut_lower) return false;
    if (value >= cut_upper) return false;
  }
  ++yield_;
  return true;
}

bool SR::PassesSelection(const map<string, float>& values) {
  const float ep = 0.000001;
  if ((kAllowDummyVars_ == 0 && GetNumberOfVariables() != values.size()) ||
      (kAllowDummyVars_ == 1 && GetNumberOfVariables()  > values.size())) {
    cout << "Number of variables to cut on != number of variables in signal region. Passed " << values.size() << ", expected " << GetNumberOfVariables() << endl;
    throw invalid_argument(srname_ + ": Number of variables to cut on != number of variables in signal region");
  }
  for (auto cit = cuts_.begin(); cit != cuts_.end(); cit++) {
    const auto vit = values.find(cit->first);
    if (vit != values.end()) { //check that we actually have bounds set for this variable
      float value = vit->second;
      float cut_lower = (cit->second).first;
      float cut_upper = (cit->second).second;
      if (value < cut_lower) return false;
      if (( abs(cut_upper + 1.0) > ep ) && (value >= cut_upper)) return false;
    }
    else if (!kAllowDummyVars_) {
      throw invalid_argument("Cut variable " + cit->first + " not found in values");
    }
  }
  ++yield_;
  return true;
}

int debug_print_count_SR_cc = 0;
const int k_debug_print_limit_SR_cc = 100;

bool SR::PassesSelectionPrintFirstFail(const map<string,float>& values) {
  const float ep = 0.000001;
  if ((kAllowDummyVars_ == 0 && GetNumberOfVariables() != values.size()) ||
      (kAllowDummyVars_ == 1 && GetNumberOfVariables()  > values.size())) {
    cout << "Number of variables to cut on != number of variables in signal region. Passed " << values.size() << ", expected " << GetNumberOfVariables() << endl;
    throw invalid_argument(srname_ + ": Number of variables to cut on != number of variables in signal region");
  }
  for (auto it = cuts_.begin(); it != cuts_.end(); it++) {
    const auto vit = values.find(it->first);
    if (vit != values.end()) { //check that we actually have bounds set for this variable
      float value = vit->second;
      float cut_lower = (it->second).first;
      float cut_upper = (it->second).second;
      // if (value < cut_lower) return false;
      // if (( abs(cut_upper + 1.0) > ep ) && (value >= cut_upper)) return false;
      const string k_suppress_vals = "met,dphijmet,mt";
      if (value < cut_lower) {
        if (k_suppress_vals.find(it->first) == string::npos && debug_print_count_SR_cc++ < k_debug_print_limit_SR_cc)
          cout << __LINE__ <<": var " << it->first << ": value= " << value << ", cut_lower= " << cut_lower << ": sr: " << srname_ << endl;
        return false;
      }
      if (( abs(cut_upper + 1.0) > ep ) && (value >= cut_upper)) {
        if (k_suppress_vals.find(it->first) == string::npos && debug_print_count_SR_cc++ < k_debug_print_limit_SR_cc)
          cout << __LINE__ <<": var " << it->first << ": value= " << value << ", cut_upper= " << cut_upper << ": sr: " << srname_ << endl;
        return false;
      }
    }
    else if (!kAllowDummyVars_) {
      throw invalid_argument("Cut variable " + it->first + " not found in values");
    }
  }
  return true;
}

bool SR::VarExists(std::string var_name) const {
  return cuts_.count(var_name);
}

bool SR::VarExists(int ivar) const {
  return (size_t(ivar) < vcuts_.size()) && (vcuts_[ivar] != kemptycut);
}

void SR::RemoveVar(string var_name) {
  if (cuts_.find(var_name) != cuts_.end()) cuts_.erase(var_name);
  else cerr << "WARNING: Variable " << var_name << " is not present in " << srname_ << ". Cannot remove!" << endl;
}

void SR::RemoveVar(int ivar) {
  if (size_t(ivar) < vcuts_.size()) vcuts_[ivar] = kemptycut;
  else cerr << "WARNING: The " << ivar << "th variable is not present in " << srname_ << ". Cannot remove!" << endl;
}

void SR::ReplaceVar(int ivar_old, int ivar_new) {
  if (!VarExists(ivar_old))
    throw invalid_argument("Variable " + to_string(ivar_old) + " not found in values. Cannot replace!");
  SetVar(ivar_new, GetLowerBound(ivar_old), GetUpperBound(ivar_old));
  RemoveVar(ivar_old);
}

void SR::Clear() {
  yield_ = 0;
  srname_ = "";
  detailname_ = "";
  cuts_.clear();
  defaultplots_.clear();
  metbins_.clear();
  kAllowDummyVars_ = false;
}


// --------------------------------------------------------------------------------
//                             Helper Functions
// --------------------------------------------------------------------------------

map<string,SRptrSet> generateSRptrSet(const vector<SR*>& SRptrVec)
{
  map<string,SRptrSet> var_to_SRset;
  map<string,set<tuple<float,float,SR*>>> master_map;

  cout << __FILE__ << __LINE__ << ": SRptrVec.size(): " << SRptrVec.size() << endl;
  for (const auto& isr : SRptrVec) {
    auto varlist = isr->GetListOfVariables();
    for (auto& var : varlist) {
      // cout << __FILE__ << __LINE__ << ": var: " << var << endl;
      float lower = isr->GetLowerBound(var), upper = isr->GetUpperBound(var);
      if (!master_map.count(var))
        master_map[var] = { {lower, upper, isr} };
      else
        master_map[var].emplace(lower, upper, isr);
    }
  }

  cout << __FILE__ << __LINE__ << ": master_map.size(): " << master_map.size() << endl;
  return var_to_SRset;

  for (auto it = master_map.begin(); it != master_map.end(); ++it) {
    set<float> edges;
    for (auto vlist : it->second) {
      edges.insert(get<0>(vlist));
      edges.insert(get<1>(vlist));
    }

    if (edges.erase(-1)) edges.insert(numeric_limits<float>::max());
    vector<float> bins(edges.begin(),edges.end());
    vector<set<SR*>> srsets(edges.size()-1);

    auto ib = bins.begin(), ie = bins.end();
    // sort(ib, ie);

    for (auto vlist : it->second) {
      auto il = lower_bound(ib, ie, get<0>(vlist));
      auto ir = lower_bound(ib, ie, get<1>(vlist));
      for (; il != ir; ++il) {
        srsets[il-ib].insert(get<2>(vlist));
      }
    }
    var_to_SRset[it->first] = {bins, srsets};
  }

  return var_to_SRset;
}


vector<SR*> findSatisfiedSRset(const map<string,float>& vars, const map<string,SRptrSet>& setsMap)
{
  vector<SR*> srset;
  bool firstvar = true;

  for (const auto& v : vars) {
    string var = v.first;
    float value = v.second;

    auto set_pair = setsMap.find(var);
    if (set_pair == setsMap.end()) throw std::invalid_argument("Var " + var + " wasn't found required in any SR of the SRvec");
    auto iset = set_pair->second;

    auto ibin_b = iset.bins.begin(), ibin_e = iset.bins.end();
    auto ibin = upper_bound(ibin_b, ibin_e, value);
    if (ibin == ibin_b || ibin == ibin_e) continue;

    unsigned int i = ibin - ibin_b - 1;
    if (firstvar) {
      srset = vector<SR*>(iset.sets[i].begin(), iset.sets[i].end());
      firstvar = false;
    } else {
      auto ibegin = srset.begin();
      auto iend = set_intersection(ibegin, srset.end(), iset.sets[i].begin(), iset.sets[i].end(), ibegin);
      if (iend == ibegin) break;
      srset.resize(iend - ibegin);
    }
  }

  return srset;
}
