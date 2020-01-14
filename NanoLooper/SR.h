#ifndef SR_h
#define SR_h

#include <string>
#include <map>
#include <set>
#include <vector>
#include <memory>

#include "TH1.h"

// class: SR [Stop Region]
// some more descriptions here...

class SR {

public:

  SR() : yield_(0), kAllowDummyVars_(0) {}
  SR(std::string name) : yield_(0), srname_(name), kAllowDummyVars_(0) {}
  SR(const SR& sr) = default;

  ~SR() {}

  void SetName(std::string sr_name);
  void SetDetailName(std::string detail_name);
  void SetVar(std::string var_name, float lower_bound, float upper_bound);
  void SetVar(int ivar, float lower_bound, float upper_bound);
  void SetMETBins(std::vector<float> met_bins);
  void SetAllowDummyVars(int val);

  std::string GetName() const;
  std::string GetDetailName() const;
  unsigned int GetYield() const;
  float GetLowerBound(std::string var_name) const;
  float GetUpperBound(std::string var_name) const;
  float GetLowerBound(int ivar) const;
  float GetUpperBound(int ivar) const;
  unsigned int GetNumberOfVariables() const;
  unsigned int GetNumberOfCuts() const;
  std::vector<std::string> GetListOfVariables() const;
  std::vector<int> GetListOfVariableIndicies() const;
  int GetNMETBins();
  float* GetMETBinsPtr();

  bool PassesSelection(const std::map<std::string,float>& values);
  bool PassesSelection(const std::vector<float>& values);
  bool VarExists(std::string var_name) const;
  bool VarExists(int ivar) const;
  void RemoveVar(std::string var_name);
  void RemoveVar(int ivar);
  void ReplaceVar(std::string var_old, std::string var_new);
  void ReplaceVar(int ivar_old, int ivar_new);
  void Clear();

  // Function for debug
  bool PassesSelectionPrintFirstFail(const std::map<std::string,float>& values);

  // used for plotting
  std::map<std::string, TH1*> histMap;

private:

  unsigned int yield_;
  std::string srname_;
  std::string detailname_;
  std::map<std::string,std::pair<float,float>> cuts_;
  std::vector<std::pair<float,float>> vcuts_;
  std::vector<std::string> defaultplots_;
  std::vector<float> metbins_;
  int kAllowDummyVars_;

};


// Helper Class & Functions <-- Under development

struct SRptrSet {
  std::vector<float> bins;
  std::vector<std::set<SR*>> sets;
};

std::map<std::string,SRptrSet> generateSRptrSet(const std::vector<SR*>& SRptrVec);

std::vector<SR*> findSatisfiedSRset(const std::map<std::string,float>& vars, const std::map<std::string,SRptrSet>& setsMap);


#endif
