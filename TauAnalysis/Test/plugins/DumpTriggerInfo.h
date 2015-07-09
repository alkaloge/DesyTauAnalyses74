#ifndef TauAnalysis_Test_DumpTriggerInfo_h
#define TauAnalysis_Test_DumpTriggerInfo_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"

#include "DQMServices/Core/interface/MonitorElement.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/Common/interface/Handle.h"

class DumpTriggerInfo : public edm::EDAnalyzer 
{
 public:

  explicit DumpTriggerInfo(const edm::ParameterSet&);
  ~DumpTriggerInfo() {}
    
  void analyze(const edm::Event&, const edm::EventSetup&);
  void beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup);

 private:

  edm::InputTag srcL1GtReadoutRecord_; 
  //edm::InputTag srcL1GtObjectMaps_;
  edm::InputTag srcHLTresults_; 

  edm::InputTag srcFilters_;

  HLTConfigProvider hltConfigProvider;

  void  printL1bits(const DecisionWord& l1GtDecision,
		    const L1GtTriggerMenu& l1GtTriggerMenu);

  void printHLTpaths(const edm::TriggerNames& triggerNames, const edm::TriggerResults& hltResults);


};

#endif
