#include "TauAnalysis/Test/plugins/DumpTriggerInfo.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DQMServices/Core/interface/DQMStore.h"



#include <TMath.h>

DumpTriggerInfo::DumpTriggerInfo(const edm::ParameterSet& cfg)
{
  srcL1GtReadoutRecord_ = cfg.getParameter<edm::InputTag>("srcL1GtReadoutRecord"); 
  srcHLTresults_        = cfg.getParameter<edm::InputTag>("srcHLTresults");
  srcFilters_           = cfg.getParameter<edm::InputTag>("srcFilters");  
}

void DumpTriggerInfo::printL1bits(const DecisionWord& l1GtDecision,
				   const L1GtTriggerMenu& l1GtTriggerMenu)
  {
    std::cerr << "Existing L1 bits:" << std::endl;
    const AlgorithmMap& l1GtAlgorithms = l1GtTriggerMenu.gtAlgorithmMap();
    for ( AlgorithmMap::const_iterator l1GtAlgorithm = l1GtAlgorithms.begin();
	  l1GtAlgorithm != l1GtAlgorithms.end(); ++l1GtAlgorithm ) {
      std::string l1BitName = l1GtAlgorithm->second.algoName();
      int index = l1GtAlgorithm->second.algoBitNumber();
      std::string triggerDecision = ( l1GtDecision[index] ) ? "passed" : "failed";
      std::cout << " L1 bit = " << l1BitName << ": " << triggerDecision << std::endl;
    }
  }

void DumpTriggerInfo::printHLTpaths(const edm::TriggerNames& triggerNames, const edm::TriggerResults& hltResults)
  {
    std::cerr << "Existing HLT paths:" << std::endl;
    for ( edm::TriggerNames::Strings::const_iterator triggerName = triggerNames.triggerNames().begin();
	  triggerName != triggerNames.triggerNames().end(); ++triggerName ) {
      unsigned int idx = triggerNames.triggerIndex(*triggerName);
      if ( idx < triggerNames.size() ) {
	std::string triggerDecision = ( hltResults.accept(idx) ) ? "passed" : "failed";
	std::cout << " HLT path = " << (*triggerName) << ": " << triggerDecision << std::endl;
	const std::vector<std::string> moduleNames = hltConfigProvider.moduleLabels(idx);
	//	unsigned int nModules = moduleNames.size();
	//for (unsigned int im=0; im<nModules; ++im) {
	//	  std::cout << "     " << moduleNames.at(im) << std::endl;
	//}
      }
    }
  }


void DumpTriggerInfo::beginRun(const edm::Run& iRun, const edm::EventSetup& iSetup)
{

  const std::string processName("HLT");
  bool changedConfig = false;

  hltConfigProvider.init(iRun, iSetup, processName,changedConfig); 

}

void DumpTriggerInfo::analyze(const edm::Event& evt, const edm::EventSetup& es)
{
  edm::Handle<L1GlobalTriggerReadoutRecord> l1GtReadoutRecord;
  evt.getByLabel(srcL1GtReadoutRecord_, l1GtReadoutRecord);
  
  edm::ESHandle<L1GtTriggerMenu> l1GtTriggerMenu;
  es.get<L1GtTriggerMenuRcd>().get(l1GtTriggerMenu);

  DecisionWord l1GtDecision = l1GtReadoutRecord->decisionWord();
  
  printL1bits(l1GtDecision, *l1GtTriggerMenu);

  edm::Handle<edm::TriggerResults> hltResults;
  evt.getByLabel(srcHLTresults_, hltResults);
  
  const edm::TriggerNames& triggerNames = evt.triggerNames(*hltResults);

  printHLTpaths(triggerNames, *hltResults);

  edm::Handle<trigger::TriggerEvent> triggerEventHandle;
  evt.getByLabel("hltTriggerSummaryAOD",triggerEventHandle);

  for (size_t iF = 0; iF < triggerEventHandle->sizeFilters(); ++iF)
    {
      const std::string filterName(triggerEventHandle->filterTag(iF).label());
      std::cout << "Filter " << filterName << std::endl;
    }


  //  std::map<std::string, std::string> hltMapping;



}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(DumpTriggerInfo);
