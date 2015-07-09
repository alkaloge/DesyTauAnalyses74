#ifndef TauAnalysis_RecoTools_VertexMultiplicityReweightProducer_h
#define TauAnalysis_RecoTools_VertexMultiplicityReweightProducer_h

/** \class VertexMultiplicityReweightProducer
 *
 * Reweight Monte Carlo events simulated with pile-up to match the
 * vertex multiplicity distribution observed in the analyzed data sample
 *
 * NOTE:
 *      o weight > 1: fraction of events which given vertex multiplicity higher in data than in (pile-up) MC
 *      o weight = 1: fraction of events which given vertex multiplicity same   in data than in (pile-up) MC
 *      o weight < 1: fraction of events which given vertex multiplicity lower  in data than in (pile-up) MC
 *
 * \authors Christian Veelken
 *
 * \version $Revision: 1.2 $
 *
 * $Id: VertexMultiplicityReweightProducer.h,v 1.2 2011/04/13 17:05:13 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TauAnalysis/BgEstimationTools/interface/ObjValExtractorBase.h"

#include <string>

class VertexMultiplicityReweightProducer : public edm::EDProducer 
{
 public:
  explicit VertexMultiplicityReweightProducer(const edm::ParameterSet&);
  ~VertexMultiplicityReweightProducer();

  void produce(edm::Event&, const edm::EventSetup&);

 private:
  std::string moduleLabel_;

  ObjValExtractorBase* extractor_;
};

#endif

