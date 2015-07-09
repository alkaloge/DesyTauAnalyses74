#include "TauAnalysis/RecoTools/plugins/VertexMultiplicityReweightProducer.h"

VertexMultiplicityReweightProducer::VertexMultiplicityReweightProducer(const edm::ParameterSet& cfg)
  : moduleLabel_(cfg.getParameter<std::string>("@module_label")),
    extractor_(0)
{
  std::string pluginType = "VertexMultiplicityReweightExtractor";
  extractor_ = ObjValExtractorPluginFactory::get()->create(pluginType, cfg);

  produces<double>("");
}

VertexMultiplicityReweightProducer::~VertexMultiplicityReweightProducer() 
{
  delete extractor_;
}

void VertexMultiplicityReweightProducer::produce(edm::Event& evt, const edm::EventSetup& es)
{
  //std::cout << "<VertexMultiplicityReweightProducer::produce (moduleLabel = " << moduleLabel_ << ")>" << std::endl;
  std::auto_ptr<double> weightPtr(new double((*extractor_)(evt)));
  evt.put(weightPtr);
}

#include "FWCore/Framework/interface/MakerMacros.h"

DEFINE_FWK_MODULE(VertexMultiplicityReweightProducer);


