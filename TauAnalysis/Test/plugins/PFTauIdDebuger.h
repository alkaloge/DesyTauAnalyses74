#ifndef TauAnalysis_Test_PFTauIdDebuger_h  
#define TauAnalysis_Test_PFTauIdDebuger_h
 
/** \class PFTauIdEffNtupleNtupleProducer2
 *
 * Produce an Ntuple of various quantities useful 
 * to check tau id. efficiencies and e/mu -> tau fake-rates
 *
 * \author Christian Veelken, LLR
 *
 * \version $Revision: 1.3 $
 *
 * $Id: PFTauIdDebuger.h,v 1.3 2012/03/08 10:31:49 veelken Exp $
 *
 */

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "FWCore/Utilities/interface/Exception.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/JetReco/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/Tau.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "TauAnalysis/CandidateTools/interface/candidateAuxFunctions.h"

#include <TMath.h>

#include <TTree.h>

#include <map>
#include <string>
#include <vector>
#include <ostream>

class PFTauIdDebuger : public edm::EDAnalyzer
{
 public:
  
  PFTauIdDebuger(const edm::ParameterSet&);
  ~PFTauIdDebuger();

  void analyze(edm::Event&, const edm::EventSetup&);
  void beginJob();
  const pat::Tau* findMatchingRecTau(const pat::TauCollection&, const reco::Candidate::LorentzVector&);
  const pat::Jet* findMatchingRecJet(const pat::JetCollection&, const reco::Candidate::LorentzVector&);

 private:

  std::string moduleLabel_;

  edm::InputTag srcGenParticles_;
  edm::InputTag srcGenJets_;
  edm::InputTag srcRecVetoElectrons_;
  edm::InputTag srcRecTaus_;
  edm::InputTag srcRecJets_;
  edm::InputTag srcVertices_;

  double minGenVisPt_;
  static int verbosity_;
};

#endif

PFTauIdDebuger::PFTauIdDebuger(const edm::ParameterSet& cfg) 
: moduleLabel_(cfg.getParameter<std::string>("@module_label")),
  ntuple_(0)
{
  srcGenParticles_ = cfg.getParameter<edm::InputTag>("srcGenParticles");
  srcGenJets_ = cfg.getParameter<edm::InputTag>("srcGenJets");
  srcRecTaus_ = cfg.getParameter<edm::InputTag>("srcRecTaus");
  srcRecJets_ = cfg.getParameter<edm::InputTag>("srcRecJets");
  srcVertices_ = cfg.getParameter<edm::InputTag>("srcVertices");

  minGenVisPt_ = 15.;
}

void PFTauIdDebuger::beginJob()
{ }

const pat::Tau* PFTauIdDebuger::findMatchingRecTau(const pat::TauCollection& recTaus, const reco::Candidate::LorentzVector& genParticleP4)
{
  const pat::Tau* recTau_matched = 0;

  double genTauDeltaR = 9.9;
  for ( pat::TauCollection::const_iterator recTau = recTaus.begin();
        recTau != recTaus.end(); ++recTau ) {

    if ( recTau->pt() < 15. ) continue;

    double dR = deltaR(recTau->p4(), genParticleP4);
    if ( dR < 0.3 && dR < genTauDeltaR ) {
      genTauDeltaR = dR;
      recTau_matched = &(*recTau);
    }
  }

  return recTau_matched;
}

const pat::Jet* PFTauIdDebuger::findMatchingRecJet(const pat::JetCollection& recJets, const reco::Candidate::LorentzVector& genParticleP4)
{
  const pat::Jet* recJet_matched = 0;

  double genTauDeltaR = 9.9;
  for ( pat::JetCollection::const_iterator recJet = recJets.begin();
	recJet != recJets.end(); ++recJet ) {

    if ( recJet->pt() < 15. ) continue;

    double dR = deltaR(recJet->p4(), genParticleP4);
    if ( dR < 0.3 && dR < genTauDeltaR ) {
      genTauDeltaR = dR;
      recJet_matched = &(*recJet);
    }
  }

  return recJet_matched;
}

void PFTauIdDebuger::analyze(edm::Event& evt, const edm::EventSetup& es) 
{
  
  edm::Handle<reco::GenParticleCollection> genParticles;
  evt.getByLabel(srcGenParticles_, genParticles);

  edm::Handle<pat::TauCollection> recTaus;
  evt.getByLabel(srcRecTaus_, recTaus);

  edm::Handle<pat::JetCollection> recJets;
  evt.getByLabel(srcRecJets_, recJets);

  //edm::Handle<reco::VertexCollection> vertices;
  //evt.getByLabel(srcVertices_, vertices);

  for ( reco::GenParticleCollection::const_iterator genParticle = genParticles->begin();
        genParticle != genParticles->end(); ++genParticle ) {

    if ( genParticle->pt() < minGenVisPt_ ) continue;

    unsigned numHypotheses = 0;

    bool isGenTau = false;
    for ( std::vector<int>::const_iterator pdgIdGenTau = pdgIdsGenTau.begin();
          pdgIdGenTau != pdgIdsGenTau.end(); ++pdgIdGenTau ) {
      if ( genParticle->status() == 2 && genParticle->pdgId() == (*pdgIdGenTau) ) isGenTau = true;
    }
    if ( isGenTau ) {
      reco::Candidate::LorentzVector genTauP4 = getVisMomentum(&(*genParticle));
      std::string genTauDecayMode_string = getGenTauDecayMode(&(*genParticle));
      int genTauDecayMode_int = -1;
      if      ( genTauDecayMode_string == "oneProng0Pi0"    ) genTauDecayMode_int = reco::PFTau::kOneProng0PiZero;
      else if ( genTauDecayMode_string == "oneProng1Pi0"    ) genTauDecayMode_int = reco::PFTau::kOneProng1PiZero;
      else if ( genTauDecayMode_string == "oneProng2Pi0"    ) genTauDecayMode_int = reco::PFTau::kOneProng2PiZero;
      else if ( genTauDecayMode_string == "threeProng0Pi0"  ) genTauDecayMode_int = reco::PFTau::kThreeProng0PiZero;
      else if ( genTauDecayMode_string == "threeProng1Pi0"  ) genTauDecayMode_int = reco::PFTau::kThreeProng1PiZero;
      //else if ( genTauDecayMode_string == "oneProngOther"   ||
      //          genTauDecayMode_string == "threeProngOther" ||
      //          genTauDecayMode_string == "rare"            ) genTauDecayMode_int = reco::PFTau::kRareDecayMode;
      if ( genTauDecayMode_int == reco::PFTau::kOneProng1PiZero && genTauP4.pt() > minGenVisPt_ ) {
        const pat::Tau* recTau_matched = findMatchingRecTau(*recTaus, genTauP4);
        double genTauDeltaR = ( recTau_matched ) ? deltaR(recTau_matched->p4(), genTauP4) : 9.9;
        bool genTauMatch = (genTauDeltaR < 0.3);
        int genTauDecayMode = genTauDecayMode_int;
	reco::Candidate::LorentzVector genTauLeadChHadP4 = getLeadChHadMomentum(&(*genParticle));
	if(!genTauMatch){
	  const pat::Jet* recJet_matched = findMatchingRecJet(*recJets, genTauP4);
	  double genTauDeltaR_Jet = ( recJet_matched ) ? deltaR(recJet_matched->p4(), genTauP4) : 9.9;
	  bool genTauMatch_Jet = (genTauDeltaR_Jet < 0.3);

	  if(genTauMatch_Jet){
	    
	    //match jet tracks to leading ch. aprticle and check its properties 
        
	std::vector<const reco::GenParticle*> genTauStableDaughters;
        findDaughters((&(*genParticle)), genTauStableDaughters, 1);
        //std::cout << "gen Tau Vis Pt "<<genTauP4.pt()<<std::endl;
        int nDaughters = genTauStableDaughters.size();
        int nDPhotons = 0;
        const reco::GenParticle* leadChParticle = 0;
        float lchPt = 0;
        for ( std::vector<const reco::GenParticle*>::const_iterator daughter = genTauStableDaughters.begin();
              daughter != genTauStableDaughters.end(); ++daughter ) {
          //std::cout << "daughter: pdgId = " << (*daughter)->pdgId() << ", Pt = " << (*daughter)->pt() << ","                                                                 
          //<< " eta = " << (*daughter)->eta() << ", phi = " << (*daughter)->phi()*180./TMath::Pi() << std::endl;
          if((*daughter)->pdgId() == 22)nDPhotons++;
          if ( (*daughter)->status() == 1 && (*daughter)->charge() != 0 ) {
            if((*daughter)->pt() > lchPt){
              lchPt = (*daughter)->pt();
              leadChParticle = (&(*(*daughter)));
            }
          }
        }
        int lchMother = 0;
        float lchDPhPt = 0;
        if(leadChParticle){
          lchMother = leadChParticle->mother()->pdgId();
          if(leadChParticle->mother()->pdgId() == leadChParticle->pdgId()){
            //std::cout<<"mother Id "<<leadChParticle->mother()->pdgId()<<" pt "<<leadChParticle->mother()->pt()<<" lch "<<leadChParticle->pdgId()<< "pt "<<leadChParticle->pt()<<std::endl;
            lchDPhPt = (leadChParticle->mother()->p4() - leadChParticle->p4()).pt();
          }
        }
        setGenTauExtraValues(lchMother, lchDPhPt, nDaughters, nDPhotons);

        ++numHypotheses;
      }
    }
  }

}
