// -*- C++ -*-
//
// Package:    ZcandHisto
// Class:      ZcandHisto
// 
/**\class ZcandHisto ZcandHisto.cc AWhitbeck/ZcandHisto/src/ZcandHisto.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Andrew Whitbeck
//         Created:  Wed Oct 16 10:10:37 CDT 2013
// $Id$
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TH1.h"
#include "TH2.h"
#include "TLorentzVector.h"

#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>

#include <vector>



class ZcandHisto : public edm::EDAnalyzer {

public:
  explicit ZcandHisto(const edm::ParameterSet&);
  ~ZcandHisto();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  virtual void beginRun(edm::Run const&, edm::EventSetup const&);
  virtual void endRun(edm::Run const&, edm::EventSetup const&);
  virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
  
  // ----------member data ---------------------------
  TH1D * Zmass_histo; 
  TH1D * Zpt_histo; 
  TH1D * Hmass_histo; 
  TH1D * Hpt_histo; 
  TH1D * telHmass_histo; 
  TH1D * telHpt_histo; 
  TH2D * pfCand_histo;

  // ---------- configurable data ----------------
  // --------------- members ---------------------
  
  std::string ZCompCollection;    // name of di-lepton collection
  std::string pfCandCollection;      // name of PF candidate collection

};


ZcandHisto::ZcandHisto(const edm::ParameterSet& iConfig):
  ZCompCollection(iConfig.getUntrackedParameter<std::string>("ZCompCollection","ZCandidate")),
  pfCandCollection(iConfig.getUntrackedParameter<std::string>("pfCandCollection","pfNoElectronPFlow"))
{
   //now do what ever initialization is needed
  edm::Service<TFileService> fs;
  Zmass_histo = fs->make<TH1D>("Zmass" , "m_{Z}"   , 200 , 0.   , 200. );
  Zpt_histo   = fs->make<TH1D>("Zpt"   , "p_{T,Z}" , 200 , 0.   , 200. );
  Hmass_histo = fs->make<TH1D>("Hmass" , "m_{H}"   , 200 , 0.   , 200. );
  Hpt_histo   = fs->make<TH1D>("Hpt" , "p_{T,H}"   , 200 , 0.   , 200. );

  telHmass_histo = fs->make<TH1D>("telHmass" , "m_{H,tel}" , 200 , 0. , 200. );
  telHpt_histo   = fs->make<TH1D>("telHpt"   , "p_{T,tel}" , 200 , 0. , 200. );

  pfCand_histo = fs->make<TH2D>("pfCand" , ";#Phi;#eta;Energy"   , 300 , -3.1415 , 3.1415,
				                                   300 , -3.     , 3.);

}


ZcandHisto::~ZcandHisto()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZcandHisto::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // fill histograms for di-lepton system

  using namespace edm;

  //  sanity check on muon collection

  Handle< View<pat::Muon> > muons;
  iEvent.getByLabel("selectedPatMuons",muons);
  
  std::cout << "number of muons: " << muons->size() << std::endl;

  Handle< View<reco::CompositeCandidate> > Zcands;
  iEvent.getByLabel(ZCompCollection,Zcands);
  
  double bestMass = -100000.0;
  double bestPt   = -100000.0;
  
  for(View<reco::CompositeCandidate>::const_iterator iCand = Zcands->begin();
      iCand != Zcands->end();
      ++iCand){
    
    if ( fabs( iCand->mass() - 91.188 ) < fabs ( bestMass - 91.188 ) ){
      bestMass = iCand->mass();
      bestPt   = iCand->pt();
    }
    
    
  }

  if(bestPt < 120.){

    std::cout << "Z candidate didn't pass pt cut" << std::endl;
    std::cout << "Z mass: " << bestMass << std::endl;
    std::cout << "Z pt:   " << bestPt   << std::endl;
    return;

  }
  
  Zmass_histo->Fill( bestMass );
  Zpt_histo->  Fill( bestPt    );

  /*
  
  // cluster jets and find jet axies 
  // with ak-5 algorithm

  Handle< View<reco::PFCandidate> > pfCands;
  iEvent.getByLabel(pfCandCollection,pfCands);

  std::vector<fastjet::PseudoJet> PFparticles;

  std::vector<fastjet::PseudoJet> jets_ak4;

  for(View<reco::PFCandidate>::const_iterator iCand = pfCands->begin();
      iCand != pfCands->end();
      ++iCand){

    PFparticles.push_back( fastjet::PseudoJet( iCand->px(), 
					      iCand->py(),
					      iCand->pz(),
					      iCand->energy()
					      ));

    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, 0.4);
    fastjet::ClusterSequence cs(PFparticles, jet_def);
    jets_ak4 = sorted_by_pt(cs.inclusive_jets());

  }

  TLorentzVector j1(jets_ak4[0].px(),
		    jets_ak4[0].py(),
		    jets_ak4[0].pz(),
		    jets_ak4[0].E()
		    );

  TLorentzVector j2(jets_ak4[1].px(),
		    jets_ak4[1].py(),
		    jets_ak4[1].pz(),
		    jets_ak4[1].E()
		    );
		    
  Hmass_histo->Fill((j1+j2).M());

  */

  Handle< View<pat::Jet> > jetColl;
  iEvent.getByLabel("selectedPatJetsPFlow",jetColl);

  std::vector<TLorentzVector> jets;

  for(View<pat::Jet>::const_iterator iJet = jetColl->begin();
      iJet != jetColl->end();
      ++iJet){

    jets.push_back( TLorentzVector(iJet->px(),
				   iJet->py(),
				   iJet->pz(),
				   iJet->energy()));
    
  }  

  //std::cout << "leading jet mass: " << j1.M() << std::endl;
  //std::cout << "sub-leading jet mass: " << j2.M() << std::endl;

  Hmass_histo->Fill( (jets[0]+jets[1]).M() );
  Hpt_histo->Fill( (jets[0]+jets[1]).Pt() );

  Handle< View<reco::PFCandidate> > pfCands;
  iEvent.getByLabel(pfCandCollection,pfCands);

  // loop over various radii
  
  for(double R=0.2; R<1.5; R+=.01){

    TLorentzVector teljet1(0.0,0.0,0.0,0.0);
    TLorentzVector teljet2(0.0,0.0,0.0,0.0);

    // sum pf-candidate 4-vectors

    for(View<reco::PFCandidate>::const_iterator iCand = pfCands->begin();
	iCand != pfCands->end();
	++iCand){
      
      // delta-R squared w/sp to leading jet
      double dR1sq = pow( jets[0].Eta() - iCand->eta() , 2 )
	+ pow( jets[0].Phi() - iCand->phi() , 2 );

      //std::cout << "cand dR1 squared: " << dR1sq << std::endl;

      // delta-R squared w/sp to sub-leading jet
      double dR2sq = pow( jets[1].Eta() - iCand->eta() , 2 )
	+ pow( jets[1].Phi() - iCand->phi() , 2 );

      //std::cout << "cand dR2 squared: " << dR2sq << std::endl;

      // check if pf-candidate is in jet cone
      if( dR1sq < R*R && dR1sq<dR2sq){

	//std::cout << "found cand for j1" << std::endl;
	teljet1+=TLorentzVector(iCand->px(),
				iCand->py(),
				iCand->pz(),
				iCand->energy());
      }

      if( dR2sq < R*R && dR2sq<dR1sq){

	//std::cout << "found cand for j2" << std::endl;
	teljet2+=TLorentzVector(iCand->px(),
				iCand->py(),
				iCand->pz(),
				iCand->energy());
      }

      pfCand_histo->Fill(iCand->phi(),iCand->eta(),iCand->energy());

    }// end loop over pf-candidates

    //std::cout << "R: " << R << "mass: " << (teljet1+teljet2).M() << std::endl;

    telHmass_histo->Fill( (teljet1+teljet2).M() );
    telHpt_histo->Fill( (teljet1+teljet2).Pt() );

  }// end loop over jet cones

  //std::cout << "default mass: " << (jets[0]+jets[1]).M() << std::endl;
  //std::cout << "jet1 (eta,phi): " << jets[0].Eta() << ", " << jets[0].Phi() << std::endl;
  //std::cout << "jet2 (eta,phi): " << jets[1].Eta() << ", " << jets[1].Phi() << std::endl; 

}


// ------------ method called once each job just before starting event loop  ------------
void 

ZcandHisto::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZcandHisto::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ZcandHisto::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ZcandHisto::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ZcandHisto::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ZcandHisto::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZcandHisto::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {

  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

}

//define this as a plug-in
DEFINE_FWK_MODULE(ZcandHisto);
