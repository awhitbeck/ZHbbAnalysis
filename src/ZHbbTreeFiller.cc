// -*- C++ -*-
//
// Package:    ZHbbTreeFiller
// Class:      ZHbbTreeFiller
// 
/**\class ZHbbTreeFiller ZHbbTreeFiller.cc AWhitbeck/ZHbbTreeFiller/src/ZHbbTreeFiller.cc
 
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
#include "TTree.h"

#include <DataFormats/PatCandidates/interface/Jet.h>
#include <DataFormats/PatCandidates/interface/Muon.h>
#include <DataFormats/PatCandidates/interface/CompositeCandidate.h>
#include <DataFormats/ParticleFlowCandidate/interface/PFCandidate.h>

#include <fastjet/ClusterSequence.hh>
#include <fastjet/PseudoJet.hh>
#include <fastjet/JetDefinition.hh>

#include <vector>



class ZHbbTreeFiller : public edm::EDAnalyzer {
    
public:
    explicit ZHbbTreeFiller(const edm::ParameterSet&);
    ~ZHbbTreeFiller();
    
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
    
    TTree* ZHbbTree;
    
    struct treeStructure {
        
        double Zmass;
        double Zpt;
        double lep1_pt;
        double lep2_pt;
        double lep1_eta;
        double lep2_eta;
        std::vector<double> Hmass;
        std::vector<double> Hpt;
        std::vector<double> jet1_pt;
        std::vector<double> jet2_pt;
        std::vector<double> jet1_eta;
        std::vector<double> jet2_eta;
        std::vector<double> R;
        double fracPass;
        
    };
    
    treeStructure myTree;
    
    // ---------- configurable data ----------------
    // --------------- members ---------------------
    
    std::string ZCompCollection;    // name of di-lepton collection
    std::string pfCandCollection;      // name of PF candidate collection
    bool isGEN_;

};


ZHbbTreeFiller::ZHbbTreeFiller(const edm::ParameterSet& iConfig):
ZCompCollection(iConfig.getUntrackedParameter<std::string>("ZCompCollection","ZCandidate")),
pfCandCollection(iConfig.getUntrackedParameter<std::string>("pfCandCollection","pfNoElectronPFlow")),
isGEN_(iConfig.getUntrackedParameter<bool>("isGen",true))
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
    
    ZHbbTree = fs->make<TTree>("ZHbbTree","ZHbbTree");
    
    ZHbbTree->Branch("Zmass",&myTree.Zmass,"Zmass/D");
    ZHbbTree->Branch("Zpt",&myTree.Zpt,"Zpt/D");
    ZHbbTree->Branch("lep1_pt",&myTree.lep1_pt,"lep1_pt/D");
    ZHbbTree->Branch("lep2_pt",&myTree.lep2_pt,"lep2_pt/D");
    ZHbbTree->Branch("lep1_eta",&myTree.lep1_eta,"lep1_eta/D");
    ZHbbTree->Branch("lep2_eta",&myTree.lep2_eta,"lep2_eta/D");
    ZHbbTree->Branch("Hmass",&myTree.Hmass);
    ZHbbTree->Branch("Hpt",&myTree.Hpt);
    ZHbbTree->Branch("jet1_pt",&myTree.jet1_pt);
    ZHbbTree->Branch("jet2_pt",&myTree.jet2_pt);
    ZHbbTree->Branch("jet1_eta",&myTree.jet1_eta);
    ZHbbTree->Branch("jet2_eta",&myTree.jet2_eta);
    ZHbbTree->Branch("R",&myTree.R);
    ZHbbTree->Branch("fracPass",&myTree.fracPass,"fracPass/D");
    
}


ZHbbTreeFiller::~ZHbbTreeFiller()
{
    
    // do anything here that needs to be done at desctruction time
    // (e.g. close files, deallocate resources etc.)
    
}


//
// member functions
//

// ------------ method called for each event  ------------
void
ZHbbTreeFiller::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////    
    // Do Z stuff
    
    using namespace edm;
    
    Handle< View<reco::CompositeCandidate> > Zcands;
    iEvent.getByLabel(ZCompCollection,Zcands);
    
    double bestMass     = -100000.0;
    double bestPt       = -100000.0;
    double bestLep1_pt  = -100000.0;
    double bestLep2_pt  = -100000.0;
    double bestLep1_eta = -100000.0;
    double bestLep2_eta = -100000.0;
    
    for(View<reco::CompositeCandidate>::const_iterator iCand = Zcands->begin();
        iCand != Zcands->end();
        ++iCand){
        
        if ( fabs( iCand->mass() - 91.188 ) < fabs ( bestMass - 91.188 ) ){
            bestMass = iCand->mass();
            bestPt   = iCand->pt();
            bestLep1_pt = iCand->daughter(0)->pt();
            bestLep1_eta = iCand->daughter(0)->eta();
            bestLep2_pt = iCand->daughter(1)->pt();
            bestLep2_eta = iCand->daughter(1)->eta();
        }
        
        
    }
    
    if(bestPt < 120.){
        
        //std::cout << "Z candidate didn't pass pt cut" << std::endl;
        //std::cout << "Z mass: " << bestMass << std::endl;
        //std::cout << "Z pt:   " << bestPt   << std::endl;
        return;
        
    }
    
    Zmass_histo->Fill( bestMass );
    Zpt_histo->  Fill( bestPt    );
    
    // jet stuff...
    //std::cout << "jet stuff..." << std::endl;
    
    // reset vectors  
    myTree.Hmass.clear();
    myTree.Hpt.clear();
    myTree.R.clear();
    myTree.jet1_pt.clear();
    myTree.jet1_eta.clear();
    myTree.jet2_pt.clear();
    myTree.jet2_eta.clear();
    
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////    
    // Get Particles
    
    std::vector< fastjet::PseudoJet > particles;
    if (isGEN_){
        Handle < reco::GenParticleRefVector > gens_nonu;
        iEvent.getByLabel("genParticlesForJetsNoNu",gens_nonu);
        
        for (unsigned int i = 0; i < gens_nonu->size(); i++){
            const reco::GenParticle& P = *(gens_nonu->at(i));
            fastjet::PseudoJet tmp_psjet( P.px(), P.py(), P.pz(), P.energy() );
            particles.push_back( tmp_psjet );
        }
    }       
    
    
    else{
        
        Handle< View<reco::PFCandidate> > pfCandsCHS;
        iEvent.getByLabel("pfNoElectronPFlow",pfCandsCHS);
        Handle< View<reco::PFCandidate> > pfCandsPU;
        iEvent.getByLabel("pfPileUpPFlow",pfCandsPU);
        
        for (unsigned int i = 0; i < pfCandsCHS->size(); i++){
            const reco::PFCandidate P = (pfCandsCHS->at(i));
            fastjet::PseudoJet tmp_psjet( P.px(), P.py(), P.pz(), P.energy() );
            particles.push_back( tmp_psjet );
        }
        for (unsigned int i = 0; i < pfCandsPU->size(); i++){
            const reco::PFCandidate P = (pfCandsPU->at(i));
            fastjet::PseudoJet tmp_psjet( P.px(), P.py(), P.pz(), P.energy() );
            particles.push_back( tmp_psjet );
        }        
    }
    
    
    ////////////////////////////////////////////////
    ////////////////////////////////////////////////    
    // A quick clustering for getting telescoping seeds
    fastjet::JetDefinition jd_ak5(fastjet::antikt_algorithm, 0.5);
    fastjet::ClusterSequence* thisClustering = new fastjet::ClusterSequence(particles, jd_ak5);
    std::vector<fastjet::PseudoJet> jets = sorted_by_pt(thisClustering->inclusive_jets(15.0));
    
    
    //    // Get jets from standard PAT workflow to use
    //    // as seed for telescoping
    //    Handle< View<pat::Jet> > jetColl;
    //    iEvent.getByLabel("selectedPatJetsPFlow",jetColl);
    //    
    //    std::vector<TLorentzVector> jets;
    //    
    //    for(View<pat::Jet>::const_iterator iJet = jetColl->begin();
    //        iJet != jetColl->end();
    //        ++iJet){
    //        
    //        jets.push_back( TLorentzVector(iJet->px(),
    //                                       iJet->py(),
    //                                       iJet->pz(),
    //                                       iJet->energy()));
    //        
    //    }  
    //    
    //    Hmass_histo->Fill( (jets[0]+jets[1]).M() );
    //    Hpt_histo->Fill( (jets[0]+jets[1]).Pt() );
    
    // this rho for simplicity
    Handle< double > rhoSW;
    iEvent.getByLabel("kt6PFJetsPFlow","rho",rhoSW);
    
    // also get rho
    double rho = *rhoSW;
    
    //std::cout << "got PFCandidate collection" << std::endl;
    
    // loop over various radii
    
    double beginR=0.3;
    double endR  =1.3;
    int numR  =100;
    
    myTree.fracPass=0.0;
    
    std::cout << "jets[0].pt() = " << jets[0].pt() << ", eta = " << jets[0].eta() << std::endl;
    std::cout << "jets[1].pt() = " << jets[1].pt() << ", eta = " << jets[1].eta() << std::endl;
    
    for(double R=beginR; R<endR; R += (endR-beginR)/numR ){
        
        //std::cout << "R: " << R << std::endl;
        
        TLorentzVector teljet1_raw(0.0,0.0,0.0,0.0);
        TLorentzVector teljet2_raw(0.0,0.0,0.0,0.0);
        
        // sum pf-candidate 4-vectors
        
        for(unsigned int i = 0; i < particles.size(); i++){
            
            // delta-R squared w/sp to leading jet
            double dR1sq = pow( jets[0].eta() - particles[i].eta() , 2 )
            + pow( jets[0].phi() - particles[i].phi() , 2 );
            
            // delta-R squared w/sp to sub-leading jet
            double dR2sq = pow( jets[1].eta() - particles[i].eta() , 2 )
            + pow( jets[1].phi() - particles[i].phi() , 2 );
            
            // check if pf-candidate is in jet cone
            if( dR1sq < R*R && dR1sq<dR2sq){
                
                //std::cout << "found cand for j1" << std::endl;
                teljet1_raw+=TLorentzVector(particles[i].px(),
                                            particles[i].py(),
                                            particles[i].pz(),
                                            particles[i].e());
            }
            
            if( dR2sq < R*R && dR2sq<dR1sq){
                
                teljet2_raw+=TLorentzVector(particles[i].px(),
                                            particles[i].py(),
                                            particles[i].pz(),
                                            particles[i].e());
            }
            
        }// end loop over pf-candidates
        
        if (fabs(teljet1_raw.Pt()) < 1e-9 || fabs(teljet2_raw.Pt()) < 1e-9){
            telHmass_histo->Fill( 0. );
            telHpt_histo->Fill( 0. );
            myTree.Hmass.push_back   ( 0.  );
            myTree.Hpt.push_back     ( 0. );
            myTree.R.push_back       (  R                     );
            myTree.jet1_pt.push_back (  0.          );
            myTree.jet1_eta.push_back(  -99.         );
            myTree.jet2_pt.push_back (  0.          );
            myTree.jet2_eta.push_back(  -99.         );
        }
        else{
            
            TLorentzVector teljet1(0.,0.,0.,0.);
            TLorentzVector teljet2(0.,0.,0.,0.);            
            
            if (isGEN_){
                teljet1.SetPxPyPzE(teljet1_raw.Px(),teljet1_raw.Py(),teljet1_raw.Pz(),teljet1_raw.E());
                teljet2.SetPxPyPzE(teljet2_raw.Px(),teljet2_raw.Py(),teljet2_raw.Pz(),teljet2_raw.E());                            
            }
            else{
                // Need to correct the Tel Jet! 
                double curArea = TMath::Pi()*R*R;
                double pTRaw_1 = teljet1_raw.Pt(); 
                double pTRaw_2 = teljet2_raw.Pt(); 
                double corrFactor1 = 1. - (rho*curArea/pTRaw_1);
                double corrFactor2 = 1. - (rho*curArea/pTRaw_2);
                teljet1.SetPxPyPzE(corrFactor1*teljet1_raw.Px(),corrFactor1*teljet1_raw.Py(),corrFactor1*teljet1_raw.Pz(),corrFactor1*teljet1_raw.E());
                teljet2.SetPxPyPzE(corrFactor2*teljet2_raw.Px(),corrFactor2*teljet2_raw.Py(),corrFactor2*teljet2_raw.Pz(),corrFactor2*teljet2_raw.E());                            
            }
            
            if (fabs(R - 0.5) < 0.001){
                std::cout << "r = " << R << "jet1_corr = " << teljet1.Pt() << ", jet2_corr = " << teljet2.Pt() << std::endl;
            }
            
            
            if( (teljet1+teljet2).M() > 110. &&
               (teljet1+teljet2).M() < 140. &&
               teljet1.Pt()          > 25.  &&
               teljet2.Pt()          > 25.  ){
                
                myTree.fracPass += 1./numR ;
                
            }
            
            telHmass_histo->Fill( (teljet1+teljet2).M() );
            telHpt_histo->Fill( (teljet1+teljet2).Pt() );
            
            myTree.Hmass.push_back   ( (teljet1+teljet2).M()  );
            myTree.Hpt.push_back     ( (teljet1+teljet2).Pt() );
            myTree.R.push_back       (  R                     );
            myTree.jet1_pt.push_back (  teljet1.Pt()          );
            myTree.jet1_eta.push_back(  teljet1.Eta()         );
            myTree.jet2_pt.push_back (  teljet2.Pt()          );
            myTree.jet2_eta.push_back(  teljet2.Eta()         );
        }
    }// end loop over jet cones
    
    //std::cout << "end loop over jet cones" << std::endl;
    
    
    myTree.Zmass    = bestMass    ;
    myTree.Zpt      = bestPt      ;
    myTree.lep1_pt  = bestLep1_pt ;
    myTree.lep2_pt  = bestLep2_pt ;
    myTree.lep1_eta = bestLep1_eta;
    myTree.lep2_eta = bestLep2_eta;
    
    ZHbbTree->Fill();    
    
    //std::cout << "ZHbbTree->Fill();" << std::endl;
    
}


// ------------ method called once each job just before starting event loop  ------------
void 

ZHbbTreeFiller::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
ZHbbTreeFiller::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
void 
ZHbbTreeFiller::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
ZHbbTreeFiller::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
ZHbbTreeFiller::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
ZHbbTreeFiller::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
ZHbbTreeFiller::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
    
    edm::ParameterSetDescription desc;
    desc.setUnknown();
    descriptions.addDefault(desc);
    
}

//define this as a plug-in
DEFINE_FWK_MODULE(ZHbbTreeFiller);
