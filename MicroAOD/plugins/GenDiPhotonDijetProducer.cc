#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/GenDiPhoton.h"
#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/GenPhotonExtra.h"

#include "TMVA/Reader.h"
#include "TMath.h"
#include "TVector3.h"
#include "TLorentzVector.h"

using namespace std;
using namespace edm;

namespace flashgg {

    class GenDiPhotonDiJetProducer : public EDProducer
    {

    public:
        GenDiPhotonDiJetProducer( const ParameterSet & );
    private:
        void produce( Event &, const EventSetup & ) override;

        EDGetTokenT<View<GenPhotonExtra> > genPhotonToken_;
        EDGetTokenT<View<reco::GenJet> > genJetToken_;
        EDGetTokenT<View<Met> > metToken_;
        EDGetTokenT<View<reco::GenJet> > genJetNuToken_;
        std::vector<edm::EDGetTokenT<View<flashgg::Jet> > >  recoJetToken_;
        EDGetTokenT<View<flashgg::DiPhotonCandidate> > diPhotonToken_;
        bool overlapRemoval_;

    };

    GenDiPhotonDiJetProducer::GenDiPhotonDiJetProducer( const ParameterSet &iConfig ) :
        genPhotonToken_( consumes<View<flashgg::GenPhotonExtra> >( iConfig.getParameter<InputTag> ( "src" ) ) ),
        genJetToken_( consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "jets" ) ) ),
        metToken_( consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag> ( "flashggMets" ) ) ),   //// adding for regresion
        genJetNuToken_( consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "flashggHiggsBJets" ) ) ), //// adding for regresion  
        //        recoJetToken_( consumes<View<flashgg::Jet> >( iConfig.getParameter<InputTag> ( "flashggJets" ) ) ),
	    diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "flashggDiPhotonTag" ) ) ),
        overlapRemoval_(false)
    {
        if( iConfig.exists("overlapRemoval") ) { 
            overlapRemoval_ = iConfig.getParameter<bool>("overlapRemoval");
        }
        produces<vector<GenDiPhoton> >(); // one per diphoton, always in same order, vector is more efficient than map
        auto jetTags = iConfig.getParameter<std::vector<edm::InputTag> > ( "flashggJets" );
        for( auto & tag : jetTags ) { recoJetToken_.push_back( consumes<edm::View<flashgg::Jet> >( tag ) ); }
    }
    
    void GenDiPhotonDiJetProducer::produce( Event &evt, const EventSetup & )
    {

        Handle<View<flashgg::GenPhotonExtra> > photons;
        evt.getByToken( genPhotonToken_, photons );
	
        Handle<View<reco::GenJet> > jets;
        evt.getByToken( genJetToken_, jets );

        Handle<View<flashgg::Met> > mets;       //// adding for regresion  
        evt.getByToken( metToken_, mets );      //// adding for regresion  

        Handle<View<reco::GenJet> > genjetsnu;       //// adding for regresion    
        evt.getByToken( genJetNuToken_, genjetsnu );   //// adding for regresion    

        Handle<View<flashgg::DiPhotonCandidate> > diPhos;
        evt.getByToken( diPhotonToken_, diPhos );

        double btag1=0.,  btag2=0.;// ref_bTag = -999., sum_bTag = 0.;
        //        int NGenJet=0, NB_GenJet=0, NRecoJet=0, Nneutrino=0, Nbquark=0;
        
        
        std::unique_ptr<vector<GenDiPhoton> > diphotons( new vector<GenDiPhoton> );
        for( size_t ii = 0 ; ii < photons->size() ; ++ii ) {
            auto pi = photons->ptrAt( ii );

            for( size_t jj = ii + 1 ; jj < photons->size() ; ++jj ) {
                auto pj = photons->ptrAt( jj );
                std::vector<edm::Ptr<reco::GenJet> > seljets;

                for( size_t ij = 0 ; ij < jets->size() ; ++ij ) {
                    auto jet = jets->ptrAt( ij );

                    if( ! overlapRemoval_ || 
                        ( reco::deltaR(*jet,pi->cand()) > 0.3 
                          && reco::deltaR(*jet,pj->cand()) > 0.3 ) ) {
                        seljets.push_back(jet);
                    }
                }

                if( seljets.size() >= 2 ) { 
                    auto jet0 = seljets[0], jet1 = seljets[1];

                    ///////// adding met and gentjets with nu for regression study
    
                    for( size_t p = 0 ; p < genjetsnu->size() ; ++p ) {
                      auto pn = genjetsnu->ptrAt(p);

                      for( size_t q = p + 1 ; q < genjetsnu->size() ; ++q ) {
                          auto qn = genjetsnu->ptrAt(q);

                            for( size_t kk = 0 ; kk < mets->size() ; ++kk ) {
                              auto met = mets->ptrAt( kk );
                              
                              //////// RECO level objects with basic analysis preselections

                              for( unsigned int candIndex = 0; candIndex < diPhos->size() ; candIndex++ ) {
                                  edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhos->ptrAt( candIndex );
                                  
                                  auto & leadPho = *(dipho->leadingPhoton());
                                  auto & subleadPho = *(dipho->subLeadingPhoton());

                                  double leadPt = leadPho.pt();
                                  double subleadPt = subleadPho.pt();
                                  
                                  if( leadPt <= 1./3. || subleadPt <= 0.25 ) { continue; }
                                  
                                  if(leadPho.passElectronVeto() < 1 || subleadPho.passElectronVeto() < 1){
                                      continue;
                                  }

                                  size_t vtx = (size_t)dipho->jetCollectionIndex();

                                  Handle<View<flashgg::Jet> > recojets;      
                                  evt.getByToken( recoJetToken_[vtx], recojets );

                                  std::vector<edm::Ptr<flashgg::Jet> > cleaned_recojets;

                                  for( size_t l = 0 ; l < recojets->size() ; ++l ) {
                                   auto lr = recojets->ptrAt(l);

				      if ( lr->pt() < 0. || fabs(lr->eta()) > 2.5 ) continue;
     
                                      btag1 = lr->bDiscriminator("pfDeepCSVJetTags:probb") + lr->bDiscriminator("pfDeepCSVJetTags:probbb"); // btaggable jets

                                      //  if ( btag1 < 0 || btag2 < 0 ) continue;
                                      if ( btag1 < 0 ) continue; 
                                      
                                      if( !lr->passesJetID  ( flashgg::Tight2017 ) ) continue;

                                      if( reco::deltaR( *lr, *(dipho->leadingPhoton()) ) > 0.4 && reco::deltaR( *lr, *(dipho->subLeadingPhoton()) ) > 0.4  ) {
                                          cleaned_recojets.push_back( lr );
                                      }
                                  }
                                  if( cleaned_recojets.size() < 2 ) { continue; }
                                  edm::Ptr<flashgg::Jet>  recoJet1, recoJet2;
                                  bool hasRecoDiJet = false;
                                  double  ref_bTag = -999.;

                                  for( size_t ijet=0; ijet < cleaned_recojets.size()-1;++ijet){
                                      auto jet_1 = cleaned_recojets[ijet];

                                      for( size_t kjet=ijet+1; kjet < cleaned_recojets.size();++kjet){
                                          auto jet_2 = cleaned_recojets[kjet];
                                          double sum_bTag=0.;

                                      //selecting DiJet with highest b-Tags score
                                          sum_bTag = jet_1->bDiscriminator("pfDeepCSVJetTags:probb") + jet_1->bDiscriminator("pfDeepCSVJetTags:probbb") + jet_2->bDiscriminator("pfDeepCSVJetTags:probb") + jet_2->bDiscriminator("pfDeepCSVJetTags:probbb");
                                      if (sum_bTag > ref_bTag){
                                          hasRecoDiJet = true;
                                          ref_bTag = sum_bTag;
                                          recoJet1 = jet_1;
                                          recoJet2 = jet_2;
                                      }
                                      //   diphotons->push_back(GenDiPhoton(pi,pj,jet0,jet1,met,pn,qn,lr,mr));
                                   }
                                  }
                                  if (!hasRecoDiJet)  continue; //selecting events with one DiJet object with highest bTag score from each event
        
                                  diphotons->push_back(GenDiPhoton(pi,pj,jet0,jet1,met,pn,qn,recoJet1,recoJet2,dipho)); 
                              }
                            }
                      }
                    }                           
                }
            }
        }
        
        evt.put( std::move( diphotons ) );
    }
}

typedef flashgg::GenDiPhotonDiJetProducer FlashggGenDiPhotonDiJetProducer;
DEFINE_FWK_MODULE( FlashggGenDiPhotonDiJetProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

