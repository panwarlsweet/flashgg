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
#include "flashgg/DataFormats/interface/GenDiJet.h"
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
        bool overlapRemoval_;

    };

    GenDiPhotonDiJetProducer::GenDiPhotonDiJetProducer( const ParameterSet &iConfig ) :
        genPhotonToken_( consumes<View<flashgg::GenPhotonExtra> >( iConfig.getParameter<InputTag> ( "src" ) ) ),
        genJetToken_( consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "jets" ) ) ),
        metToken_( consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag> ( "flashggMets" ) ) ),   //// adding for regresion
        genJetNuToken_( consumes<View<reco::GenJet> >( iConfig.getParameter<InputTag> ( "flashggHiggsBJets" ) ) ), //// adding for regresion  
        overlapRemoval_(false)
    {
        if( iConfig.exists("overlapRemoval") ) { 
            overlapRemoval_ = iConfig.getParameter<bool>("overlapRemoval");
        }
        produces<vector<GenDiPhoton> >(); // one per diphoton, always in same order, vector is more efficient than map
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
        

        // std::unique_ptr<vector<GenDiPhoton> > diphotons( new vector<GenDiPhoton> );
        std::unique_ptr<vector<GenDiJet> > dijets( new vector<GenDiJet> );
        //for( size_t ii = 0 ; ii < photons->size() ; ++ii ) {
        //  auto pi = photons->ptrAt( ii );
        //  for( size_t jj = ii + 1 ; jj < photons->size() ; ++jj ) {
        //      auto pj = photons->ptrAt( jj );
                std::vector<edm::Ptr<reco::GenJet> > seljets;
                for( size_t ij = 0 ; ij < jets->size() ; ++ij ) {
                    auto jet = jets->ptrAt( ij );
                    // if( ! overlapRemoval_ || 
                    //  ( reco::deltaR(*jet,pi->cand()) > 0.3 
                    //    && reco::deltaR(*jet,pj->cand()) > 0.3 ) ) {
                        seljets.push_back(jet);
                        //}
                }
                if( seljets.size() >= 2 ) { 
                    auto jet0 = seljets[0], jet1 = seljets[1];

                    ///////// adding met and gentjets with nu for regression study
                    //if (genjetsnu->size() >= 2){
                    //   for( size_t p = 0 ; p < genjetsnu->size() ; ++p ) {
                    //  auto pn = genjetsnu->ptrAt(p);
                    //  for( size_t q = p + 1 ; q < genjetsnu->size() ; ++q ) {
                    //      auto qn = genjetsnu->ptrAt(q);
                            for( size_t kk = 0 ; kk < mets->size() ; ++kk ) {
                              auto met = mets->ptrAt( kk );
                              
                              dijets->push_back(GenDiJet(jet0,jet1,met,jet0,jet1));

                          }
                        }
                  }
                //   }
            //}
    //}
        
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

