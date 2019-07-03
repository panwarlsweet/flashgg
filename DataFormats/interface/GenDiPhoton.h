#ifndef FLASHgg_GenDiPhoton_h
#define FLASHgg_GenDiPhoton_h

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "flashgg/DataFormats/interface/GenPhotonExtra.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/WeightedObject.h"
#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"

namespace flashgg {

    class GenDiPhoton : public DiPhotonTagBase, public reco::LeafCandidate//, public WeightedObject
    {
    public:
        GenDiPhoton();
        GenDiPhoton( edm::Ptr<flashgg::GenPhotonExtra>, edm::Ptr<flashgg::GenPhotonExtra>);
        GenDiPhoton( edm::Ptr<flashgg::GenPhotonExtra>, edm::Ptr<flashgg::GenPhotonExtra>, edm::Ptr<reco::GenJet>, edm::Ptr<reco::GenJet>, std::vector<edm::Ptr<pat::PackedGenParticle> >, double, double, edm::Ptr<flashgg::Met>, edm::Ptr<reco::GenJet>, edm::Ptr<reco::GenJet>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::DiPhotonCandidate>);
        ~GenDiPhoton();

        virtual GenDiPhoton *clone() const { return ( new GenDiPhoton( *this ) ); }
        
        const flashgg::GenPhotonExtra::cand_type* leadingPhoton() const { return   &(leadingPhoton_.cand()); };
        const flashgg::GenPhotonExtra::cand_type* subLeadingPhoton() const { return &(subLeadingPhoton_.cand()); }
        
        const reco::GenJet & leadingJet() const { return *leadingJet_; };
        const reco::GenJet & subLeadingJet() const { return *subLeadingJet_; }

        const reco::GenJet & leadGenJetNu() const { return *leadGenJetNu_; };            //// extra object for bregression study 
        const reco::GenJet & subleadGenJetNu() const { return *subleadGenJetNu_; }       //// extra object for bregression study 

        const flashgg::GenPhotonExtra& leadingExtra() const { return leadingPhoton_; };
        const flashgg::GenPhotonExtra& subLeadingExtra() const { return subLeadingPhoton_; }
            
        const std::vector<edm::Ptr<pat::PackedGenParticle> >& Nus() const { return Nus_; }  //// extra object for bregression study

        const flashgg::Met& Met() const { return Met_; }  //// extra object for bregression study 
        
        const flashgg::Jet & recoJet1() const { return *leadJet_; } 
        const flashgg::Jet & recoJet2() const { return *subleadJet_; } 
        
        float sumPt() const { return  leadingPhoton_.cand().pt() + subLeadingPhoton_.cand().pt(); }
        bool operator <( const GenDiPhoton &b ) const { return ( sumPt() < b.sumPt() ); }
        bool operator >( const GenDiPhoton &b ) const { return ( sumPt() > b.sumPt() ); }

        float sumHt() const { return  sumPt() + leadingJet().pt() + subLeadingJet().pt(); }

        double sumNuPt() const{return sumNuPt_;}

        double sumNuPhi() const {return sumNuPhi_;}
        
        LorentzVector dijet() const;
        LorentzVector dijetNu() const; //// extra object for bregression study 
        const LorentzVector &reco_dijet() const { return reco_dijet_; }

        void setTag(const std::string & tag) { tag_ = tag; }
        std::string tag() const { return tag_; }
        bool isTagged(const std::string & cmp) const { return tag_ == cmp; }
        void setTagObj(const edm::Ptr<DiPhotonTagBase> recoTagObj) { recoTagObj_ = recoTagObj; }
        const edm::Ptr<DiPhotonTagBase> recoTagObj() const { return recoTagObj_; }
        
        
        void setCategoryNumber(int cat) { cat_ = cat; }
        int categoryNumber() const { return cat_; }
        operator int() const { return categoryNumber(); }
        
    private:
        void computeP4AndOrder();

        flashgg::GenPhotonExtra leadingPhoton_;
        flashgg::GenPhotonExtra subLeadingPhoton_;
        
        edm::Ptr<reco::GenJet> leadingJet_;
        edm::Ptr<reco::GenJet> subLeadingJet_;

        edm::Ptr<reco::GenJet> leadGenJetNu_;   //// extra object for bregression study 
        edm::Ptr<reco::GenJet> subleadGenJetNu_; //// extra object for bregression study 
        
        std::vector<edm::Ptr<pat::PackedGenParticle> > Nus_; //// extra object for bregression study
        double  sumNuPt_;
        double  sumNuPhi_;
        flashgg::Met Met_;  //// extra object for bregression study

        edm::Ptr<flashgg::Jet> leadJet_, subleadJet_;
        LorentzVector reco_dijet_;

        int cat_;
        std::string tag_;
        edm::Ptr<DiPhotonTagBase> recoTagObj_;
    };
    



}


#endif
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

