#ifndef FLASHgg_GenDiJet_h
#define FLASHgg_GenDiJet_h

#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Common/interface/Ptr.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/WeightedObject.h"


namespace flashgg {

    class GenDiJet : public reco::LeafCandidate, public WeightedObject
    {
    public:
        GenDiJet();
        GenDiJet(edm::Ptr<reco::GenJet>, edm::Ptr<reco::GenJet>, edm::Ptr<flashgg::Met>, edm::Ptr<reco::GenJet>, edm::Ptr<reco::GenJet>);
        ~GenDiJet();

        virtual GenDiJet *clone() const { return ( new GenDiJet( *this ) ); }
                      
        const reco::GenJet & leadingJet() const { return *leadingJet_; };
        const reco::GenJet & subLeadingJet() const { return *subLeadingJet_; }

        const reco::GenJet & leadGenJetNu() const { return *leadGenJetNu_; };            //// extra object for bregression study 
        const reco::GenJet & subleadGenJetNu() const { return *subleadGenJetNu_; }       //// extra object for bregression study 
       
        const flashgg::Met& Met() const { return Met_; }  //// extra object for bregression study 
        
        
        LorentzVector dijet() const;
        LorentzVector dijetNu() const; //// extra object for bregression study 
        
        void setTag(const std::string & tag) { tag_ = tag; }
        std::string tag() const { return tag_; }
        bool isTagged(const std::string & cmp) const { return tag_ == cmp; }
        
        void setCategoryNumber(int cat) { cat_ = cat; }
        int categoryNumber() const { return cat_; }
        operator int() const { return categoryNumber(); }
        
    private:
        void computeP4AndOrder();
        
        edm::Ptr<reco::GenJet> leadingJet_;
        edm::Ptr<reco::GenJet> subLeadingJet_;

        edm::Ptr<reco::GenJet> leadGenJetNu_;   //// extra object for bregression study 
        edm::Ptr<reco::GenJet> subleadGenJetNu_; //// extra object for bregression study 
        
        flashgg::Met Met_;  //// extra object for bregression study 
        int cat_;
        std::string tag_;

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

