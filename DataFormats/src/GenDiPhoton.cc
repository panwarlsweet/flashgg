#include "flashgg/DataFormats/interface/GenDiPhoton.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"


using namespace flashgg;

GenDiPhoton::GenDiPhoton(): DiPhotonTagBase::DiPhotonTagBase() {}

GenDiPhoton::~GenDiPhoton() {}

GenDiPhoton::GenDiPhoton( edm::Ptr<flashgg::GenPhotonExtra> photon1, edm::Ptr<flashgg::GenPhotonExtra> photon2)
{
    leadingPhoton_ = *photon1;
    subLeadingPhoton_ = *photon2;
    computeP4AndOrder();
}

GenDiPhoton::GenDiPhoton( edm::Ptr<flashgg::GenPhotonExtra> photon1, edm::Ptr<flashgg::GenPhotonExtra> photon2, edm::Ptr<reco::GenJet> jet1, edm::Ptr<reco::GenJet> jet2, edm::Ptr<flashgg::Met> Met, edm::Ptr<reco::GenJet> jetNu1, edm::Ptr<reco::GenJet> jetNu2, edm::Ptr<flashgg::Jet> recoJet1, edm::Ptr<flashgg::Jet> recoJet2, edm::Ptr<flashgg::DiPhotonCandidate> diPho)
    : GenDiPhoton(photon1,photon2) 
{
    dipho_ = diPho;
    leadingJet_ = jet1;
    subLeadingJet_ = jet2;
    leadJet_ = recoJet1;
    subleadJet_ = recoJet2;
    //// extra object for bregression study
    leadGenJetNu_ = jetNu1;
    subleadGenJetNu_ = jetNu2;
    Met_ = *Met;
    reco_dijet_ = leadJet_->p4() + subleadJet_->p4();
    if( leadingJet_->pt() < subLeadingJet_->pt() ) { std::swap(leadingJet_,subLeadingJet_); }
    //if( leadGenJetNu_->pt() < subleadGenJetNu_->pt() ) { std::swap(leadGenJetNu_,subleadGenJetNu_); }
}

void GenDiPhoton::computeP4AndOrder()
{
    if( leadingPhoton()->pt() < subLeadingPhoton()->pt() ) {
        std::swap( leadingPhoton_, subLeadingPhoton_ );
    }
    this->setP4( leadingPhoton()->p4() + subLeadingPhoton()->p4() );
}

GenDiPhoton::LorentzVector GenDiPhoton::dijet() const
{
    return leadingJet().p4() + subLeadingJet().p4();
}
//// extra object for bregression study
GenDiPhoton::LorentzVector GenDiPhoton::dijetNu() const
{
    return leadGenJetNu().p4() + subleadGenJetNu().p4();
}


// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4

