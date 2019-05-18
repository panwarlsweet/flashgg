#include "flashgg/DataFormats/interface/GenDiJet.h"
#include "CommonTools/CandUtils/interface/AddFourMomenta.h"


using namespace flashgg;

GenDiJet::GenDiJet() {}

GenDiJet::~GenDiJet() {}

GenDiJet::GenDiJet( edm::Ptr<reco::GenJet> jet1, edm::Ptr<reco::GenJet> jet2 )
{
    leadingJet_ = *jet1;
    subLeadingJet_ = *jet2;
}

GenDiJet::GenDiJet( edm::Ptr<reco::GenJet> jet1, edm::Ptr<reco::GenJet> jet2, edm::Ptr<flashgg::Met> Met, edm::Ptr<reco::GenJet> jetNu1, edm::Ptr<reco::GenJet> jetNu2 )
    : GenDiJet(jet1,jet2)
{
    //// extra object for bregression study
    leadGenJetNu_ = jetNu1;
    subleadGenJetNu_ = jetNu2;
    Met_ = *Met;

    if( leadingJet_->pt() < subLeadingJet_->pt() ) { std::swap(leadingJet_,subLeadingJet_); }
    if( leadGenJetNu_->pt() < subleadGenJetNu_->pt() ) { std::swap(leadGenJetNu_,subleadGenJetNu_); }
}


GenDiJet::LorentzVector GenDiJet::dijet() const
{
    return leadingJet().p4() + subLeadingJet().p4();
}
//// extra object for bregression study
GenDiJet::LorentzVector GenDiJet::dijetNu() const
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

