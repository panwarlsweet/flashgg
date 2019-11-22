#ifndef flashgg_DoubleHTag
#define flashgg_DoubleHTag

#include "TLorentzVector.h"

#include "flashgg/DataFormats/interface/DiPhotonTagBase.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "DataFormats/Candidate/interface/LeafCandidate.h"
#include "DataFormats/Math/interface/deltaR.h"

#include "flashgg/Taggers/interface/FunctionHelpers.h"

namespace flashgg {

    class DoubleHTag: public DiPhotonTagBase, public reco::LeafCandidate
    {
    public:
        DoubleHTag();
        ~DoubleHTag();

        DoubleHTag( edm::Ptr<DiPhotonCandidate>, edm::Ptr<flashgg::Jet>, edm::Ptr<flashgg::Jet> );
        virtual DoubleHTag *clone() const override;
        /// DiPhotonTagBase::tag_t tagEnum() const override { return DiPhotonTagBase::kDoubleH; }

        void setMVA(double x) { mva_ = x; }
        double MVA() const { return mva_; }
      //  void setMVAprob(std::vector<float> &x) const { mva_prob_ = x; }
      //  std::vector<float> MVAprob() const { return mva_prob_; }
        void setMX(double x) { MX_ = x; }
        double MX() const { return MX_; }
        void setGenMhh(double x) { genMhh_ = x; }
        double genMhh() const { return genMhh_; }
        float ttHScore() const { return ttHScore_; }
        double diphotonPtOverM() const {return diPhoton()->pt()/mass(); }
        double dijetPtOverM() const {return dijet().pt()/mass(); }

        const flashgg::Jet & leadJet() const { return *leadJet_; } 
        const flashgg::Jet & subleadJet() const { return *subleadJet_; } 
        
        const LorentzVector & dijet() const { return dijet_; }

        float getCosThetaStar_CS() const;
        float getCosThetaStar_CS_old(float ebeam) const;
        std::vector<float> CosThetaAngles() const;
        float HelicityCosTheta( TLorentzVector Booster, TLorentzVector Boosted) const;
        float getPhoJetMinDr() const;
        float getSigmaMDecorr() const;
        float getSigmaMOverMJets() const;
        void  setSigmaMDecorrTransf( DecorrTransform* transfEBEB, DecorrTransform* transfNotEBEB){ transfEBEB_= transfEBEB; transfNotEBEB_=transfNotEBEB;}
        LorentzVector getdiHiggsP4() const {return p4();}
        void setBenchmarkReweight(std::vector<float> x) { benchmark_reweights_ = x; }
        float getBenchmarkReweight(int targetNode) const { return benchmark_reweights_[targetNode]; }
        void setEventNumber(double x) { eventNumber_ = x; }
        double eventNumber() const { return eventNumber_; }

        float ttHScore_;
        float sumET_, MET_, phiMET_, dPhi1_, dPhi2_, PhoJetMinDr_, njets_, Xtt0_, Xtt1_, pte1_, pte2_, ptmu1_, ptmu2_, ptdipho_, etae1_, etae2_, etamu1_, etamu2_, etadipho_, phie1_, phie2_, phimu1_, phimu2_, phidipho_, fabs_CosThetaStar_CS_, fabs_CosTheta_bb_, mjj_, ptjet1_, ptjet2_, etajet1_, etajet2_, phijet1_, phijet2_, ttDecay_ID_, deepjetCdiscr_jet1_, deepjetCdiscr_jet2_, deepjetCvsLdiscr_jet1_, deepjetCvsLdiscr_jet2_, deepcsvCdiscr_jet1_, deepcsvCdiscr_jet2_, deepcsvCvsLdiscr_jet1_, deepcsvCvsLdiscr_jet2_,  MT_leadpho_met_, MT_subleadpho_met_, MT_dipho_met_, sumPT_Had_Act_, bgenJetNoNu_1_pt_, bgenJetNoNu_2_pt_, bgenJetNu_1_pt_, bgenJetNu_2_pt_, mbbNu_, mbbNoNu_; 
        float sumET() const {return sumET_;}
        float MET() const {return MET_;}
        float phiMET() const {return phiMET_;}
        float dPhi1() const {return dPhi1_;}
        float dPhi2() const {return dPhi2_;}
        float PhoJetMinDr() const {return PhoJetMinDr_;}
        float njets() const {return njets_;}
        float mjj() const {return mjj_;}
        float Xtt0() const {return Xtt0_;}
        float Xtt1() const {return Xtt1_;}
        float pte1() const {return pte1_;}
        float pte2() const {return pte2_;}
        float ptmu1() const {return ptmu1_;}
        float ptmu2() const {return ptmu2_;}
        float ptdipho() const {return ptdipho_;}
        float ptjet1() const {return ptjet1_;}
        float ptjet2() const {return ptjet2_;}
        float etajet1() const {return etajet1_;}
        float etajet2() const {return etajet2_;}
        float phijet1() const {return phijet1_;}
        float phijet2() const {return phijet2_;}
        float etae1() const {return etae1_;}
        float etae2() const {return etae2_;}
        float etamu1() const {return etamu1_;}
        float etamu2() const {return etamu2_;}
        float etadipho() const {return etadipho_;}
        float phie1() const {return phie1_;}
        float phie2() const {return phie2_;}
        float phimu1() const {return phimu1_;}
        float phimu2() const {return phimu2_;}
        float phidipho() const {return phidipho_;}
        float fabs_CosThetaStar_CS() const {return fabs_CosThetaStar_CS_;}
        float fabs_CosTheta_bb() const {return fabs_CosTheta_bb_;}
        int ttDecay_ID() const {return ttDecay_ID_;}
        float deepjetCdiscr_jet1() const {return deepjetCdiscr_jet1_;}
        float deepjetCdiscr_jet2() const {return deepjetCdiscr_jet2_;}
        float deepjetCvsLdiscr_jet1() const {return deepjetCvsLdiscr_jet1_;}
        float deepjetCvsLdiscr_jet2() const {return deepjetCvsLdiscr_jet2_;}
        float deepcsvCdiscr_jet1() const {return deepcsvCdiscr_jet1_;}
        float deepcsvCdiscr_jet2() const {return deepcsvCdiscr_jet2_;}
        float deepcsvCvsLdiscr_jet1() const {return deepcsvCvsLdiscr_jet1_;}
        float deepcsvCvsLdiscr_jet2() const {return deepcsvCvsLdiscr_jet2_;}
        float MT_leadpho_met() const {return MT_leadpho_met_;}
        float MT_subleadpho_met() const {return MT_subleadpho_met_;}
        float MT_dipho_met() const {return MT_dipho_met_;}
        float sumPT_Had_Act() const {return sumPT_Had_Act_;}
        float reg_mbbNu() const {return mbbNu_;}
        float reg_mbbNoNu() const {return mbbNoNu_;}
        float reg_bgenJetNoNu_1_pt() const {return bgenJetNoNu_1_pt_;}
        float reg_bgenJetNoNu_2_pt() const {return bgenJetNoNu_2_pt_;}
        float reg_bgenJetNu_1_pt() const {return bgenJetNu_1_pt_;}
        float reg_bgenJetNu_2_pt() const {return bgenJetNu_2_pt_;}

    private:
        double mva_, MX_, genMhh_;
        vector<float> benchmark_reweights_;
 //       std::vector<float> mva_prob_;
         long eventNumber_;
        edm::Ptr<flashgg::Jet> leadJet_, subleadJet_;
        LorentzVector dijet_;
        DecorrTransform* transfEBEB_;
        DecorrTransform* transfNotEBEB_;
        
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

