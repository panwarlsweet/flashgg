#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "flashgg/DataFormats/interface/DiPhotonCandidate.h"
#include "flashgg/DataFormats/interface/Jet.h"
#include "flashgg/DataFormats/interface/Met.h"
#include "flashgg/DataFormats/interface/Electron.h"
#include "flashgg/DataFormats/interface/Muon.h"

#include "flashgg/DataFormats/interface/DiPhotonMVAResult.h"
#include "flashgg/DataFormats/interface/DoubleHTag.h"
#include "flashgg/DataFormats/interface/TagTruthBase.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "flashgg/Taggers/interface/LeptonSelection.h"
#include "flashgg/MicroAOD/interface/MVAComputer.h"
#include "flashgg/DataFormats/interface/DoubleHttHTagger.h"

#ifdef CMSSW9
#include "PhysicsTools/TensorFlow/interface/TensorFlow.h"
#elif CMSSW8
#include "DNN/Tensorflow/interface/Graph.h"
#include "DNN/Tensorflow/interface/Tensor.h"
#endif

#include <vector>
#include <algorithm>
#include "TGraph.h"

using namespace std;
using namespace edm;

namespace flashgg {

    class DoubleHTagProducer : public EDProducer
    {

        public:
            typedef math::XYZPoint Point;

            DoubleHTagProducer( const ParameterSet & );
        private:
            void produce( Event &, const EventSetup & ) override;
            float EvaluateNN();
            int chooseCategory( float mva, float mx );
            bool isclose(double a, double b, double rel_tol, double abs_tol);        
            EDGetTokenT<View<DiPhotonCandidate> > diPhotonToken_;
            std::vector<edm::EDGetTokenT<edm::View<flashgg::Jet> > > jetTokens_;
            EDGetTokenT<View<reco::GenParticle> > genParticleToken_;
            string systLabel_;

            double minLeadPhoPt_, minSubleadPhoPt_;
            bool scalingPtCuts_, doPhotonId_, doMVAFlattening_, doCategorization_, dottHTagger_;
            double photonIDCut_;
            double vetoConeSize_;         
            unsigned int doSigmaMDecorr_;
            edm::FileInPath sigmaMDecorrFile_;
            std::vector<int> photonElectronVeto_;

            DecorrTransform* transfEBEB_;
            DecorrTransform* transfNotEBEB_;


            double minJetPt_;
            double maxJetEta_;
            vector<double>mjjBoundaries_;
            vector<double>mjjBoundariesLower_;
            vector<double>mjjBoundariesUpper_;
            std::string bTagType_;
            bool       useJetID_;
            string     JetIDLevel_;        

            flashgg::MVAComputer<DoubleHTag> mvaComputer_;
            vector<double> mvaBoundaries_, mxBoundaries_;
            int multiclassSignalIdx_;

            edm::FileInPath MVAFlatteningFileName_;
            TFile * MVAFlatteningFile_;
            TGraph * MVAFlatteningCumulative_;

            DoubleHttHTagger tthKiller_;
            float ttHTagger;
            edm::EDGetTokenT<edm::View<flashgg::Met> > METToken_;
            edm::EDGetTokenT<edm::View<flashgg::Electron> > electronToken_;
            edm::EDGetTokenT<edm::View<flashgg::Muon> > muonToken_;
            edm::EDGetTokenT<edm::View<reco::Vertex> > vertexToken_;
            edm::EDGetTokenT<double> rhoToken_;

            std::vector<double> HLF_VectorVar_;
            std::vector<std::vector<double>> PL_VectorVar_;
            string ttHWeightfileName_ ;
#ifdef CMSSW9
            tensorflow::Session* session_ttH;
#elif CMSSW8
            dnn::tf::Graph NNgraph_;
#endif
    };

    DoubleHTagProducer::DoubleHTagProducer( const ParameterSet &iConfig ) :
        diPhotonToken_( consumes<View<flashgg::DiPhotonCandidate> >( iConfig.getParameter<InputTag> ( "DiPhotonTag" ) ) ),
        genParticleToken_( consumes<View<reco::GenParticle> >( iConfig.getParameter<InputTag> ( "GenParticleTag" ) ) ),
        systLabel_( iConfig.getParameter<string> ( "SystLabel" ) ),
        minLeadPhoPt_( iConfig.getParameter<double> ( "MinLeadPhoPt" ) ),
        minSubleadPhoPt_( iConfig.getParameter<double> ( "MinSubleadPhoPt" ) ),
        scalingPtCuts_( iConfig.getParameter<bool> ( "ScalingPtCuts" ) ),
        vetoConeSize_( iConfig.getParameter<double> ( "VetoConeSize" ) ),
        minJetPt_( iConfig.getParameter<double> ( "MinJetPt" ) ),
        maxJetEta_( iConfig.getParameter<double> ( "MaxJetEta" ) ),
        bTagType_( iConfig.getUntrackedParameter<std::string>( "BTagType") ),
        useJetID_( iConfig.getParameter<bool>   ( "UseJetID"     ) ),
        JetIDLevel_( iConfig.getParameter<string> ( "JetIDLevel"   ) ),
        mvaComputer_( iConfig.getParameter<edm::ParameterSet>("MVAConfig") ),
        ttHWeightfileName_( iConfig.getUntrackedParameter<std::string>("ttHWeightfile"))
    {
        mjjBoundaries_ = iConfig.getParameter<vector<double > >( "MJJBoundaries" ); 
        mvaBoundaries_ = iConfig.getParameter<vector<double > >( "MVABoundaries" );
        mxBoundaries_ = iConfig.getParameter<vector<double > >( "MXBoundaries" );
        mjjBoundariesLower_ = iConfig.getParameter<vector<double > >( "MJJBoundariesLower" ); 
        mjjBoundariesUpper_ = iConfig.getParameter<vector<double > >( "MJJBoundariesUpper" ); 
        multiclassSignalIdx_ = (iConfig.getParameter<edm::ParameterSet>("MVAConfig")).getParameter<int>("multiclassSignalIdx"); 

        auto jetTags = iConfig.getParameter<std::vector<edm::InputTag> > ( "JetTags" ); 
        for( auto & tag : jetTags ) { jetTokens_.push_back( consumes<edm::View<flashgg::Jet> >( tag ) ); }

        assert(is_sorted(mvaBoundaries_.begin(), mvaBoundaries_.end()) && "mva boundaries are not in ascending order (we count on that for categorization)");
        assert(is_sorted(mxBoundaries_.begin(), mxBoundaries_.end()) && "mx boundaries are not in ascending order (we count on that for categorization)");
        doPhotonId_ = iConfig.getUntrackedParameter<bool>("ApplyEGMPhotonID");        
        photonIDCut_ = iConfig.getParameter<double>("PhotonIDCut");

        doMVAFlattening_ = iConfig.getParameter<bool>("doMVAFlattening"); 
        doCategorization_ = iConfig.getParameter<bool>("doCategorization"); 
        dottHTagger_ = iConfig.getParameter<bool>("dottHTagger"); 
        photonElectronVeto_=iConfig.getUntrackedParameter<std::vector<int > >("PhotonElectronVeto");
        //needed for HHbbgg MVA
        if(doMVAFlattening_){
            MVAFlatteningFileName_=iConfig.getUntrackedParameter<edm::FileInPath>("MVAFlatteningFileName");
            MVAFlatteningFile_ = new TFile((MVAFlatteningFileName_.fullPath()).c_str(),"READ");
            MVAFlatteningCumulative_ = (TGraph*)MVAFlatteningFile_->Get("cumulativeGraph"); 
        }

        doSigmaMDecorr_ = iConfig.getUntrackedParameter<unsigned int>("DoSigmaMDecorr");
        if(doSigmaMDecorr_){
            sigmaMDecorrFile_ = iConfig.getUntrackedParameter<edm::FileInPath>("SigmaMDecorrFile");
            TFile* f_decorr = new TFile((sigmaMDecorrFile_.fullPath()).c_str(), "READ");
            TH2D* h_decorrEBEB_ = (TH2D*)f_decorr->Get("hist_sigmaM_M_EBEB"); 
            TH2D* h_decorrNotEBEB_ = (TH2D*)f_decorr->Get("hist_sigmaM_M_notEBEB");

            if(h_decorrEBEB_ && h_decorrNotEBEB_){
                transfEBEB_ = new DecorrTransform(h_decorrEBEB_ , 125., 1, 0);
                transfNotEBEB_ = new DecorrTransform(h_decorrNotEBEB_ , 125., 1, 0);

            } else {
                throw cms::Exception( "Configuration" ) << "The file "<<sigmaMDecorrFile_.fullPath()<<" provided for sigmaM/M decorrelation does not contain the expected histograms."<<std::endl;
            }
        }

        if(dottHTagger_){
            tthKiller_.initializeSelectionThresholds();
            tthKiller_.initializeMVAVariables();
            tthKiller_.setupMVA(tthKiller_.ttHMVAWeights_,tthKiller_.ttHMVAVars_);//FIXME do configurable
        }

        HLF_VectorVar_.resize(9);  // High-level features. 9 at the moment
        PL_VectorVar_.resize(6);
        for (int i = 0; i < 6; i++)
            PL_VectorVar_[i].resize(7); // List of particles. 6 objects. Each object has 7 attributes.


        //needed for ttH MVA
        METToken_= consumes<View<flashgg::Met> >( iConfig.getParameter<InputTag> ("METTag") ) ;
        electronToken_ = consumes<edm::View<flashgg::Electron> >( iConfig.getParameter<edm::InputTag> ("ElectronTag") );
        muonToken_ = consumes<edm::View<flashgg::Muon> >( iConfig.getParameter<edm::InputTag>("MuonTag") );
        vertexToken_ = consumes<edm::View<reco::Vertex> >( iConfig.getParameter<edm::InputTag> ("VertexTag") );
        rhoToken_ = consumes<double>( iConfig.getParameter<edm::InputTag>( "rhoTag" ) );

        produces<vector<DoubleHTag>>();
        produces<vector<TagTruthBase>>();
    }

    int DoubleHTagProducer::chooseCategory( float mvavalue, float mxvalue)
    {
        //// should return 0 if mva above all the numbers, 1 if below the first, ..., boundaries.size()-N if below the Nth, ...
        //this is for mva, then you have mx
        if (!doCategorization_) {
            return 0;
        }
        int mvaCat=-1;
        for( int n = 0 ; n < ( int )mvaBoundaries_.size() ; n++ ) {
            if( ( double )mvavalue > mvaBoundaries_[mvaBoundaries_.size() - n - 1] ) {
                mvaCat = n;
                break;
            }
        }

        if (mvaCat==-1) return -1;// Does not pass, object will not be produced

        int mxCat=-1;
        for( int n = 0 ; n < ( int )mxBoundaries_.size() ; n++ ) {
            if( ( double )mxvalue > mxBoundaries_[mxBoundaries_.size() - n - 1] ) {
                mxCat = n;
                break;
            }
        }


        if (mxCat==-1) return -1;// Does not pass, object will not be produced

        int cat=-1;
        cat = mvaCat*mxBoundaries_.size()+mxCat;

        //the schema is like this: (different from HHbbgg_ETH)
        //            "cat0 := MXbin0 * MVAcat0",   #High MX, High MVA
        //            "cat1 := MXbin1 * MVAcat0",   #High but lower MX , High MVA
        //            "cat2 := MXbin2 * MVAcat0",
        //            "cat3 := MXbin0 * MVAcat1",
        //            "cat4 := MXbin1 * MVAcat1",
        // [.................................]

        return cat;

    }

    bool DoubleHTagProducer::isclose(double a, double b, double rel_tol=1e-09, double abs_tol=0.0)
    {
        return fabs(a-b) <= max(rel_tol * max(fabs(a), fabs(b)), abs_tol);
    }

    void DoubleHTagProducer::produce( Event &evt, const EventSetup & )
    {
        // read diphotons
        Handle<View<flashgg::DiPhotonCandidate> > diPhotons;
        evt.getByToken( diPhotonToken_, diPhotons );

        // prepare output
        std::unique_ptr<vector<DoubleHTag> > tags( new vector<DoubleHTag> );
        std::unique_ptr<vector<TagTruthBase> > truths( new vector<TagTruthBase> );
        edm::RefProd<vector<TagTruthBase> > rTagTruth = evt.getRefBeforePut<vector<TagTruthBase> >();

        // MC truth
        TagTruthBase truth_obj;
        if( ! evt.isRealData() ) {
            Handle<View<reco::GenParticle> > genParticles;
            evt.getByToken( genParticleToken_, genParticles );
            Point higgsVtx(0.,0.,0.);
            for( unsigned int genLoop = 0 ; genLoop < genParticles->size(); genLoop++ ) {
                int pdgid = genParticles->ptrAt( genLoop )->pdgId();
                if( pdgid == 25 || pdgid == 22 ) {
                    higgsVtx = genParticles->ptrAt( genLoop )->vertex();
                    break;
                }
            }
            truth_obj.setGenPV( higgsVtx );
            truths->push_back( truth_obj );
        }

        // loop over diphotons
        for( unsigned int candIndex = 0; candIndex < diPhotons->size() ; candIndex++ ) {
            edm::Ptr<flashgg::DiPhotonCandidate> dipho = diPhotons->ptrAt( candIndex );

            // kinematic cuts on diphotons
            auto & leadPho = *(dipho->leadingPhoton());
            auto & subleadPho = *(dipho->subLeadingPhoton());

            double leadPt = leadPho.pt();
            double subleadPt = subleadPho.pt();
            if( scalingPtCuts_ ) {
                leadPt /= dipho->mass();
                subleadPt /= dipho->mass();
            }
            if( leadPt <= minLeadPhoPt_ || subleadPt <= minSubleadPhoPt_ ) { continue; }
            //apply egm photon id with given working point
            if(doPhotonId_){
                if(leadPho.userFloat("EGMPhotonMVA")<photonIDCut_ || subleadPho.userFloat("EGMPhotonMVA")<photonIDCut_){
                    continue;
                }
            }
            //electron veto
            if(leadPho.passElectronVeto()<photonElectronVeto_[0] || subleadPho.passElectronVeto()<photonElectronVeto_[1]){
                continue;
            }


            // find vertex associated to diphoton object
            size_t vtx = (size_t)dipho->jetCollectionIndex();
            // size_t vtx = (size_t)dipho->vertexIndex();
            //  if( vtx >= jetTokens_.size() ) { vtx = 0; }
            // and read corresponding jet collection
            edm::Handle<edm::View<flashgg::Jet> > jets;
            evt.getByToken( jetTokens_[vtx], jets);

            // photon-jet cross-cleaning and pt/eta/btag/jetid cuts for jets
            std::vector<edm::Ptr<flashgg::Jet> > cleaned_jets;
            for( size_t ijet=0; ijet < jets->size(); ++ijet ) {//jets are ordered in pt
                auto jet = jets->ptrAt(ijet);
                if (jet->pt()<minJetPt_ || fabs(jet->eta())>maxJetEta_)continue;
                if (jet->bDiscriminator(bTagType_)<0) continue;//FIXME threshold might not be 0?
                if( useJetID_ ){
                    if( JetIDLevel_ == "Loose" && !jet->passesJetID  ( flashgg::Loose ) ) continue;
                    if( JetIDLevel_ == "Tight" && !jet->passesJetID  ( flashgg::Tight ) ) continue;
                }
                if( reco::deltaR( *jet, *(dipho->leadingPhoton()) ) > vetoConeSize_ && reco::deltaR( *jet, *(dipho->subLeadingPhoton()) ) > vetoConeSize_ ) {
                    cleaned_jets.push_back( jet );
                }
            }
            if( cleaned_jets.size() < 2 ) { continue; }
            //dijet pair selection. Do pair according to pt and choose the pair with highest b-tag
            double sumbtag_ref = -999;
            bool hasDijet = false;
            edm::Ptr<flashgg::Jet>  jet1, jet2;
            for( size_t ijet=0; ijet < cleaned_jets.size()-1;++ijet){
                auto jet_1 = cleaned_jets[ijet];
                for( size_t kjet=ijet+1; kjet < cleaned_jets.size();++kjet){
                    auto jet_2 = cleaned_jets[kjet];
                    auto dijet_mass = (jet_1->p4()+jet_2->p4()).mass(); 
                    if (dijet_mass<mjjBoundaries_[0] || dijet_mass>mjjBoundaries_[1]) continue;
                    double sumbtag = jet_1->bDiscriminator(bTagType_) + jet_2->bDiscriminator(bTagType_);
                    if (sumbtag > sumbtag_ref) {
                        hasDijet = true;
                        sumbtag_ref = sumbtag;
                        jet1 = jet_1;
                        jet2 = jet_2;
                    }
                }
            }
            if (!hasDijet) continue;             

            auto & leadJet = jet1; 
            auto & subleadJet = jet2; 

            // prepare tag object
            DoubleHTag tag_obj( dipho, leadJet, subleadJet );
            tag_obj.setDiPhotonIndex( candIndex );
            tag_obj.setSystLabel( systLabel_ );

            if (tag_obj.dijet().mass()<mjjBoundaries_[0] || tag_obj.dijet().mass()>mjjBoundaries_[1]) continue;

            // compute extra variables here
            tag_obj.setMX( tag_obj.p4().mass() - tag_obj.dijet().mass() - tag_obj.diPhoton()->mass() + 250. );

            if(doSigmaMDecorr_){
                tag_obj.setSigmaMDecorrTransf(transfEBEB_,transfNotEBEB_);
            }

            // eval MVA discriminant
            std::vector<float> mva_vector = mvaComputer_.predict_prob(tag_obj);
            double mva = mva_vector[multiclassSignalIdx_];
            if(doMVAFlattening_){
                mva = MVAFlatteningCumulative_->Eval(mva);
            }

            tag_obj.setMVA( mva );
            //tag_obj.setMVAprob( mva_vector );

            //tth Tagger
            if(dottHTagger_){
                float sumEt=0.,njets=0.;
                std::map<std::string, float> ttHVars;
                std::vector<flashgg::Jet> cleanedDR_jets;
                for( size_t ijet=0; ijet < cleaned_jets.size();++ijet){
                    auto jet = cleaned_jets[ijet];
                    if( reco::deltaR(*jet,*leadJet)< vetoConeSize_) continue;
                    if( reco::deltaR(*jet,*subleadJet)< vetoConeSize_) continue;
                    sumEt+=jet->p4().pt();
                    njets+=1;
                    cleanedDR_jets.push_back(*jet);
                }
                ttHVars["sumEt"] = sumEt;
                edm::Handle<View<flashgg::Met> > METs;
                evt.getByToken( METToken_, METs );
                if( METs->size() != 1 )
                { std::cout << "WARNING number of MET is not equal to 1" << std::endl; }
                Ptr<flashgg::Met> theMET = METs->ptrAt( 0 );
                auto p4MET=theMET->p4();
                ttHVars["MET"]=p4MET.pt();
                ttHVars["phiMET"]=p4MET.phi();

                ttHVars["dPhi1"] = reco::deltaPhi(p4MET.Phi(), leadJet->p4().phi());
                ttHVars["dPhi2"] = reco::deltaPhi(p4MET.Phi(), subleadJet->p4().phi());
                ttHVars["PhoJetMinDr"] = tag_obj.getPhoJetMinDr();
                ttHVars["njets"] = njets;
                std::vector<flashgg::Jet> DiJet;
                DiJet.push_back(tag_obj.leadJet());
                DiJet.push_back(tag_obj.subleadJet());
                std::vector<float> Xtt = tthKiller_.XttCalculation(cleanedDR_jets,DiJet);
                if(Xtt.size()>1){
                    ttHVars["Xtt0"] = Xtt[0];
                    ttHVars["Xtt1"] = Xtt[1];
                }else{
                    ttHVars["Xtt0"] = 1000;
                    ttHVars["Xtt1"] = 0;
                }
                Handle<View<flashgg::Electron> > theElectrons;
                evt.getByToken( electronToken_, theElectrons );

                Handle<View<reco::Vertex> > vertices;
                evt.getByToken( vertexToken_, vertices );
                edm::Handle<double>  rho;
                evt.getByToken(rhoToken_,rho);

                ttHVars["ptdipho"] = dipho->p4().pt();
                ttHVars["etadipho"] = dipho->p4().eta();
                ttHVars["phidipho"] = dipho->p4().phi();

                std::vector<edm::Ptr<flashgg::Electron> > selectedElectrons = selectStdAllElectrons( theElectrons->ptrs(), vertices->ptrs(), tthKiller_.looseLeptonPtThreshold, tthKiller_.elecEtaThresholds, tthKiller_.useElecMVARecipe, tthKiller_.useElecLooseId, *rho, evt.isRealData() );
                std::vector<edm::Ptr<flashgg::Electron> > tagElectrons = tthKiller_.filterElectrons( selectedElectrons, *tag_obj.diPhoton(), leadJet->p4(), subleadJet->p4(), tthKiller_.dRPhoLeptonThreshold, tthKiller_.dRJetLeptonThreshold);

                if (tagElectrons.size() > 0) 
                {
                    ttHVars["pte1"] = tagElectrons.at( 0 )->p4().pt();
                    ttHVars["etae1"] = tagElectrons.at( 0 )->p4().eta();
                    ttHVars["phie1"] = tagElectrons.at( 0 )->p4().phi();
                }
                else 
                {
                    ttHVars["pte1"] = 0.;
                    ttHVars["etae1"] = 0.;
                    ttHVars["phie1"] = 0.;
                }
                if (tagElectrons.size() > 1) 
                {
                    ttHVars["pte2"] = tagElectrons.at( 1 )->p4().pt();     
                    ttHVars["etae2"] = tagElectrons.at( 1 )->p4().eta();     
                    ttHVars["phie2"] = tagElectrons.at( 1 )->p4().phi();     
                }
                else 
                {
                    ttHVars["pte2"] = 0.;
                    ttHVars["etae2"] = 0.;
                    ttHVars["phie2"] = 0.;
                } 
                Handle<View<flashgg::Muon> > theMuons;
                evt.getByToken( muonToken_, theMuons );
                std::vector<edm::Ptr<flashgg::Muon> > selectedMuons = selectAllMuons( theMuons->ptrs(), vertices->ptrs(), tthKiller_.muEtaThreshold, tthKiller_.looseLeptonPtThreshold, tthKiller_.muPFIsoSumRelThreshold);
                std::vector<edm::Ptr<flashgg::Muon> > tagMuons = tthKiller_.filterMuons( selectedMuons, *tag_obj.diPhoton(), leadJet->p4(), subleadJet->p4(), tthKiller_.dRPhoLeptonThreshold, tthKiller_.dRJetLeptonThreshold);

                if (tagMuons.size() > 0) 
                {
                    ttHVars["ptmu1"] = tagMuons.at( 0 )->p4().pt();
                    ttHVars["etamu1"] = tagMuons.at( 0 )->p4().eta();
                    ttHVars["phimu1"] = tagMuons.at( 0 )->p4().phi();
                }
                else 
                {
                    ttHVars["ptmu1"] = 0.;
                    ttHVars["etamu1"] = 0.;
                    ttHVars["phimu1"] = 0.;
                }
                if (tagMuons.size() > 1) 
                {
                    ttHVars["ptmu2"] = tagMuons.at( 1 )->p4().pt();    
                    ttHVars["etamu2"] = tagMuons.at( 1 )->p4().eta();    
                    ttHVars["phimu2"] = tagMuons.at( 1 )->p4().phi();    
                }
                else 
                {
                    ttHVars["ptmu2"] = 0.;
                    ttHVars["etamu2"] = 0.;
                    ttHVars["phimu2"] = 0.;
                }

                ttHVars["fabs_CosThetaStar_CS"] = abs(tag_obj.getCosThetaStar_CS(6500));//FIXME don't do hardcoded
                ttHVars["fabs_CosTheta_bb"] = abs(tag_obj.CosThetaAngles()[1]);

                //FIXME: Add normalization

                //9 HLFs: 'sumEt','dPhi1','dPhi2','PhoJetMinDr','njets','Xtt0',
                //'Xtt1','fabs_CosThetaStar_CS','fabs_CosTheta_bb'
                HLF_VectorVar_[0] = ttHVars["sumET"];
                HLF_VectorVar_[1] = ttHVars["dPhi1"];
                HLF_VectorVar_[2] = ttHVars["dPhi2"];
                HLF_VectorVar_[3] = ttHVars["PhoJetMinDr"];
                HLF_VectorVar_[4] = ttHVars["njets"];
                HLF_VectorVar_[5] = ttHVars["Xtt0"];
                HLF_VectorVar_[6] = ttHVars["Xtt1"];
                HLF_VectorVar_[7] = ttHVars["fabs_CosThetaStar_CS"];
                HLF_VectorVar_[8] = ttHVars["fabs_CosTheta_bb"];

                // 6 objects: ele1, ele2, mu1, mu2, dipho, MET
                // Each object has 7 attributes: pt, eta, phi, isele, ismuon, isdipho, isMET
                //
                // 0: leading ele
                PL_VectorVar_[0][0] = ttHVars["pte1"];
                PL_VectorVar_[0][1] = ttHVars["etae1"];
                PL_VectorVar_[0][2] = ttHVars["phie1"];
                PL_VectorVar_[0][3] = (isclose(ttHVars["pte1"],0)) ? 0 : 1; // isEle
                PL_VectorVar_[0][4] = 0; // isMuon
                PL_VectorVar_[0][5] = 0; // isDiPho
                PL_VectorVar_[0][6] = 0; // isMET

                // 1: subleading ele
                PL_VectorVar_[1][0] = ttHVars["pte2"];
                PL_VectorVar_[1][1] = ttHVars["etae2"];
                PL_VectorVar_[1][2] = ttHVars["phie2"];
                PL_VectorVar_[1][3] = (isclose(ttHVars["pte2"],0)) ? 0 : 1; // isEle
                PL_VectorVar_[1][4] = 0; // isMuon
                PL_VectorVar_[1][5] = 0; // isDiPho
                PL_VectorVar_[1][6] = 0; // isMET

                // 2: leading muon
                PL_VectorVar_[2][0] = ttHVars["ptmu1"];
                PL_VectorVar_[2][1] = ttHVars["etamu1"];
                PL_VectorVar_[2][2] = ttHVars["phimu1"];
                PL_VectorVar_[2][3] = 0; // isEle
                PL_VectorVar_[2][4] = (isclose(ttHVars["ptmu1"],0)) ? 0 : 1; // isMuon
                PL_VectorVar_[2][5] = 0; // isDiPho
                PL_VectorVar_[2][6] = 0; // isMET

                // 3: subleading muon
                PL_VectorVar_[3][0] = ttHVars["ptmu2"];
                PL_VectorVar_[3][1] = ttHVars["etamu2"];
                PL_VectorVar_[3][2] = ttHVars["phimu2"];
                PL_VectorVar_[3][3] = 0; //isEle
                PL_VectorVar_[3][4] = (isclose(ttHVars["ptmu2"],0)) ? 0 : 1; // isMuon
                PL_VectorVar_[3][5] = 0; // isDiPho
                PL_VectorVar_[3][6] = 0; // isMET

                // 4: dipho
                PL_VectorVar_[4][0] = ttHVars["ptdipho"];
                PL_VectorVar_[4][1] = ttHVars["etadipho"];
                PL_VectorVar_[4][2] = ttHVars["phidipho"];
                PL_VectorVar_[4][3] = 0; // isEle
                PL_VectorVar_[4][4] = 0; // isMuon
                PL_VectorVar_[4][5] = (isclose(ttHVars["ptdipho"],0)) ? 0 : 1; // isDiPho
                PL_VectorVar_[4][6] = 0; // isMET

                // 5: MET
                PL_VectorVar_[5][0] = ttHVars["MET"];
                PL_VectorVar_[5][1] = 0; // MET eta
                PL_VectorVar_[5][2] = ttHVars["phiMET"];
                PL_VectorVar_[5][3] = 0; //isEle
                PL_VectorVar_[5][4] = 0; // isMuon
                PL_VectorVar_[5][5] = 0; // isDiPho
                PL_VectorVar_[5][6] = (isclose(ttHVars["MET"],0)) ? 0 : 1; // isMET

                // Sort by pT
                std::sort(PL_VectorVar_.rbegin(), PL_VectorVar_.rend()); 


#ifdef CMSSW9
                tensorflow::GraphDef* graphDef_ttH = tensorflow::loadGraphDef(ttHWeightfileName_.c_str());
                session_ttH = tensorflow::createSession(graphDef_ttH);

#elif CMSSW8
                NNgraph_ = *(new dnn::tf::Graph(ttHWeightfileName_.c_str()));
#endif 

                float ttHScore = EvaluateNN();

                tag_obj.ttHScore_ = ttHScore;

            }

            // choose category and propagate weights
            int catnum = chooseCategory( tag_obj.MVA(), tag_obj.MX() );
            tag_obj.setCategoryNumber( catnum );
            tag_obj.includeWeights( *dipho );
            tag_obj.includeWeights( *leadJet );
            tag_obj.includeWeights( *subleadJet );

            if (catnum>-1){
                if (doCategorization_) {
                    if (tag_obj.dijet().mass()<mjjBoundariesLower_[catnum] || tag_obj.dijet().mass()>mjjBoundariesUpper_[catnum]) continue;
                }
                tags->push_back( tag_obj );
                // link mc-truth
                if( ! evt.isRealData() ) {
                    tags->back().setTagTruth( edm::refToPtr( edm::Ref<vector<TagTruthBase> >( rTagTruth, 0 ) ) );                 
                }
            }
        }
        evt.put( std::move( truths ) );
        evt.put( std::move( tags ) );
    }

    
    float DoubleHTagProducer::EvaluateNN(){

        unsigned int shape = HLF_VectorVar_.size();
        unsigned int plshape1 = PL_VectorVar_.size();
        unsigned int plshape2 = PL_VectorVar_[0].size();
        std::printf("shape = %d %d\n", plshape1, plshape2); 
#ifdef CMSSW9
        tensorflow::Tensor HLFinput(tensorflow::DT_FLOAT, {1,shape});
        for (unsigned int i = 0; i < shape; i++){
            HLFinput.matrix<float>()(0,i) =  float(HLF_VectorVar_[i]);
        }
        tensorflow::Tensor PLinput(tensorflow::DT_FLOAT, tensorflow::TensorShape({1,plshape1, plshape2}));
        for (unsigned int i = 0; i < plshape1; i++)
            for (unsigned int j = 0; j < plshape2; j++)
            {
                PLinput.matrix<float>()(0, i, j) = float(PL_VectorVar_[i][j]);
            }
        std::vector<tensorflow::Tensor> outputs;
        tensorflow::run(session_ttH, { {"input_216:0", HLFinput}, {"input_215:0", PLinput} }, { "dense_302/Sigmoid:0" }, &outputs);
        float NNscore = outputs[0].matrix<float>()(0, 0);
#elif CMSSW8
        dnn::tf::Shape xShape[] = { 1, shape };
        dnn::tf::Tensor* hlf_input = DNNmodel_->defineInput(new dnn::tf::Tensor("input_216:0", 2, xShape));
        dnn::tf::Shape plShape[] = { 1, plshape1, plshape2 };
        dnn::tf::Tensor* pl_input = DNNmodel_->defineInput(new dnn::tf::Tensor("input_215:0", 3, plShape));

        dnn::tf::Tensor* y = DNNmodel_->defineOutput(new dnn::tf::Tensor("dense_302/Sigmoid:0"));
        for (unsigned int i = 0; i < hlf_input->getShape(1); i++){
            hlf_input->setValue<double>(0, i, HLF_VectorVar_[i]);
        }
        for (unsigned int i = 0; i < pl_input->getShape(1); i++)
            for (unsigned int j = 0; j < pl_input->getShape(2); j++)
                pl_input->setValue<double>(0, i, j, PL_VectorVar_[i][j]);
        DNNmodel_->eval();
        float NNscore = y->getValue<float>(0,0);
#endif         
        std::cout << "Output score is " << NNscore << std::endl;
        return NNscore;
    }//end EvaluateNN
    


}

typedef flashgg::DoubleHTagProducer FlashggDoubleHTagProducer;
DEFINE_FWK_MODULE( FlashggDoubleHTagProducer );
// Local Variables:
// mode:c++
// indent-tabs-mode:nil
// tab-width:4
// c-basic-offset:4
// End:
// vim: tabstop=4 expandtab shiftwidth=4 softtabstop=4
