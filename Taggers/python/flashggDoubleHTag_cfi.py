import os
import FWCore.ParameterSet.Config as cms

import flashgg.Taggers.dumperConfigTools as cfgTools

from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag
import flashgg.Taggers.PUJID_wps as pujid

from flashgg.Taggers.globalVariables_cff import globalVariables

import flashgg.Systematics.settings as settings
year = settings.year
#default values first
year_norm = 0
jetPUID = 'Loose'
weightsFile="flashgg/Taggers/data/HHTagger/training_with_10_12_2018_commonTraining_2016.weights.xml"# path to TMVA weights
MVAscalingValue=cms.double(1.)#scale MVA output before the cumulative transformation for 2017(2016 kept unchanged for simplicity, we will probably change that once we have all 3 years.)

if year == "2016":
    year_norm = 0
    jetPUID = 'Loose'
    weightsFile="flashgg/Taggers/data/HHTagger/training_with_10_12_2018_commonTraining_2016.weights.xml", 
    MVAscalingValue=1.
elif year == "2017":
    year_norm = 1
    jetPUID = 'Tight2017'
    weightsFile="flashgg/Taggers/data/HHTagger/training_with_10_12_2018_commonTraining_2017.weights.xml", 
    MVAscalingValue=1.011026



flashggDoubleHTag = cms.EDProducer("FlashggDoubleHTagProducer",
                                   DiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'), # diphoton collection (will be replaced by systematics machinery at run time)
                                   JetTags= UnpackedJetCollectionVInputTag, # one jet per vertex
                                   GenParticleTag = cms.InputTag( "flashggPrunedGenParticles" ), # to compute MC-truth info
                                   SystLabel      = cms.string(""), # used by systematics machinery
                                   
                                   VetoConeSize   = cms.double(0.4),
                                   MinLeadPhoPt   = cms.double(1./3.),
                                   MinSubleadPhoPt   = cms.double(0.25),
                                   ScalingPtCuts = cms.bool(True),
                                   DoSigmaMDecorr =cms.untracked.uint32(1),#transformation of sigmaM/M
                                   SigmaMDecorrFile = cms.untracked.FileInPath("flashgg/Taggers/data/diphoMVA_sigmaMoMdecorr_split_Mgg40_180.root"),
                                   ApplyEGMPhotonID = cms.untracked.bool(False),
                                   PhotonIDCut = cms.double(0.2),#this is loose id for 2016
                                   PhotonElectronVeto =cms.untracked.vint32(1, 1), #0: Pho1, 1: Pho2

                                   MinJetPt   = cms.double(25.),
                                   MaxJetEta   = cms.double(2.5),
                                   MJJBoundaries = cms.vdouble(70.,190.),
                                  # BTagType = cms.untracked.string('pfCombinedInclusiveSecondaryVertexV2BJetTags'), #string for btag algorithm
                                 #c  BTagType = cms.untracked.string('pfDeepCSVJetTags:probb'), #string for btag algorithm
                                   BTagType = cms.vstring('pfDeepCSVJetTags:probb','pfDeepCSVJetTags:probbb'), #string for btag algorithm
                                   UseJetID = cms.bool(True),
                                   JetIDLevel = cms.string(jetPUID),

                                   MVABoundaries  = cms.vdouble(0.29,0.441, 0.724), # category boundaries for MVA
                                   MXBoundaries   = cms.vdouble(250., 354., 478., 560.), # .. and MX
                                   MJJBoundariesLower = cms.vdouble(98.0,95.0,97.0,96.0,95.0,95.0,95.0,95.0,95.0,95.0,95.0,95.0),#for each category following the convention cat0=MX0 MVA0, cat1=MX1 MVA0, cat2=MX2 MVA0....
                                   MJJBoundariesUpper = cms.vdouble(150.0,150.0,143.0,150.0,150.0,150.0,150.0,145.0,155.0,142.0,146.0,152.0),#for each category following the convention cat0=MX0 MVA0, cat1=MX1 MVA0, cat2=MX2 MVA0....
                                   MVAConfig = cms.PSet(variables=cms.VPSet(), # variables are added below
                                                        classifier=cms.string("BDT::bdt"), # classifier name
                                                        weights=cms.FileInPath("%s"%weightsFile), 
                                                        regression=cms.bool(False), # this is not a regression
                                                        multiclass=cms.bool(True), # this is multiclass 
                                                        multiclassSignalIdx=cms.int32(2), # this is multiclass index for Signal
                                                        ),

                                   doMVAFlattening=cms.bool(True),#do transformation of cumulative to make it flat
                                   MVAscaling=cms.double(MVAscalingValue),
                                   doCategorization=cms.bool(False),#do categorization based on MVA x MX or only fill first tree with all events
                                   MVAFlatteningFileName=cms.untracked.FileInPath("flashgg/Taggers/data/HHTagger/cumulativeTransformation_20181210_common_2016_2017.root"),#FIXME, this should be optional, is it?
                                   globalVariables=globalVariables,
                                   
                                   dottHTagger=cms.bool(True), #whether to do ttH killer. 
                                    
                                   ElectronTag=cms.InputTag('flashggSelectedElectrons'),
                                   MuonTag=cms.InputTag('flashggSelectedMuons'),
                                   VertexTag=cms.InputTag('offlineSlimmedPrimaryVertices'),
                                   METTag=cms.InputTag('flashggMets'),
                                   rhoTag = cms.InputTag('fixedGridRhoFastjetAll'),
                                   looseLeptonPtThreshold = cms.double(10.),
                                   muonEtaThreshold = cms.double(2.4),
                                   muPFIsoSumRelThreshold = cms.double(0.25),
                                   deltaRPhoElectronThreshold = cms.double(1.),
                                   deltaRPhoMuonThreshold = cms.double(0.5),
                                   deltaRJetLepThreshold = cms.double(0.4),
                                   useElectronMVARecipe = cms.bool(False),
                                   useElectronLooseID = cms.bool(True),
                                   electronEtaThresholds=cms.vdouble(1.4442,1.566,2.5),
                                   ttHWeightfile2016 = cms.untracked.string(os.environ["CMSSW_BASE"]+"/src/flashgg/Taggers/data/ttHKiller/2016model.pb"),
                                   ttHWeightfile2017 = cms.untracked.string(os.environ["CMSSW_BASE"]+"/src/flashgg/Taggers/data/ttHKiller/2017model.pb"), 

                                   ttHScoreThreshold2016 = cms.double(0.0), #to be updated
                                   ttHScoreThreshold2017 = cms.double(0.0), #to be updated

                                   # For standardization
                                   mean2016 = cms.vdouble(2.97010674e+02,  6.73602509e+01,  1.53733458e-02,  2.91760738e-03,
                                            -4.50615632e-03,  1.28529658e+00,  5.28815498e+00,  9.27775937e+00,
                                              1.58128232e+02,  6.66740444e+01,  3.38588398e+01,  6.54014733e+01,
                                                3.69811443e+01,  1.32720750e+02,  5.12574639e-02,  2.62416583e-01,
                                                  8.39057188e-03,  6.40663447e-02,  6.53202679e-04,  8.78594807e-02,
                                                   -2.89199316e-01,  4.77685875e-02, -3.23971631e-01, -2.53331744e-03,
                                                     5.10083859e-01,  4.94531853e-01),
                                   std2016 = cms.vdouble(2.21478393e+02, 5.73321613e+01, 1.89617878e+00, 1.85361320e+00,
                                            1.82010584e+00, 5.88030990e-01, 1.72380424e+00, 2.07288648e+01,
                                             3.42620549e+02, 5.16832362e+01, 2.18013696e+01, 5.02077374e+01,
                                              2.27158339e+01, 1.01620307e+02, 1.10418092e+00, 1.26526355e+00,
                                               1.08953531e+00, 1.20423568e+00, 1.34930679e+00, 1.82733989e+00,
                                                1.97811838e+00, 1.84003858e+00, 1.78236270e+00, 1.79919968e+00,
                                                 2.89647229e-01, 2.85520706e-01),

                                   mean2017 = cms.vdouble(2.85986095e+02,  7.13652588e+01,  1.50018609e-01, -1.15857046e-03,
                                             7.31171337e-03,  1.27996563e+00,  5.00714430e+00,  1.08883721e+01,
                                               1.98869805e+02,  6.70324122e+01,  3.64951560e+01,  6.55364012e+01,
                                                 3.73749351e+01,  1.32349654e+02,  3.10983225e-03, -3.24629632e-02,
                                                  -1.50654180e-02, -3.20541466e-02, -1.89750793e-03, -2.03327023e-02,
                                                   -2.10177238e-01,  3.62470208e-02,  8.61734016e-02,  5.79434733e-03,
                                                     5.15481843e-01,  4.92488659e-01),

                                   std2017 = cms.vdouble(2.08976832e+02, 5.80302348e+01, 1.90812647e+00, 1.85441428e+00,
                                            1.82586757e+00, 5.87443460e-01, 1.61323090e+00, 2.31589605e+01,
                                             3.77645905e+02, 5.21242994e+01, 1.98463133e+01, 5.07578154e+01,
                                              2.33425594e+01, 1.02008742e+02, 1.11125066e+00, 1.16451668e+00,
                                               1.08169242e+00, 1.10331921e+00, 1.36011930e+00, 1.79742985e+00,
                                                1.84847216e+00, 1.82114239e+00, 1.75615738e+00, 1.81673458e+00,
                                                 2.90625874e-01, 2.85115622e-01)
                                  ) 

cfgTools.addVariables(flashggDoubleHTag.MVAConfig.variables,
                      # here the syntax is VarNameInTMVA := expression
                      [#"leadingJet_bDis := leadJet().bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",#FIXME make the btag type configurable?
                       #"subleadingJet_bDis := subleadJet().bDiscriminator('pfCombinedInclusiveSecondaryVertexV2BJetTags')",
                       "leadingJet_DeepCSV := leadJet().bDiscriminator('pfDeepCSVJetTags:probb')+leadJet().bDiscriminator('pfDeepCSVJetTags:probbb')",#FIXME make the btag type configurable?
                       "subleadingJet_DeepCSV := subleadJet().bDiscriminator('pfDeepCSVJetTags:probb')+subleadJet().bDiscriminator('pfDeepCSVJetTags:probbb')",
                       "absCosThetaStar_CS := abs(getCosThetaStar_CS(6500))",#FIXME get energy from somewhere?
                       "absCosTheta_bb := abs(CosThetaAngles()[1])",
                       "absCosTheta_gg := abs(CosThetaAngles()[0])",
                       "diphotonCandidatePtOverdiHiggsM := diphotonPtOverM()",
                       "dijetCandidatePtOverdiHiggsM := dijetPtOverM()",
                       "customLeadingPhotonIDMVA := diPhoton.leadingView.phoIdMvaWrtChosenVtx",
                       "customSubLeadingPhotonIDMVA := diPhoton.subLeadingView.phoIdMvaWrtChosenVtx",
                       "leadingPhotonSigOverE := diPhoton.leadingPhoton.sigEOverE",
                       "subleadingPhotonSigOverE := diPhoton.subLeadingPhoton.sigEOverE",
                     #  "sigmaMOverMDecorr := getSigmaMDecorr()",
                       "sigmaMOverM := sqrt(0.5*(diPhoton.leadingPhoton.sigEOverE*diPhoton.leadingPhoton.sigEOverE + diPhoton.subLeadingPhoton.sigEOverE*diPhoton.subLeadingPhoton.sigEOverE))",
                       "PhoJetMinDr := getPhoJetMinDr()",
                       "rho := global.rho",
                       "(leadingJet_bRegNNResolution*1.4826) := leadJet().userFloat('bRegNNResolution')*1.4826",
                       "(subleadingJet_bRegNNResolution*1.4826) := subleadJet().userFloat('bRegNNResolution')*1.4826",
                       "(sigmaMJets*1.4826) := getSigmaMOverMJets()*1.4826",
                     #  "year := %d"%year_norm,
                       ]
                      )


