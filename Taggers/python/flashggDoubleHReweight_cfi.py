import FWCore.ParameterSet.Config as cms

import flashgg.Systematics.settings as settings
year = settings.year
#default values first
weightsFile = "flashgg/Taggers/data/HHTagger/HHreweight_2016nodes_04022018.root "
if year == "2016":
    weightsFile = "flashgg/Taggers/data/HHTagger/HHreweight_2016nodes_04022018.root"
elif year == "2017":
    weightsFile = "flashgg/Taggers/data/HHTagger/HHreweight_2017nodes_04022018.root"

flashggDoubleHReweight = cms.EDProducer("FlashggDoubleHReweighter",
                                        GenParticleTag = cms.InputTag( "flashggPrunedGenParticles" ), # to compute MC-truth info
                                        targetNode = cms.int32(-1), 
                                        weightsFile = cms.FileInPath("%s"%weightsFile), 
                                        NCOEFFSA = cms.uint32(15), 
                                        A_13TeV_SM = cms.vdouble(2.09078, 10.1517, 0.282307, 0.101205, 1.33191, -8.51168, -1.37309, 2.82636, 1.45767, -4.91761, -0.675197, 1.86189, 0.321422, -0.836276, -0.568156),# coefficients of 15 operators for SM point
                                        benchmarks_map = cms.PSet( #values of kl,k5,c2,cg,c2g for 12 benchmarks as in arxic 1507.02245 v4.
                                             kl = cms.vdouble(7.5,  1.0,  1.0, -3.5, 1.0, 2.4, 5.0, 15.0, 1.0, 10.0, 2.4, 15.0 ), 
                                             kt = cms.vdouble(1.0,  1.0,  1.0,  1.5, 1.0, 1.0, 1.0, 1.0,  1.0, 1.5 , 1.0, 1.0 ), 
                                             c2 = cms.vdouble(-1.0, 0.5, -1.5, -3.0, 0.0, 0.0, 0.0, 0.0,  1.0, -1.0, 0.0, 1.0), 
                                             cg = cms.vdouble(0.0, -0.8,  0.0,  0.0, 0.8 ,0.2, 0.2,-1.0, -0.6, 0.0, 1.0, 0.0), 
                                             c2g = cms.vdouble(0.0, 0.6, -0.8,  0.0,-1.0,-0.2,-0.2, 1.0,  0.6, 0.0, -1.0, 0.0),
                                        )
                                        )

reweight_producer = "flashggDoubleHReweight"
reweight_names = ["benchmark0","benchmark1","benchmark2","benchmark3","benchmark4","benchmark5","benchmark6","benchmark7","benchmark8","benchmark9","benchmark10","benchmark11"]
