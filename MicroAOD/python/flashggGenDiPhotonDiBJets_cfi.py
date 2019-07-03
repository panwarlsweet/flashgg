import FWCore.ParameterSet.Config as cms
from flashgg.Taggers.flashggTags_cff import UnpackedJetCollectionVInputTag

flashggGenDiPhotonDiBJets = cms.EDProducer('FlashggGenDiPhotonDiJetProducer',
                                     src=cms.InputTag('flashggGenPhotonsExtra'),
                                     jets=cms.InputTag('slimmedGenJets'),
                                     flashggMets=cms.InputTag('flashggMets'),
			             flashggGenNus=cms.InputTag('flashggGenNeutrinos'),
                                     flashggHiggsBJets=cms.InputTag('flashggHiggsBJets'),
                                     flashggDiPhotonTag=cms.InputTag('flashggPreselectedDiPhotons'),
                                     flashggJets=UnpackedJetCollectionVInputTag,
			             )

