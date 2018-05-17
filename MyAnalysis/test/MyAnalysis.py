#!/usr/bin/env python

from __future__ import print_function
import sys, os, shutil, re, subprocess
import ROOT
from DataFormats.FWLite import Events, Handle
import numpy as np

def main():

  from FWCore.ParameterSet.VarParsing import VarParsing
  
  options = VarParsing ('analysis')
  options.inputFiles  = ["/afs/cern.ch/work/l/lata/HH_bbyy/CMSSW_9_4_2/src/flashgg/root_files/myMicroAODOutput_6.root"]
  options.maxEvents = -1
  options.outputFile = "fout.root"
  options.parseArguments()

###### save hisots in root file ######

  fout = ROOT.TFile("output_3.root", 'RECREATE')
  fout.cd()

##### define Histograms ######

  h_H_pt			 = ROOT.TH1F("gen_H_pT" ,";Higgs p_{T} [GeV];Events;;", 50, 0., 500.)
  h_H_eta                        = ROOT.TH1F("gen_H_eta",";Higgs #eta;Events;;", 50, -5., 5.)
  h_b_pt                         = ROOT.TH1F("gen_b_pT" ,";p_{T}of b [GeV];Events;;", 50, 0., 500.)
  h_b_eta                        = ROOT.TH1F("gen_b_eta",";#eta of b;Events;;", 50, -5., 5.)
  h_y_pt                         = ROOT.TH1F("gen_y_pT" ,";p_{T}of #gamma[GeV];Events;;", 50, 0., 500)
  h_y_eta                        = ROOT.TH1F("gen_y_eta",";#eta of #gamma;Events;;", 50, -5., 5.)

  h_jet_pt 	                 = ROOT.TH1F("jet_pt",";jet p_{T} [GeV];Events;;", 50, 0., 500.)
  h_jet_eta                      = ROOT.TH1F("jet_eta",";jet #eta;Events;;", 50,-5., 5.)
  h_jet_deepcsv                  = ROOT.TH1F("jet_deepcsv",";jet DeepCSV [GeV];Events;;", 50, 0., 1.)
  h_btaggedjet0_pt               = ROOT.TH1F("btagged_jet0_pt",";b-tagged leading jet p_{T} [GeV];Events;;", 50, 0., 500.)
  h_btaggedjet0_eta              = ROOT.TH1F("btagged_jet0_eta",";b-tagged leading jet #eta;Events;;", 50,-5., 5.)
  h_btaggedjet1_pt               = ROOT.TH1F("btagged_jet1_pt",";b-tagged subleading jet p_{T} [GeV];Events;;",50, 0, 500)
  h_btaggedjet1_eta              = ROOT.TH1F("btagged_jet1_eta",";b-tagged subleading jet #eta;Events;;", 50,-5, 5)

  h_photon_pt			 = ROOT.TH1F("photon_pt",";#gamma p_{T} [GeV];Events;;", 50, 0., 500.)
  h_photon_eta                   = ROOT.TH1F("photon_eta",";#gamma #eta [GeV];Events;;", 50, -5, 5)
  h_photon0_pt                   = ROOT.TH1F("photon0_pt",";leading #gamma p_{T} [GeV];Events;;", 50, 0., 500.)
  h_photon0_eta                  = ROOT.TH1F("photon0_eta",";leading #gamma #eta [GeV];Events;;", 50, -5, 5)
  h_photon1_pt                   = ROOT.TH1F("photon1_pt",";subleading #gamma p_{T} [GeV];Events;;", 50, 0., 500)
  h_photon1_eta                  = ROOT.TH1F("photon1_eta",";subleading #gamma #eta [GeV];Events;;", 50, -5, 5)
  h_njet			 = ROOT.TH1D("njet",";Njet;Events;;",15, 0 ,15)
  h_nphoton     		 = ROOT.TH1D("nphoton",";Nphoton;Events;;",5, 0 ,5)

  h_dipho_mass			 = ROOT.TH1F("dipho_mass ",";diphoton mass [GeV];Events;;", 50, 0., 200.)
  h_reconsM_photon		 = ROOT.TH1F("photon_recomass",";M_{#gamma#gamma} [GeV];Events;;", 50, 0., 200.)
  h_reconsM_bb			 = ROOT.TH1F("bb_recomass",";M_{bb} [GeV];Events;;", 50, 0., 200)

  events = Events(options.inputFiles)
  nevents = 0
  for event in events:

    if options.maxEvents > 0 and nevents > options.maxEvents: break

    ### Find gen particles
    h_prunedGenpar = Handle("std::vector<reco::GenParticle>")
    event.getByLabel("flashggPrunedGenParticles", h_prunedGenpar)

    particles = h_prunedGenpar.product()
    for p in particles:
      print('prunedgenparticles P pid = {} status = {}'.format(p.pdgId(), p.status()))
      if p.pdgId() == 25:
        h_H_pt.Fill(p.pt())
        h_H_eta.Fill(p.eta())
      if abs(p.pdgId()) == 5:
        h_b_pt.Fill(p.pt())
        h_b_eta.Fill(p.eta())
      if p.pdgId() == 22:
        h_y_pt.Fill(p.pt())
        h_y_eta.Fill(p.eta())

    h_photons = Handle("std::vector<flashgg::Photon>")
    event.getByLabel("flashggRandomizedPhotons", h_photons)
    print("N(photons) = {}".format(len(h_photons.product())))
    photons = h_photons.product()
    h_nphoton.Fill(len(photons))
    if(len(photons)>=2):
	    p4_g0 = ROOT.TLorentzVector()
            p4_g1 = ROOT.TLorentzVector()
            p4_g0.SetPtEtaPhiM(photons[0].pt(),photons[0].eta(),photons[0].phi(),photons[0].mass())
            p4_g1.SetPtEtaPhiM(photons[1].pt(),photons[1].eta(),photons[1].phi(),photons[1].mass())
            h_photon0_pt.Fill(photons[0].pt())
            h_photon1_pt.Fill(photons[1].pt())
            h_photon0_eta.Fill(photons[0].eta())
            h_photon1_eta.Fill(photons[1].eta())
	    if photons[0].pt() > 20 and photons[1].pt() > 20: 
		    h_reconsM_photon.Fill((p4_g0+p4_g1).M())
    for gam in photons:
      print('photon pt = {}'.format(gam.pt()))
      h_photon_eta.Fill(gam.eta())
      h_photon_pt.Fill(gam.pt())

    h_packedGenpar = Handle("std::vector<pat::PackedGenParticle>")
    event.getByLabel("flashggGenPhotons", h_packedGenpar)
    particles = h_packedGenpar.product()
    for p in particles:
      print('packedgenparticles P pid = {} status = {}'.format(p.pdgId(), p.status()))

    h_jets = Handle("std::vector<std::vector<flashgg::Jet> >")
    event.getByLabel("flashggFinalJets", h_jets)

    print( "N(jets) = %i" % len(h_jets.product()) )
    if len(h_jets.product()) <= 0: continue
    jets = h_jets.product()[0]
    h_njet.Fill(len(jets))
    if(len(jets)>=2):
	disc1=jets[0].bDiscriminator("pfDeepCSVJetTags:probb")+jets[0].bDiscriminator("pfDeepCSVJetTags:probbb")
	disc2=jets[1].bDiscriminator("pfDeepCSVJetTags:probb")+jets[1].bDiscriminator("pfDeepCSVJetTags:probbb")
	pt0=jets[0].pt()
        pt1=jets[1].pt()
        eta0=jets[0].eta()
        eta1=jets[1].eta()
	p4_b0 = ROOT.TLorentzVector()
        p4_b1 = ROOT.TLorentzVector()
        p4_b0.SetPtEtaPhiM(jets[0].pt(),jets[0].eta(),jets[0].phi(),jets[0].mass())
        p4_b1.SetPtEtaPhiM(jets[1].pt(),jets[1].eta(),jets[1].phi(),jets[1].mass())
	if disc1>0.4941 and disc2>0.4941:
             h_btaggedjet0_pt.Fill(pt0)
	     h_btaggedjet1_pt.Fill(pt1)
             h_btaggedjet0_eta.Fill(eta0)
             h_btaggedjet1_eta.Fill(eta1)
             if pt0 > 25 and pt1 > 25 and abs(eta0) <= 2.5 and abs(eta1) <= 2.5:
                 h_reconsM_bb.Fill((p4_b0+p4_b1).M())
    for jet in jets:
      print( "jet pt = %f DeepCSVBDisc = %f" % (jet.pt(), jet.bDiscriminator("pfDeepCSVJetTags:probb")+jet.bDiscriminator("pfDeepCSVJetTags:probbb")) )
      
      deepCSV = jet.bDiscriminator("pfDeepCSVJetTags:probb")+jet.bDiscriminator("pfDeepCSVJetTags:probbb")
      h_jet_pt.Fill(jet.pt())
      h_jet_eta.Fill(jet.eta())
      h_jet_deepcsv.Fill(deepCSV)
      
    h_dipho = Handle("std::vector<flashgg::DiPhotonCandidate>")
    event.getByLabel("flashggDiPhotons", h_dipho)

    diphos = h_dipho.product()
    print( "N(diphotons) = %i" %len(diphos) )
    for dipho in diphos:
      print( "Diphoton pt = %f" % dipho.pt() )
      h_dipho_mass.Fill(dipho.mass())

    nevents+=1
  fout.Write()
  fout.Close()

if __name__ == "__main__":
    main()
