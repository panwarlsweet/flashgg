#!/usr/bin/env python

from __future__ import print_function
import sys, os, shutil, re, subprocess
import ROOT
from DataFormats.FWLite import Events, Handle
import numpy as np
import json

#jsonfile=sys.argv[1]                                                                         
f = open('../../MetaData/data/Era2016_RR-17Jul2018_v2_p2/datasets_NMSSM_0.json')             
data = json.load(f)                                                                          
dasfile="/NMSSM_XToYHTo2b2g_MX-400_TuneCUETP8M1_13TeV-madgraph-pythia8/lata-Era2016_RR-17Jul2018_v2_p2-v2_p2-v0-RunIISummer16MiniAODv3-PUMoriond17_RP_94X_mcRun2_asymptotic_v3-v1-58814a4778ef9be1eca50f9befd9e8ad/USER"
mass="400"

def main():

  from FWCore.ParameterSet.VarParsing import VarParsing
  
  options = VarParsing ('analysis')
  #options.inputFiles  = ["root://cms-xrd-global.cern.ch//store/user/lata/t3store2/NMSSM_microAODs/2016/NMSSM_XToYHTo2b2g_MX-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/Era2016_RR-17Jul2018_v2_p2-v2_p2-v0-RunIISummer16MiniAODv3-PUMoriond17_RP_94X_mcRun2_asymptotic_v3-v1/200416_111337/0000/myMicroAODOutputFile_121.root"]
  #f="/store/user/lata/t3store2/NMSSM_microAODs/2016/NMSSM_XToYHTo2b2g_MX-1000_TuneCUETP8M1_13TeV-madgraph-pythia8/Era2016_RR-17Jul2018_v2_p2-v2_p2-v0-RunIISummer16MiniAODv3-PUMoriond17_RP_94X_mcRun2_asymptotic_v3-v1/200416_111337/0000/myMicroAODOutputFile_121.root"
  options.inputFiles  = ["root://cms-xrd-global.cern.ch/"]
  options.maxEvents = -1
  options.outputFile = "fout.root"
  options.parseArguments()

  ###### save hisots in root file ######
  for i,fi in enumerate(data[dasfile]["files"]):
    fout = ROOT.TFile("output_M%s_%i.root"%(mass,i), 'RECREATE')
    fout.cd()
    ##### define Histograms ######
    h_gen_Mdiphoton = ROOT.TH1F("genmgg",";gen m_{H} [GeV];Events;;",1000, 124.,126.)
    h_gen_Mbb = ROOT.TH1F("genmbb",";gen m_{Y} [GeV];Events;;", 1500, 50.,1550.)
    
    events = Events("root://cms-xrd-global.cern.ch/"+str(fi["name"]))
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
          h_gen_Mdiphoton.Fill(p.mass())
        elif p.pdgId() == 35:
          h_gen_Mbb.Fill(p.mass())
        else: continue

      nevents+=1
    fout.Write()
    fout.Close()

if __name__ == "__main__":
    main()
