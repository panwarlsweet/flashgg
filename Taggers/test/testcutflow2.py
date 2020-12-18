import ROOT
import os,sys
sum_=0
list_=[]
#!/usr/bin/env python                                                                          
import sys, os, shutil, re, subprocess
import ROOT
from DataFormats.FWLite import Events, Handle
import numpy as np
import json

mass=sys.argv[1]
json_no=sys.argv[2]
year=sys.argv[3]
year=int(year)
if year == 2016:
  campaign="Era2016_RR-17Jul2018"
  dasfile="/NMSSM_XToYHTo2b2g_MX-"+mass+"_TuneCUETP8M1_13TeV-madgraph-pythia8/lata-Era2016_RR-17Jul2018_v2_p2-v2_p2-v0-RunIISummer16MiniAODv3-PUMoriond17_RP_94X_mcRun2_asymptotic_v3-v1-58814a4778ef9be1eca50f9befd9e8ad/USER" 
elif year == 2017:
  campaign="Era2017_RR-31Mar2018"
  dasfile="/NMSSM_XToYHTo2b2g_MX-"+mass+"_TuneCP5_13TeV-madgraph-pythia8/lata-Era2017_RR-31Mar2018-v2_p11-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_RP_94X_mc2017_realistic_v14-v1-7ffdea9265ec39789c1251c92052ae33/USER"
elif year == 2018:
  campaign="Era2018_RR-17Sep2018"
  dasfile="/NMSSM_XToYHTo2b2g_MX-"+mass+"_TuneCP5_13TeV-madgraph-pythia8/lata-Era2018_RR-17Sep2018-v2_p12-v0-RunIIAutumn18MiniAOD-RP_102X_upgrade2018_realistic_v15-v1-7b6feb31f1c715e59da0e42d28340658/USER"

f = open('flashgg/MetaData/data/'+str(campaign)+'_v2_p2/datasets_NMSSM_'+json_no+'.json')
data = json.load(f)

from FWCore.ParameterSet.VarParsing import VarParsing

for i,fi in enumerate(data[dasfile]["files"]):
  neve_f = fi["events"]
  f_cutflow=ROOT.TFile.Open("cutflow_"+str(year)+"_NMSSM_XToYHTo2b2g_MX-"+mass+"_%i.root"%i)
  h_cutflow = f_cutflow.Get("h_cutflow")
  if h_cutflow.GetBinContent(1) != neve_f:
    print(i, h_cutflow.GetBinContent(1), neve_f)
    list_.append(i)
    #print(i, h_cutflow.GetBinContent(1), genmbb.Integral())
    #print(i, end=" ")
    #sum_+=genmbb.Integral() - h_cutflow.GetBinContent(1)
  else: continue
print(sum_)
print(list_, len(list_))
  
