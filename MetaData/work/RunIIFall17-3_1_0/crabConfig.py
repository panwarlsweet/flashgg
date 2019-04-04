from WMCore.Configuration import Configuration
config = Configuration()
import os
config.section_("General")
config.General.requestName = "DiHiggs_SM"
config.General.transferLogs = False
config.section_("JobType")
config.JobType.pluginName = "Analysis"
config.JobType.psetName = "microAODstd.py"
## to include local file in the sendbox, this will put the file in the directory where cmsRun runs
config.JobType.inputFiles = ['Fall17_17Nov2017BCDEF_V6_DATA.db','Fall17_17Nov2017_V6_MC.db','QGL_cmssw8020_v2.db']
## incrase jobs time wall, maximum 2750 minutes (~46 hours)
config.JobType.maxJobRuntimeMin = 2750
config.JobType.maxMemoryMB = 2500 # For memory leaks. NB. will block jobs on many sites
## config.JobType.scriptExe = "cmsWrapper.sh"
config.JobType.pyCfgParams = ['datasetName=/GluGluToHHTo2B2G_node_SM_13TeV-madgraph_bReg/RunIIFall17MiniAODv2-PU202017_12Apr2018_94X_mc2017_realistic_v14-v2', 'processType=bkg']
config.JobType.sendPythonFolder = True
config.JobType.sendExternalFolder = True
config.section_("Data")
config.Data.inputDataset = "/GluGluToHHTo2B2G_node_SM_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM"
config.Data.inputDBS = 'global'
config.Data.splitting = 'EventAwareLumiBased'
config.Data.unitsPerJob = 20000
config.Data.totalUnits   =  -1
config.Data.publication = True
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'RunIIFall17-3_1_0-v0-RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2'
config.Data.outLFNDirBase = "/store/group/phys_higgs/HiggsExo/HH_bbgg/RunIIFall17-3_1_0/3_1_0"
config.section_("Site")
config.Site.storageSite = "T2_CH_CERN"

