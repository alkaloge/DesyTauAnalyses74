voms-proxy-init -voms cms -valid 100:0

rm -r /afs/cern.ch/work/v/veelken/CMSSW_7_0_x/crab/TauAnalysis_Skimming/crabdir_copyToCastor_simZplusJets20PUat25ns

crab -create -cfg crabCopyToCastor.cfg
crab -submit -c /afs/cern.ch/work/v/veelken/CMSSW_7_0_x/crab/TauAnalysis_Skimming/crabdir_copyToCastor_simZplusJets20PUat25ns

