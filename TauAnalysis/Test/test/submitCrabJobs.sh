voms-proxy-init -voms cms -valid 100:0

rm -r /afs/cern.ch/work/v/veelken/CMSSW_7_0_x/crab/TauAnalysis_Test/crabdir_producePFTauIdEffNtuple2_simZplusJets20PUat25ns

crab -create -cfg crabProducePFTauIdEffNtuple_ZplusJets20PUat25ns.cfg
crab -submit 1-500 -c /afs/cern.ch/work/v/veelken/CMSSW_7_0_x/crab/TauAnalysis_Test/crabdir_producePFTauIdEffNtuple2_simZplusJets20PUat25ns
crab -submit 501-1000 -c /afs/cern.ch/work/v/veelken/CMSSW_7_0_x/crab/TauAnalysis_Test/crabdir_producePFTauIdEffNtuple2_simZplusJets20PUat25ns
