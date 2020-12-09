```
export DISPLAY=localhost:0.0
ssh -XY lxplus6
cd ~/public/mySimDir/CMSSW_10_4_0_pre2/src/BRIL_ITsim
cmsenv
scram b -j8
cmsRun ITdigiExporter/python/ITdigiExporter.py
```