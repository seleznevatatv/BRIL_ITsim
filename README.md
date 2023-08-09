This fork contains a specially modified code for the purpose of studying the linearity of the pixel cluster counting algorithm.
The Inner Tracker simulation data prepared by Cristina Oropeza in 2023 is used.

To obtain input files for rd53stream, in your working directory run the following commands:
```
cmsrel CMSSW_13_0_0_pre4
cd CMSSW_13_0_0_pre4/src
cmsenv
git clone git@github.com:seleznevatatv/BRIL_ITsim.git
```

Change input and output pathes in itdigiexporter_v2/python/ITdigiExporter.py
```
scram b -j8
cmsRun itdigiexporter_v2/python/ITdigiExporter.py
```
