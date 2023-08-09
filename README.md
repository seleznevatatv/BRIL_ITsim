This fork contains a specially modified code for the purpose of studying the linearity of the pixel cluster counting algorithm.
The Inner Tracker simulation data prepared by Cristina Oropeza in 2023 is used.

To obtain input files for rd53stream, in your working directory run the following commands:
```
cmsrel CMSSW_13_0_0_pre4
cd CMSSW_13_0_0_pre4/src
cmsenv
git clone
```



# The _not yet_ ultimate guide to BRIL Inner Tracker Simulations

### Preface
This guide is intended for people working on the BRIL Phase II upgrade, more specifically simulations of the Inner Tracker for Lumi measurements. It gives detailed instructions on how to set up and run simulations with custom IT geometries and varying levels of pileup. It is sectioned in the following parts:

1. Setting up an appropriate CMSSW environment
1. Using custom geometry files exported from [TkLayout](http://tklayout.web.cern.ch)
1. Generating Minimum Bias Events using the custom geometry to use as Pileup input for the rest of the simulation (**see Note below**)
1. Running the full simulation (**see Note below**)
    * locally using the runSim.sh script that wraps cmsDriver commands and a custom config file for more convenience
    * on the HT Condor batch system using the runSim.sh script and a submission file
1. Running the BRIL ITclusterAnalyzer on the output EDM file to populate BRIL relevant histograms for further analysis
    * locally
    * using the HT Condor batch system

**Note:** these sections describe how to produce private samples. However, there are centrally produced samples available that include not only the IT, but also information from all CMS sub-systems. If you want to use those files (**recommended**) skip through to the last section after setting up the environment and customizing the geometry files. 

### Setting up an Environment
It is preferrable to run these steps on a machine with access to the CERN cvmfs read-only filesystem as this allows to use different CMSSW releases -- this usually means lxplus. Before you get started, make sure that you have a few GB of available space in the directory you want to run in. This could either be /afs or /eos. The release this guide is based on is `CMSSW_10_6_0_patch2` that you set up like so:

```sh
source $VO_CMS_SW_DIR/cmsset_default.sh
mkdir mySimDir
cd mySimDir
cmsrel CMSSW_10_6_0_patch2
cd CMSSW_10_6_0_patch2/src/
```

This creates a CMSSW working area in the directory mySimDir/ where a directory tree is created in the `CMSSW_10_6_0_patch2/` folder. All your code and packages are always going to the `/src` subdirectory - it is the only place where you will be working.

```sh
cmsenv
```

This sets up the working environment, the paths and the right compiler. In the next step we need to check out some packages from the CMSSW github repo that we want to modify. Specifically, these will contain the modified geometry files. 

```sh
git cms-addpkg Configuration/PyReleaseValidation
git cms-addpkg Configuration/StandardSequences
git cms-addpkg Geometry/TrackerCommonData
git cms-addpkg Geometry/TrackerRecoData
git cms-addpkg Geometry/TrackerSimData
git cms-addpkg SLHCUpgradeSimulations/Geometry
```

### Using custom Geometries
Once done with the previous step, you have to get the custom geometry files and put them in the appropriate place. First, you have to check out this repo in your CMSSW work area `mySimDir/CMSSW_10_6_0_patch2/src` or, equivalently, `$CMSSW_BASE/src`.

```sh
cd $CMSSW_BASE/src
git clone https://github.com/gauzinge/BRIL_ITsim.git
cd BRIL_ITsim
```

or better yet fork the repo to your github account! Now you can call the `copyGeo.sh` with the following arguments: `-s=mysourceDir -d=mydestDir` to copy all relevant geometry files from the source (most likely `$CMSSW_BASE/src/BRIL_ITsim/ITGeometries/OT614_200_IT613` in this repo as source and `$CMSSW_BASE/src` as destination). Alternatively you can hardcode the paths in the `copyGeo.sh` script. If you use command line parsing you need to provide the absolute path. At the time of writing, IT geometry 6.1.3 is the latest and greatest. Full details [here](http://ghugo.web.cern.ch/ghugo/layouts/T15/OT616_200_IT613/layoutpixel.html).

```sh
source copyGeo.sh -s=$CMSSW_BASE/src/BRIL_ITsim/ITGeometries/OT614_200_IT613 -d=$CMSSW_BASE/src
```

The next step is crucial, so pay attention: since you modified the geometry files you have to re-build your added CMSSW sources. This is done using the scram command

```sh
scram b -j8
```

for 8 core compilation. Now everything you do will be based on the custom geometry.


### Using non-standard Pileup scenarios

The new `runSim.sh` script uses a cmsDriver config file that combines all the steps that `runTheMaxtrix.py` would call in a single python file - therefore it is not necessary to change any files for custom Pileup scenarios as the PU number is directly passed to the config from the `runSim.sh` script. I strongly advise against using the `runTheMatrix.py` command for anything at this stage.


### Generating Minimum Bias Events to use as Pileup Input in the actual Simulation

The next step is to generate a large enough sample of Minimum Bias events for your scenario that can later be used as Pileup input for the actual simulation. To do this, you need to run either the first step of the workflow generated by the CMS `runTheMatrix.py` command (**not recommended**) or the `runMinbias.sh` script (**recommended**). Just be sure to edit the script to have the appropriate output paths set correctly.

```sh
.runMinBias.sh NEVENTS JOBID
```

The first argument is the number of events to generate and the second is the jobid for batch processing (for now you can set it to 0).

Alternatively (**preferably**) you can use any of these files: `/afs/cern.ch/work/g/gauzinge/public/minBiasFiles/` which use geometry IT613. There are many such files and they can be provided to the pileup input as a list in the next stage.

#### Using runTheMatrix (not recommended)

It is important to choose the right scenario which in our case is the Phase II upgrade. The name is **2023D42**.

```sh
runTheMatrix.py --what upgrade -l 21240 -ne
```

will just dump the commands to your terminal but I suggest running the actual cmsDriver command instead:

```sh
cmsDriver.py MinBias_14TeV_pythia8_TuneCUETP8M1_cfi  --conditions auto:phase2_realistic -n 2000 --era Phase2 --eventcontent FEVTDEBUG --relval 90000,100 -s GEN,SIM --datatier GEN-SIM --beamspot HLLHC14TeV --geometry Extended2023D42 --nThreads 4
```

this generates 2000 Minimum Bias events in a file called ` MinBias_14TeV_pythia8_TuneCUETP8M1_cfi_GEN_SIM.root` -- put it somewhere where you won't accidentially delete it:

```sh
cd $CMSSW_BASE/src
mkdir myMinBiasSample
mv MinBias_14TeV_pythia8_TuneCUETP8M1_cfi_GEN_SIM.* myMinBiasSample/
```

Congratulations, you are almost done! In order to use these events you have to edit a file in `$CMSSW_BASE/src/Configuration/PyReleaseValidation/python/relval_steps.py`. The easiest is to search for the string `PUData` in your editor. Then edit these lines:

```python
PUDataSets={}
for ds in defaultDataSets:
    key='MinBias_14TeV_pythia8_TuneCUETP8M1'+'_'+ds
    name=baseDataSetReleaseBetter[key]
#    if '2017' in name or '2018' in name:
    if '2017' in name:
        PUDataSets[ds]={'-n':10,'--pileup':'AVE_35_BX_25ns','--pileup_input':'das:/RelValMinBias_13/%s/GEN-SIM'%(name,)}
    elif '2018' in name:
        PUDataSets[ds]={'-n':10,'--pileup':'AVE_50_BX_25ns','--pileup_input':'das:/RelValMinBias_13/%s/GEN-SIM'%(name,)}
    else:
        PUDataSets[ds]={'-n':10,'--pileup':'AVE_35_BX_25ns','--pileup_input':'das:/RelValMinBias_14TeV/%s/GEN-SIM'%(name,)}
```

change the last line to point to your Minimum Bias sample file (be sure to use an absolute path from your home direcotory - otherwise batched simulations won't work):

```python
PUDataSets[ds]={'-n':10,'--pileup':'AVE_35_BX_25ns','--pileup_input':'file:myMinBiasDataFile.root'}
```

this will tell CMSSW to use your Minimum Bias Data sample as Pileup input instead of some files with a standard geometry from the database. This is really important. Next, you have to once more compile your work environment to make it pick up the changes:

```sh
scram b -j8
```

Congratulations, you are all set for running the actual simulation!


### Running the full Simulation
There are 3 possible ways to run the full simulation. The normal way using the `runTheMatrix.py` command, a wrapper script called `runSim.sh` that is part of this repo or in batch mode on the CERN HT Condor batch service. The details will be described in the following sections!


#### runTheMatrix - just to validate your working environment
Now, to test your working environment, you can use the `runTheMatrix.py` script provided with CMSSW. Again, this will launch a 3 step process of generating events of a certain type, simulating the detector response and reconstruction of the events. More specifically, the pixel clustering happens in step3 (the RECO step). The scenario we want to use for our purposes is a single Neutrino (so an empty detector) overlaid with a variable number of pileup. The scenario has the label **21461**. So in your CMMSW `src` directory run:

```sh
runTheMatrix.py --what upgrade -l 21461 --command "-n 100"
```

The `--command "-n 100"` specifies 100 events. Change it if you want but be patient. Oh, and the default pileup number is 35. You can't easily change it in this workflow but it is essentially just a test to test the environment. Wait for the process to complete and then check your working directory: you should see files step1.root, step2.root and step3.root. The file you are interested in is **always** `step3.root`. Try verifying that it contains a collection of type `SiPixelCluster`:

```sh
edmDumpEventContent step3.root | grep SiPixelCluster
```

See something? Great, you are done for this part. 


#### Running simulations locally using the runSim.sh script (recommended)

Since the above workflow is tedious and does not really give you fine control over the options (plus it runs a bunch of processes that we don't need) it is usually better to use the `runSim.sh` script provided with this repo. It has some filepaths hardcoded (for example the path to put the output files and the pileup input files) so you want to open it in your editor and fix all the paths at the top of the script - you may also want to change the number of threads to something reasonable for your job. 

Before you continue however, you need to sandbox your CMSSW working environment (you could use the sandbox of this repo if you did only make the changes listed above). This is important for batch processing but we'll also use it for this case. To do that, create a new directory somewhere on your /afs or /eos and copy the `CMSSW_10_6_0_patch2` directory there. Alternatively you can just use the sandbox provided with the repo but it's the geometry mentioned above so be prudent. You'll need to repeat the steps below for each new geometry.

```sh
cd
mkdir test
cp -rf $CMSSW_BASE test
```

next, cd to that directory and remove all unnecessary files (don't be shy, you are working on a copy).

```sh
cd ~/test/CMSSW_10_6_0_patch2/src/
rm -rf BRIL_ITsim
rm -rf myMinBiasSample
```

Next, tar up your sandbox and copy it to your working area `$CMSSW_BASE/src/BRIL_ITsim/`:

```sh
cd test
tar -jcf sandbox.tar.bz2 CMSSW_10_6_0_patch2/
cp sandbox.tar.bz2 $CMSSW_BASE/src/BRIL_ITsim
```

Now you are ready to run the script. The `runSim.sh` script takes 3 command line arguments

```sh
./runSim.sh PU NEVENTS JOBID
```

The first one is the PU you want to run and any number can be specified, the second one is the number of events and the third one is a jobid (irrelevant for this case but important in batch mode) - please always provide all three. For now you can set the jobid to 0. This runs all three steps of the `runTheMatrix.py` command in a single process and in addition slims the output `step3.root` file down to only contain pixel relevant collections. The resulting output file is called `step3_pixel_PU_${PU}.${JOBID}.root` where ${PU} and ${JOBID} are replaced with your command line arguments. This is the file we will use as input to the ITclusterAnalyzer later.


#### Running simulations on CERN HT Condor batch service (recommended)

Before you do anything, first read the HT Condor [guide](http://batchdocs.web.cern.ch/batchdocs/tutorial/introduction.html) including chapter 3. 

Done? Good, now you know about submission files. Have a look at the `generatePU.sub` file provided with this repo. Change the PU parameter variable, the `NEvents` variable, the `request_cpu` number to match your number of threads in the `runSim.sh` script and the number of parallel jobs. Try with a small number of events and only one job first to verify that things work ok.

```sh
condor_submit generatePU.sub
```

and be patient and watch as your quota goes away... Happy simulating!

There is also a wrapper script for the generation of the MinBias sample called `runMinBias.sh` and a corresponding `generateMinBias.sub` submission file. Functionality is exactly the same except that `runMinBias.sh` takes only the number of events as command line argument. As always, check all variables and paths in these scripts before running.

Depending on the number of events you want to simulate, you should generate N minBias events where: *N = Nevents * PU * OOT_PU_range*

The PU is your max PU number so most likely 200, the *OOT_PU_range* is the number of BX that are simulated around the actual event to account for out-of-time Pileup. The default value is 15. Since N in the above quickly explodes, in practice it should be ok to generate a factor of 10 more minBias events than events you want to have as statistics. So for 10k events, use 100k minBias events.

**Note:** Private samples for various PU values are available in: `/eos/user/g/gauzinge/PUdata/`. These contain only pixel information and were simulated using geometry **IT613** and geometry scenario **2023D21**.


#### Centrally produced samples

They contain the full detector information and are available in the Grid. The full list can be obtained from the CMS DAS [here](https://cmsweb.cern.ch/das/request?view=list&limit=50&instance=prod%2Fglobal&input=dataset+status%3D*+dataset%3D%2FRelValNuGun%2FCMSSW_10_6_0_patch2*2023D42*%2F*).

**Important:** since these samples are stored in the Grid, you need to get a Grid certificate. See [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookStartingGrid) for instructions on how to do so.


###  Running the BRIL IT Cluster analyzer

The last step for now is running the ITclusterAnalyzer package written for CMSSW. It runs over the simulated files and creates some histograms relevant for BRIL studies. It is a standard EDAnalyzer class in CMSSW and the source code can be found in `$CMSSW_BASE/src/BRIL_ITsim/ITclusterAnalyzer/plugins`. Have a look. 

#### Running locally

The config file is in `$CMSSW_BASE/src/BRIL_ITsim/ITclusterAnalyzer/python/ITclusterAnalyzer_cfg.py`. You should edit this file to contain the right path for the file you want to analyze (to access the files from the centrally produced samples you need to use [xrootd](https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookXrootdService)) and you can modify the cuts for coincidence searches if you want - you can also disable those. The file should be pretty self explanatory.

For the _official_ samples, you need to create a proxy before attempting to run because otherwise you will not be able to access them and the job will fail.

```sh
voms-proxy-init -voms cms --valid 168:00
```

Then you can do:

```sh
cd BRIL_ITsim
cmsenv
scram b -j8
cmsRun ITclusterAnalyzer/python/ITclusterAnalyzer_cfg.py
```

This will produce a `summary.root` file that you can inspect like so:

```sh
cmsenv
root summary.root
new TBrowser
```

Have a look around, try to understand the generated histograms and do with them what you like.

In case you want to process multiple files you can also use the scripts provided in this repo. The script `runAnalysis.sh` works for our private samples, whereas the script `new_runAnalysis.sh` should be used for the centrally produced samples. As usual, make sure the paths are ok in the file (meaning you need to change them!). It will loop over all input files and generate the summary.root files for each of them and merge them in the end. The only command line argument is the Pileup step to process (pay attention to provide a double, so for PU 10 use CMD line argument 10.0):

```sh
./new_runAnalysis.sh 10.0
```

to process all files for PU 10 and generate a single summary file for this pileup value.

**Warning:** currently, for the centrally produced samples, one needs to manually enter the name of the corresponding input text file containing the names of all files for a given PU value. This will be fixed shortly!

#### Running the analysis on the HT Condor batch

**Under construction!**

If you have any bug reports or questions, contact me!
