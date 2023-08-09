# IT digi extraction for rate studies

## Setup

```bash
CMSSW_VERSION="CMS_13_0_10"
cmsrel "${CMSSW_VERSION}"
cd "${CMSSW_VERSION}/src"
cmsenv
mkdir BRIL_IT_Analysis
cd BRIL_IT_Analysis
git clone https://gitlab.cern.ch/mbensink/itdigiexporter.git
scram b -j8
```

## Running

In short:

```bash
cmsRun itdigiexporter/python/ITdigiExporter.py \
       dataset=PU_75_cbarrera.txt \
       output=/path/to/output_file.root
```

Some notes:

- The dataset name must match to a file name in `./data/datasets`. The file should contain file paths to the input files for a given data set.
- When defining a new data set name, note that the format should always be `<PhysicsProcess>_<CMSSWVersion>_<DetectorGeometryTag>_<PileUpTag>`
- The choice of cabling map is made automatically based on the `DetectorGeometryTag` (D41 in the example above). Look at the `GEOMETRY_TO_CABLING` dictionary in the python config to see / change.
- Note that the cabling maps are automatically loaded from `./data/cabling/`. Make sure that your desired map is available there. If you need to create a new cabling sqlite file, go [here](https://github.com/AndreasAlbert/cmssw/tree/2021-09-30_itrate_cabling_maps/CondTools/SiPhase2Tracker/test).
- Make sure to use a CMSSW version that contains the geometry used in the MC sample being analyzed. The geometry files are periodically added / removed from CMSSW releases, so any releases that is either much newer or much older than the sample of interest will lead to problems. Additionally, CMSSW EDM files are not guaranteed to be readable by versions older than the one used to write the file.
- Form large datasets, use crab to process. See examples in the crab/ folder.
