import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import re
import os


def get_files(dataset, redirector):
    with open("itdigiexporter/data/datasets/" + dataset + ".txt", "r") as infile:
        files = [redirector + x.strip() for x in infile.readlines()]
    print("Processing the following input files: ", files)
    return files


## Translation of overall geometry tags (D*)
## to tracker geometry tags identifying the CABLING
## NB: The tracker geometry tag for CABLING is not
## necessarily the same as the tracker geometry tag
## included in the definition of the overall geometry here:
## See https://github.com/cms-sw/cmssw/blob/master/Configuration/Geometry/README.md
## This happens in cases where two tracker geometry tags have the
## same cabling map, but that map is officially associated
## to only one of them.
GEOMETRY_TO_CABLING = {
    "D21": "OT616_200_IT613",  # PLACEHOLDER, I DO NOT KNOW THE REAL ONE YET
    "D41": "OT616_200_IT613",  # T14
    "D76": "OT800_IT701",  # T21, actual geometry name is 'OT800_IT615', but cabling is equivalent
    "D79": "OT800_IT701",  # T23, actual geometry name is 'OT800_IT700', but cabling is equivalent
    "D80": "OT800_IT701",  # T25, actual geometry name is 'OT800_IT702', but cabling is equivalent
    "D81": "OT800_IT701",  # T25, actual geometry name is 'OT800_IT702', but cabling is equivalent
    "D88": "OT800_IT701",  # T25, actual geometry name is 'OT800_IT702', but cabling is equivalent
    "D91": "OT801_IT640",  # T25, actual geometry name is 'OT800_IT702', but cabling is equivalent
}


def get_global_geometry_tag(dataset):
    """
    Extracts global geometry tag (e.g. D41) from data set name.
    """
    m = re.match(".*_(D\d+).*", dataset)
    assert m, "Could not find geometry tag in data set name: " + dataset
    assert len(m.groups()) == 1, (
        "Found multiple matches for geometry tag in data set name: " + dataset
    )
    return m.groups()[0]


def get_cabling_map(dataset):
    """
    Returns the path to the cabling SQLite file for a given data set name.

    The file is found based on the geometry tag in the dataset name.
    Refer to the GEOMETRY_TO_CABLING dictionary to see the mapping.
    """
    fname = None

    global_geometry_tag = get_global_geometry_tag(dataset)
    assert global_geometry_tag in GEOMETRY_TO_CABLING, (
        "Cabling geometry unknown for global geometry " + global_geometry_tag
    )

    tracker_geo_for_cabling = GEOMETRY_TO_CABLING[global_geometry_tag]

    fname = "itdigiexporter/data/cabling/ITDTCCablingMap_{TKGEO}.db".format(
        TKGEO=tracker_geo_for_cabling
    )
    if fname is None:
        raise RuntimeError("Failed to find cabling map for data set name: " + dataset)

    assert os.path.exists(fname), "Cabling map sqlite file does not exist: " + fname

    return "sqlite_file:" + fname


def get_geometry_import(dataset):
    global_geometry_tag = get_global_geometry_tag(dataset)
    return "Configuration.Geometry.GeometryExtended2026{TAG}Reco_cff".format(
        TAG=global_geometry_tag
    )


# create a new CMS process
process = cms.Process("ITdigiExporter")

# Command line arguments
options = VarParsing.VarParsing("analysis")

options.register(
    "dataset",
    "TTBar14TeV_10_6_0_patch2_D41_PU200",
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.string,  # string, int, or float
    "Dataset to process. Must correspond to file list name in ./itdigiexporter/data/datasets/",
)

options.register(
    "xrdredirector",
    "root://cmsxrootd.fnal.gov///",
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.string,  # string, int, or float
    "XRootD redirector string for file access.",
)
options.register(
    "output",
    "output.root",
    VarParsing.VarParsing.multiplicity.singleton,  # singleton or list
    VarParsing.VarParsing.varType.string,  # string, int, or float
    "Path where output file will be created. Directory will be created if it does not exist. {DATASET} tag will be replaced by dataset name.",
)

options.parseArguments()
options.inputFiles = get_files(options.dataset, options.xrdredirector)

# Load CondDB service
process.load("CondCore.CondDB.CondDB_cfi")

cabling_map_file = get_cabling_map(options.dataset)
print("Picked cabling file: " + cabling_map_file)
process.CondDB.connect = cabling_map_file

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))
process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

# the input file
process.source = cms.Source(
    "PoolSource",
    fileNames=cms.untracked.vstring(options.inputFiles),
    duplicateCheckMode=cms.untracked.string("checkEachFile"),
)

# number of events offset to process multiple input files and avoid duplicates
offset = 100

# Pixel Geometry
process.load(get_geometry_import(options.dataset))
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("FWCore.MessageLogger.MessageLogger_cfi")


process.PoolDBESSource = cms.ESSource(
    "PoolDBESSource",
    process.CondDB,
    DumpStat=cms.untracked.bool(True),
    toGet=cms.VPSet(
        cms.PSet(
            record=cms.string("TrackerDetToDTCELinkCablingMapRcd"),
            tag=cms.string("DTCCablingMapProducerUserRun"),
        )
    ),
)

# the config of my analyzer
process.BRIL_IT_Analysis = cms.EDAnalyzer(
    "ITdigiExporter",
    fileSize=cms.int32(offset),
    digis=cms.InputTag("simSiPixelDigis", "Pixel", "RECO"),
    clusters=cms.InputTag("siPixelClusters"),
)

process.TFileService = cms.Service("TFileService", fileName=cms.string(options.output))

process.p = cms.Path(process.BRIL_IT_Analysis)
