# imports
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import numpy as np

# create a new CMS process
process = cms.Process("ITdigiExporter")

# set up the options
options = VarParsing.VarParsing('analysis')
#set up the defaults

# Read file list
# files = [ "root://cms-xrd-global.cern.ch///" + x for x in np.loadtxt("NuGun_D41PU200.txt",dtype=str) ]
# options.inputFiles = files
options.inputFiles = 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.2TkOnly.root'


options.maxEvents = 1000 #all events

#get and parse command line arguments
options.parseArguments()

# fileIndex = int(str(options.inputFiles)[-8])
# offset = fileIndex*1000;
# print("This should be file",str(fileIndex),"for the given PU number and yield an offset of",str(offset))
offset=0

# load the geomtry that i modified
# process.load('Configuration.Geometry.GeometryExtended2022D21Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D63Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# process.MessageLogger.cerr.threshold = 'INFO'
# process.MessageLogger.categories.append('ITdigiExporter')
# process.MessageLogger.cerr.INFO = cms.untracked.PSet(
#     limit=cms.untracked.int32(-1)
# )

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False))

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

# the input file
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(options.inputFiles)
                            )

# the config of my analyzer
process.BRIL_IT_Analysis = cms.EDAnalyzer('ITdigiExporter',
                                         eventNoOffset=cms.int32(offset),
                                         clusters=cms.InputTag("siPixelClustersPreSplitting"),
                                         digis=cms.InputTag("simSiPixelDigis", "Pixel", "FULLSIM"),
                                         simlinks=cms.InputTag("simSiPixelDigis", "Pixel", "FULLSIM")
                                         )


process.TFileService = cms.Service("TFileService", fileName = cms.string("file:/eos/user/g/gauzinge/PUdata/itdata_PU200_1.root") )
process.p = cms.Path(process.BRIL_IT_Analysis)
