# imports
import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing

# create a new CMS process
process = cms.Process("binsize")

# set up the options
options = VarParsing.VarParsing('analysis')
#set up the defaults
options.inputFiles = 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.0TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.2TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.3TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.4TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.5TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.6TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.7TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.8TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.12TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.13TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.14TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.15TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.16TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.17TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.18TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.19TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.20TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_200.0.21TkOnly.root'
# options.inputFiles = 'file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.0TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.4TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.5TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.6TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.7TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.8TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.9TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.10TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.11TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.12TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.13TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.14TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.15TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.16TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.17TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.18TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.20TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.21TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.22TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.23TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.24TkOnly.root,file:/eos/cms/store/group/dpg_bril/comm_bril/phase2-sim/tkOnly/step3_pixel_PU_2.0.25TkOnly.root'
options.outputFile='binsize200.root'
options.maxEvents = -1 #all events

#get and parse command line arguments
options.parseArguments()

# load the geomtry that i modified
# process.load('Configuration.Geometry.GeometryExtended2023D21Reco_cff')
# process.load('Configuration.Geometry.GeometryExtended2023D42Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2026D63Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')

# initialize MessageLogger and output report
process.load("FWCore.MessageLogger.MessageLogger_cfi")
# process.MessageLogger.cerr.threshold = 'INFO'
# process.MessageLogger.categories.append('binsize')
# process.MessageLogger.cerr.INFO = cms.untracked.PSet(
    # limit=cms.untracked.int32(-1)
# )
process.MessageLogger.cerr.FwkReport.reportEvery = cms.untracked.int32(200)

process.options = cms.untracked.PSet(wantSummary=cms.untracked.bool(False)
                                    # ,SkipEvent = cms.untracked.vstring('ProductNotFound')
                                    )

process.maxEvents = cms.untracked.PSet(input=cms.untracked.int32(options.maxEvents))

# the input file
process.source = cms.Source("PoolSource",
                            fileNames=cms.untracked.vstring(options.inputFiles),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck')
                            )

process.content = cms.EDAnalyzer("EventContentAnalyzer")
# the config of my analyzer
process.BRIL_IT_Analysis = cms.EDAnalyzer('binsize',
                                         clusters=cms.InputTag("siPixelClustersPreSplitting"),
                                         # clusters=cms.InputTag("siPixelClusters"),
                                         # digis=cms.InputTag("simSiPixelDigis", "Pixel", "FULLSIM"),
                                         # digis=cms.InputTag("simSiPixelDigis"),
                                         maxBinTEPX=cms.untracked.uint32(50000),
                                         maxBinD4R1=cms.untracked.uint32(250000),
                                         # maxBinD4R1=cms.untracked.uint32(50000),
                                         bxperorbit=cms.untracked.uint32(2500),
                                         trgfrqTEPX=cms.untracked.uint32(75),
                                         trgfrqD4R1=cms.untracked.uint32(825),
                                         )

# the TFIleService that produces the output root files
process.TFileService = cms.Service("TFileService",
                                   fileName=cms.string(options.outputFile)
                                   )


# process.p = cms.Path( ... process.content * ...  )
process.p = cms.Path(process.BRIL_IT_Analysis)
