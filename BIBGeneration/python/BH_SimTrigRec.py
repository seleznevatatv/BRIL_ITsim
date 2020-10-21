# Auto generated configuration file
# using: 
# Revision: 1.163 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: Configuration/Generator/python/TTbar_cfi.py -s GEN,SIM,DIGI,L1,DIGI2RAW,HLT --eventcontent RAWSIM --datatier GEN-SIM-RAW --conditions FrontierConditions_GlobalTag,IDEAL_V11::All -n 10 --no_exec
import FWCore.ParameterSet.Config as cms

process = cms.Process('RECO')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
#process.load('Configuration.StandardSequences.GeometryExtended_cff')
process.load('Configuration.StandardSequences.GeometryExtendedNewPipe_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Sim_cff')
process.load('Configuration.StandardSequences.Digi_cff')
process.load('Configuration.StandardSequences.DigiToRaw_cff')
process.load('Configuration.StandardSequences.SimL1Emulator_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.RawToDigi_cff')
process.load('Configuration.StandardSequences.ReconstructionCosmics_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')
process.load("Configuration.StandardSequences.SimulationRandomNumberGeneratorSeeds_cff")

process.load('L1TriggerConfig.L1GtConfigProducers.L1GtTriggerMaskTechTrigConfig_cff')
process.load('HLTrigger/HLTfilters/hltLevel1GTSeed_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(100)
)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

# Input source
process.source = cms.Source("PoolSource",
                            noEventSort = cms.untracked.bool(True),
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            fileNames = cms.untracked.vstring('file:BeamHalo_GEN.root'),
                            skipEvents = cms.untracked.uint32(0)
)

# Output definition
process.output = cms.OutputModule("PoolOutputModule",
                                  splitLevel = cms.untracked.int32(0),
                                  outputCommands = cms.untracked.vstring( "keep *_*_*_*"),
                                  fileName = cms.untracked.string('file:BeamHalo_RECO_al.root'),
                                  dataset = cms.untracked.PSet(
    dataTier = cms.untracked.string('GEN-SIM-RAW'),
    filterName = cms.untracked.string('')
    )
)


process.RECOSIMEventContent.outputCommands.extend(
     ['keep *_*_*_GEN'])
process.RECOSIMEventContent.outputCommands.extend(
     ['keep *_trackingtruthprod_*_*'])
process.RECOSIMEventContent.outputCommands.extend(
     ['keep *_mergedtruth_MergedTrackTruth_*'])
process.RECOSIMEventContent.outputCommands.extend(
     ['keep *_g4SimHits_*_*'])
process.RECOSIMEventContent.outputCommands.extend(
     ['keep *_TriggerResults_*_*'])

# Additional output definition
import SimG4Core.Application.g4SimHits_cfi
process.g4SimHits.Generator.ApplyEtaCuts  = cms.bool(False)
process.g4SimHits.Generator.MinPCut  = cms.double(0.0001) #100keV
process.g4SimHits.Generator.BeamBkgdEvent = cms.untracked.bool(True)
process.g4SimHits.StackingAction.SaveFirstLevelSecondary = cms.untracked.bool(True)
process.g4SimHits.StackingAction.SavePrimaryDecayProductsAndConversionsInCalo = cms.untracked.bool(True)
process.g4SimHits.StackingAction.SavePrimaryDecayProductsAndConversionsInMuon = cms.untracked.bool(True)

#not present in CMSSW_11_2_X
#process.g4SimHits.SteppingAction.KillBeamPipe = cms.bool(False)

# Other statements
process.GlobalTag.globaltag = 'START42_V11::All'

# Path and EndPath definitions
process.simulation_step   = cms.Path(process.psim*process.mix)
process.digitisation_step = cms.Path(process.pdigi)
process.L1simulation_step = cms.Path(process.SimL1Emulator)
process.digi2raw_step     = cms.Path(process.DigiToRaw)

# Path and EndPath definitions
process.raw2digi_step       = cms.Path(process.RawToDigi)
process.localreco_step      = cms.Path(process.localReconstructionCosmics)
process.bhalo_step          = cms.Path(process.beamhaloTracksSeq)
process.reco_track_step     = cms.Path(process.trackerlocalreco)
process.reco_calo_step      = cms.Path(process.calolocalreco)
process.reco_muon_step      = cms.Path(process.muonsLocalRecoCosmics)
process.endjob_step         = cms.Path(process.endOfProcess)
process.out_step            = cms.EndPath(process.output)
    


# Schedule definition
#process.schedule = cms.Schedule(process.raw2digi_step,process.reco_track_step,process.reco_calo_step,process.reco_muon_step,process.endjob_step,process.out_step)




process.endjob_step       = cms.Path(process.endOfProcess)
process.out_step = cms.EndPath(process.output)

# Schedule definition

# Simulation and L1

process.schedule = cms.Schedule(process.simulation_step,
                                process.digitisation_step,
                                process.L1simulation_step,
                                process.digi2raw_step)

# High level trigger

#process.schedule.extend(process.HLTSchedule)


# Reconstruction

process.schedule.extend([process.raw2digi_step,process.localreco_step,process.bhalo_step])


process.schedule.extend([process.endjob_step,process.out_step])



