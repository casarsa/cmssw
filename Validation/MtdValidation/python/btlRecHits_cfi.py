import FWCore.ParameterSet.Config as cms
from Validation.MtdValidation.btlRecHitsDefault_cfi import btlRecHitsDefault as _btlRecHitsDefault
btlRecHits = _btlRecHitsDefault.clone()

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(btlRecHits, useCrossingFrame = cms.bool(False), simHitsTag = cms.InputTag("g4SimHits", "FastTimerHitsBarrel"))
