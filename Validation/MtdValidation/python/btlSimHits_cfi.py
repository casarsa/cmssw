import FWCore.ParameterSet.Config as cms
from Validation.MtdValidation.btlSimHitsDefault_cfi import btlSimHitsDefault as _btlSimHitsDefault
btlSimHits = _btlSimHitsDefault.clone()

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(btlSimHits, useCrossingFrame = cms.bool(False), inputTag = cms.InputTag("g4SimHits", "FastTimerHitsBarrel"))
