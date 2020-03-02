import FWCore.ParameterSet.Config as cms
from Validation.MtdValidation.etlSimHitsDefault_cfi import etlSimHitsDefault as _etlSimHitsDefault
etlSimHits = _etlSimHitsDefault.clone()

from Configuration.ProcessModifiers.premix_stage2_cff import premix_stage2
premix_stage2.toModify(etlSimHits, useCrossingFrame = cms.bool(False), inputTag = cms.InputTag("g4SimHits", "FastTimerHitsEndcap"))
