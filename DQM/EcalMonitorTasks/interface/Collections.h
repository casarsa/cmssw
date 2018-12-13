#ifndef EcalDQMonitorTaskCollections_H
#define EcalDQMonitorTaskCollections_H

#include <string>

namespace ecaldqm {

  enum Collections {
    kSource,
    kEcalRawData,
    kEBGainErrors,
    kEEGainErrors,
    kEBChIdErrors,
    kEEChIdErrors,
    kEBGainSwitchErrors,
    kEEGainSwitchErrors,
    kTowerIdErrors,
    kBlockSizeErrors,
    kMEMTowerIdErrors,
    kMEMBlockSizeErrors,
    kMEMChIdErrors,
    kMEMGainErrors,
    kEBSrFlag,
    kEESrFlag,
    kEBDigi,
    kEEDigi,
    kPnDiodeDigi,
    kTrigPrimDigi,
    kTrigPrimEmulDigi,
    kEBUncalibRecHit,
    kEEUncalibRecHit,
    kEBLaserLedUncalibRecHit,
    kEELaserLedUncalibRecHit,
    kEBTestPulseUncalibRecHit,
    kEETestPulseUncalibRecHit,
    kEBRecHit,
    kEERecHit,
    kEBReducedRecHit,
    kEEReducedRecHit,
    kEBBasicCluster,
    kEEBasicCluster,
    kEBSuperCluster,
    kEESuperCluster,
    nCollections
  };

  std::string const collectionName[nCollections] = {
    "Source",
    "EcalRawData",
    "EBGainErrors",
    "EEGainErrors",
    "EBChIdErrors",
    "EEChIdErrors",
    "EBGainSwitchErrors",
    "EEGainSwitchErrors",
    "TowerIdErrors",
    "BlockSizeErrors",
    "MEMTowerIdErrors",
    "MEMBlockSizeErrors",
    "MEMChIdErrors",
    "MEMGainErrors",
    "EBSrFlag",
    "EESrFlag",
    "EBDigi",
    "EEDigi",
    "PnDiodeDigi",
    "TrigPrimDigi",
    "TrigPrimEmulDigi",
    "EBUncalibRecHit",
    "EEUncalibRecHit",
    "EBLaserLedUncalibRecHit",
    "EELaserLedUncalibRecHit",
    "EBTestPulseUncalibRecHit",
    "EETestPulseUncalibRecHit",
    "EBRecHit",
    "EERecHit",
    "EBReducedRecHit",
    "EEReducedRecHit",
    "EBBasicCluster",
    "EEBasicCluster",
    "EBSuperCluster",
    "EESuperCluster"
  };

}

#endif
