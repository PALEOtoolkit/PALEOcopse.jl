# PALEOcopse Reactions


```@meta
CurrentModule = PALEOcopse.Forcings
```
## Forcings
```@docs
COPSEForcings.ReactionForce_CK_Solar
COPSEForcings.ReactionForce_UDWEbergman2004
COPSEForcings.COPSEForcings.ReactionForce_Bbergman2004
COPSEForcings.ReactionForce_CPlandrelbergman2004
COPSEForcings.ReactionForce_spreadsheet
LipForcing.ReactionForce_LIPs
```

```@meta
CurrentModule = PALEOcopse.COPSE
```
## COPSE
```@docs
ModelBergman2004.ReactionModelBergman2004
CIsotopes.ReactionCIsotopes
MapAtmOceanReservoirs.ReactionAtmOcean_O
MapAtmOceanReservoirs.ReactionAtmOcean_A
```

```@meta
CurrentModule = PALEOcopse.Global
```
## Global
```@docs
Temperature.ReactionGlobalTemperatureCK1992
Temperature.ReactionGlobalTemperatureBerner
```

## BioGeoChem
```@meta
CurrentModule = PALEOcopse.BioGeoChem
```

### Isotope systems

#### Strontium
```@docs
Strontium
Strontium.ReactionSrMantleCrust
Strontium.ReactionSrSed
Strontium.ReactionSrLand
Strontium.ReactionSrOceanfloor
```

## Land
```@meta
CurrentModule = PALEOcopse.Land
```
```@docs
LandBergman2004.ReactionLandBergman2004
LandArea.ReactionLandArea
LandCOPSEReloaded.ReactionLandBiota
LandCOPSEReloaded.ReactionLandWeatheringRates
LandCOPSEReloaded.ReactionLandWeatheringFluxes
```

## Ocean
```@meta
CurrentModule = PALEOcopse.Ocean
```

```@docs
OceanCOPSE.ReactionMarineBiotaCOPSE
OceanCOPSE.ReactionOceanBurialCOPSE
```

## Ocean floor
```@meta
CurrentModule = PALEOcopse.Oceanfloor
```

```@docs
CarbBurial.ReactionCarbBurialAlk
Weathering.ReactionSeafloorWeathering
```

## Sedcrust
```@meta
CurrentModule = PALEOcopse.SedCrust
```

```@docs
SedCrustCOPSE.ReactionSedCrustCOPSE
```
