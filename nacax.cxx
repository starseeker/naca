/*
MODULE NacaAuxilary
  ------------------------------------------------------------------------------
  PURPOSE -

  AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

  REVISION HISTORY
    DATE  VERS RERSON  STATEMENT OF CHANGES
  17Jun99  0.5   RLC   Original coding
   6Sep99  0.6   RLC   Added LeadingEdgeRadius function
   9Dec99  0.7   RLC   Moved all epsilon, psi tables to new module
  13Feb01  0.8   RLC   Added CombineThicknessAndCamber
  16Feb01  0.81  RLC   Added ParametrizeAirfoil, Combo6seriesMeanLine
  14Mar01  0.85  RLC   Added Polynomial,ScaleFactor;revised SetSixDigitPoints
  15Mar01  0.86  RLC   Added spline interpolation to 6-series thickness
  19Sep01  0.9   RLC   Revised LoadX - coarse is same as Abbot & vonDoenhoeff
  21Sep01  0.91  RLC   Added MEDIUM and FINE spacing to LoadX
  25Sep01  0.92  RLC   Added GetRk1 and GetRk1k2
  03Oct01  0.93  RLC   Corrected two errors in 3-digit reflex mean line
  13Oct01  0.94  RLC   Revised CalA1,CalA2,CalA3
  16Nov01  0.95  RLC   Getting the interpolation right
  24Nov01  0.96  RLC   Optional argument coding
  30Nov01  0.97  RLC   Added LeRadius4 and LeRadius4M
  06Dec01  0.98  RLC   Added LeRadius6
  26Dec01  0.99  RLC   Made sure we could do points dense at l.e.
   4Jan02  1.00  RLC   Final cleanup for PDAS 7
  16Jan09  1.05  RLC   Small cosmetic improvements
  07Nov10  1.10  RLC   Fixed bug in GetRk1 (thanks to Robert Stone)
*/
