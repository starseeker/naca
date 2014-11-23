/*
  NACA++ - C++ conversion of NACA456
  ------------------------------------------------------------------------------
  AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software

  REVISION HISTORY
    DATE  VERS PERSON  STATEMENT OF CHANGES
   2May01  0.1   RLC   Built first combined program from existing components
   3May01  0.11  RLC   xxx
   6Sep01  0.2   RLC   Fixed MakeFinalAirfoil; started ProcessGraphics
   9Sep01  0.3   RLC   Added FitSymmetric and FitAsymmetric
  19Sep01  0.4   RLC   Revised spline logic
  25Sep01  0.5   RLC   Removed the combo CL mean line stuff
  10Oct01  0.6   RLC   Started using SplineZero as the zero solver
  15Nov01  0.7   RLC   Changed input variable from tableSize to denCode
  24Nov01  0.8   RLC   Fixed coding of optional variables
  27Nov01  0.85  RLC   Reversed order of arguments on MeanLine2
  05Dec01  0.9   RLC   Replace MakeFinalAirfoil with calls to NACAauxilary
  10Dec01  0.91  RLC   Worked on making output look nice
  22Dec01  0.92  RLC   More rearrangement of output
   4Jan02  1.0   RLC   Final cleanup for release of PDAS 7
  07Nov10  1.1   RLC   Fixed typo in GetRk1
*/

#include "splprocs.cxx"
#include "epspsi.cxx"
#include "nacax.cxx"

int
main(int argc, const char **argv)
{
}
