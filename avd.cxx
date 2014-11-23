/*
!+
! PROGRAM AbbottVonDoenhoff
! ---------------------------------------------------------------------------
! PURPOSE - Make tables of the data in Abbott and von Doenhoff book
!   entitled "Theory of Airfoil Sections".
!   The output will be HTML 5 files for viewing in a browser.
! AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
!
!     REVISION HISTORY
!   DATE  VERS PERSON  STATEMENT OF CHANGES
! 20Nov01  0.1   RLC   Original coding (used program in \atmos as guide )
! 21Nov01  0.2   RLC   Added the actual airfoil calculations
! 23Nov01  0.3   RLC   Cambered sections both interp and non-interp
! 27Nov01  0.4   RLC   Rearranged tables at head of each section
! 30Nov01  0.5   RLC   Changed output to four figures; added leRadius
! 02Dec01  0.6   RLC   Renamed some procedures. Added remainder of profiles
! 07Dec01  0.7   RLC   Finally got leRadius, leSlope right
! 20Jan02  0.8   RLC   Put #top at top of page
! 26Jan02  0.9   RLC   Separated sections into 3 groups
! 07Apr02  1.0   RLC   Added the missing profiles and sections. Noted charset
! 08Feb03  1.1   RLC   Made the resulting files XHTML-1.1 STRICT 
! 16Feb03  1.2   RLC   Fixed items that caused non-compliance with XHTML 1.1
! 15Jan09  1.3   RLC   Additional fixes for XHTML 1.1 and file names are *.html
! 16Jan09  1.4   RLC   Added WriteBanner,WriteCrumb,WriteFooter
! 17Jan09  1.5   RLC   Added MainPage; removed Welcome and debug file 
! 07Nov10  1.6   RLC   Revised to create HTML 5
!
! NOTES-
*/

#include "splprocs.cxx"
#include "epspsi.cxx"
#include "nacax.cxx"
