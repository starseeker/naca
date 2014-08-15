*+
      PROGRAM NACA4
*   --------------------------------------------------------------------
*     PURPOSE - To calculate the ordinates and surface slope for both
*            symmetrical and cambered airfoils in the NACA 4-digit,
*            5-digit, and 16-series airfoil families.
*
*     AUTHORS - Charles L.Ladson and Cuyler W. Brooks, NASA Langley
*               Liam Hardy, NASA Ames
*               Ralph Carmichael, Public Domain Aeronautical Software
*
*     REVISION HISTORY
*   DATE  VERS PERSON  STATEMENT OF CHANGES
*   Nov75  1.0 CLL&CWB Original publication in NASA TM X-3284
*    ?     1.1   LH    NAMELIST input; replaced multiple GOTOs with
*                         IF-THEN-ELSE coding; PROFILE interface
*                         Some other modifications.
*  2Dec94  1.2   RLC   Temporarily disabled PROFILE interface
*                         Input from a named file
*  5Dec94  1.3   RLC   Removed the logic that seemed to provide an
*                         alternate approach to computing d1. From Ames??
*  6Dec94  1.4   RLC   Added FIG file for plotting; IMPLICIT NONE
* 24Apr95  1.5   RLC   Renames E to EPS and made it a parameter
*                      Saved scaled output for later 80 column printout
* 28Apr95  1.6   RLC   Restored ABS on tests of CMB for printout. These
*                         were removed at Ames for some reason
*                      Changed namelist name to INPUT4 to avoid 
*                         confusion with datasets for the 6-series 
*                         airfoil program.
*                      Added 'NONE' as a camber option.
*                      Recoded the equation for radius of curvature.
*  8Aug95  1.7   RLC   Additional items on output
*
*     NOTES - This program has been known by the names ANALIN and FOURDIGIT
*
*     REFERENCES-  NASA Technical Memorandum TM X-3284 (November, 1975)
*        "Development of a Computer Program to Obtain Ordinates for
*         NACA 4-Digit, 4-Digit Modified, 5-Digit and 16-Digit Airfoils"
*        by Charles L. Ladson and Cuyler W. Brooks, Jr.
*        NASA Langley Research Center
*
*        "Theory of Wing Sections" by Ira Abbott and Albert Von Doenhoff
*
*     BUG LIST -
*        If DX is too small, you can exceed 200 points and exceed the
*           bounds of the arrays. Not much defense in general against
*           faulty input.
      IMPLICIT NONE
*
************************************************************************
*     C O N S T A N T S                                                *
************************************************************************
      REAL PI,EPS,BIG
      INTEGER LUNRD  ! logical unit number of the input file
      INTEGER LUNOUT ! logical unit number of the output file
      INTEGER LUNFIG ! logical unit number of the plot file
      CHARACTER GREETING*60, AUTHOR*60, VERSION*30, FAREWELL*60
      CHARACTER MODIFIER*60
      PARAMETER (LUNRD=1, LUNOUT=2, LUNFIG=3)
      PARAMETER (GREETING=' NACA4 - Compute airfoil ordinates')
      PARAMETER (AUTHOR=
     &      ' Charles Ladson & Cuyler Brooks, NASA Langley')
      PARAMETER (MODIFIER=' ')  ! put your name here if you make changes
      PARAMETER (VERSION = '1.71 (14Aug95)' )
      PARAMETER (FAREWELL=
     & ' Files naca4.out and naca4.plt added to your directory.')
      PARAMETER (PI=3.14159265, EPS=1.0E-6, BIG=1.0/EPS)
************************************************************************
*     V A R I A B L E S                                                *
************************************************************************
      REAL A0,A1,A2,A3  ! thickness coeff. for forward portion
      REAL AMXL,AMXL1,AMXL2   !
      REAL AS   ! sum of A terms for composite camber
      CHARACTER CAMBER*10   ! type of camber line; valid values are
*                             2-DIGIT, 3-DIGIT, 3DIGITREF,
*                             6-SERIES, 6A-SERIES and NONE
      REAL CHD   ! Model chord used for listing ordinates in
*                         dimensional units.
      REAL CLIS  ! sum of CLI terms for composite camber
      REAL CM    ! Location of maximum camber position (e.g., 0.4).
*                         Note: CMB=K1 and CM=r for 5-digit airfoils, where
*                               K1 and r are defined in NASA TM X-3284.
      REAL CMB   ! Maximum camber ordinate-to-chord ratio (e.g., .04).
      REAL CMBNMR !  number of loadings for 6-series camber line
      REAL D0,D2,D3     ! thickness coeff. for rearward portion
      REAL D1    ! Trailing-edge slope of 4-digit modified airfoils.
      REAL DX    ! Basic chordwise increment in x/c for computing
*                         ordinates.  Usually set to 0.01.
      INTEGER errCode
      REAL F     ! 1 + tan(theta)**2
      REAL G,H   ! sub-expressions in eq. 4-27 of Abbott & Von Doenhoff
      CHARACTER fileName*80   ! holds the name of the input file
      REAL FRAC      ! size of the next step in x
      INTEGER I   ! counts the number of points along each surface
      INTEGER J,L
      INTEGER ICKY   ! integer value of CMBNMR
      REAL K2OK1 ! Defined in NASA TM X-3284.
      INTEGER KON    ! not a clue
*
      CHARACTER NAME*60     ! title on printed and plotted output
      REAL OMXL,OMXL1,OMXL2 ! sub-expressions in eq. 4-27 of A&V
      REAL PIA
      CHARACTER PROFILE*10  ! type of profile; valid values are
*                             4-DIGIT and 4-DIGITMOD
      REAL Q    ! sub-expression in eq. 4-27 of Abbott & Von Doenhoff
      REAL RC        ! radius of curvature at XM. Printed but not used
      REAL RLE   ! Leading-edge radius. Not used with 4-digit but
*                         may be used with 4-digit modified.
      REAL SINTH,COSTH  ! sine,cosine of theta
      REAL THP
      REAL XLL
      REAL XSV
      REAL TOC   ! Thickness-chord ratio of airfoil.
      REAL U,V
      REAL X,Y,YP,YPP
      REAL XC,YC
      REAL XM    ! Nondimensional chordwise location of maximum thick-
*                         ness.  Used for 4-digit modified airfoils only.
      REAL XUC,YUC,XLC,YLC
      REAL YUP,YLP
      REAL Z           ! term in [ ] in eq. 4-27 of Abbot & Von Doenhoff
      REAL Z1,Z2  ! sub-expressions in eq. 4-27 of Abbott & Von Doenhoff
************************************************************************
*     A R R A Y S                                                      *
************************************************************************
      REAL A(10)   ! chordwise position of end of const. load region
      REAL CLI(10) ! CL of const portion of load
      INTEGER IF6XA(10)
      REAL TANTH0(10), YCMB(10), TANTH(10), YCP2(10)
      REAL XU(200), XL(200), YU(200), YL(200)
!!!      REAL COEFFS(4)
      REAL XAU(200),YAU(200),XAL(200),YAL(200)  ! never used??
      REAL YUPR(200),YLPR(200)                  ! never used??
      REAL xScaledU(200),yScaledU(200),xScaledL(200),yScaledL(200)
!!!      COMPLEX ROOTS(3),TEMP(8)      ! apparently never referenced
      EQUIVALENCE (CLI(1),CMB) 
************************************************************************
*     N A M E L I S T   D E F I N I T I O N S                          *
************************************************************************
      NAMELIST /INPUT4/ NAME, PROFILE, CAMBER, TOC, RLE, DX, CHD, XM,
     &                  D1, CMB, CM, K2OK1, CLI, A, CMBNMR

C  NAMELIST "INPUT4":
C
C   VAR     DIM   TYPE    DESCRIPTION
C  NAME      -    C*80    Title desired on printed and plotted output
C  PROFILE   -    C*10    NACA airfoil family (4-digit)
C  CAMBER    -    C*10    Camber line (2-digit or 3-digit; upper case;
C                         left justified.  E.g. '4-DIGIT   '
C  TOC       -     R      Thickness-chord ratio of airfoil.
C  RLE       -     R      Leading-edge radius. Not used with 4-digit but
C                         may be used with 4-digit modified.
C  DX        -     R      Basic chordwise increment in x/c for computing
C                         ordinates.  Usually set to 0.01.
C  CHD       -     R      Model chord used for listing ordinates in
C                         dimensional units.
C  XM        -     R      Nondimensional chordwise location of maximum thick-
C                         ness.  Used for 4-digit modified airfoils only.
C  D1        -     R      Trailing-edge slope of 4-digit modified airfoils.
C  CMB       -     R      Maximum camber ordinate-to-chord ratio (e.g., .04).
C  CM        -     R      Location of maximum camber position (e.g., 0.4).
C                         Note: CMB=K1 and CM=r for 5-digit airfoils, where
C                               K1 and r are defined in NASA TM X-3284.
C  K2OK1     -     R      Defined in NASA TM X-3284.
C
*  CLI      10     R
*  A        10     R
*  CMBNMR    -     R      really an integer; number of camber lines
************************************************************************
*     C O M M O N   B L O C K   D E F I N I T I O N S                  *
************************************************************************
!!!      COMMON I,XU,XL,YU,YL,YUPR(200),YLPR(200)
!!!      COMMON /MAIN/ YSTART(3),CHD,KON
! not sure why this these are defined
************************************************************************
*     D A T A   I N I T I A L I Z A T I O N                            *
************************************************************************
!!!   DATA YSTART /1.0, 4.0, 7.0/
      DATA NAME/'NACAxxxx'/
      DATA PROFILE/'4-DIGIT'/
      DATA CAMBER/'NONE'/
      DATA TOC/0.2/
      DATA RLE/0.02/
      DATA XM/0.5/
!!!      DATA CMB/0.0/
      DATA D1/1.0/
      DATA DX/0.01/
      DATA CHD/1.0/
      DATA CMBNMR/1.0/
      DATA K2OK1/0.0/
      DATA CM/0.5/
      DATA CLI/10*0.0/
      DATA A/10*1.0/
*-----------------------------------------------------------------------
*
      WRITE(*,*) GREETING
      WRITE(*,*) AUTHOR 
      IF (MODIFIER.NE.' ') WRITE(*,*) 'Modified by '//MODIFIER
      WRITE(*,*) 'Version '//VERSION
*
    5 WRITE (*,*) 'Enter the name of the input file: '
      READ  (*, '(A)') fileName
      IF (fileName .EQ. ' ') STOP
      OPEN(UNIT=LUNRD, FILE=fileName, STATUS='OLD', IOSTAT=errCode)
      IF (errCode .NE. 0) THEN
        OPEN(UNIT=LUNRD, 
     &         FILE=fileName//'.inp', STATUS='OLD', IOSTAT=errCode)
        IF (errCode .NE. 0) THEN
          WRITE(*,*) 'Unable to open this file. Try again.'
          GO TO 5
        END IF
      END IF

      KON=0                        ! can't see any possible use for this
*
*
   20 CONTINUE   ! not really needed
!!! if using Watcom Fortran77, change NML to FMT in the READ statements
      READ (LUNRD, NML=INPUT4, END=999)
*..... Set up the files for print and plot output ......................

      OPEN(UNIT=LUNOUT, FILE='naca4.out', STATUS='UNKNOWN')
      OPEN(UNIT=LUNFIG, FILE='naca4.plt', STATUS='UNKNOWN')
      WRITE(LUNFIG,*) '#Created by NACA4 from '//fileName
*      You may want to delete the line above (depends on how you plot)
*
      WRITE(LUNOUT,*) 'NACA4 - A PROGRAM TO COMPUTE NACA '//
     & '4-DIGIT AND 4-DIGIT MODIFIED SECTIONS'
      WRITE(LUNOUT,*) '   INCLUDES 5-DIGIT AND 16-SERIES AIRFOILS'
      WRITE(LUNOUT,*)
      WRITE(LUNOUT,*) 'SUMMARY OF INPUT DATA AND DEFAULTS'
      WRITE(LUNOUT,*)
      WRITE(LUNOUT,NML=INPUT4)                  ! prints all input data
      WRITE(LUNOUT,5000)   !--------------------------------------------
*
      IF (CHD .NE. 1.0) THEN
        WRITE(LUNOUT,*) '   In addition to non-dimensional output, '
        WRITE(LUNOUT, '(A,F7.3)' ) 
     &   '    an airfoil will be printed with chord length = ', CHD
        WRITE(LUNOUT,5000)   !------------------------------------------
      ENDIF
*
      KON=KON+1
      ICKY=INT(CMBNMR)       ! real to integer;  CMBNMR never used again
      FRAC=1.0
      IF (CAMBER.EQ.'6-SERIES  '  .OR.  CAMBER.EQ.'6A-SERIES ' )
     &   CMB=CLI(1)                         ! already equivalenced ???

      DO 50 I=1,10
         IF6XA(I)=0
   50 CONTINUE
      IF (PROFILE.EQ.'4-DIGIT   '   ) THEN
        WRITE(LUNOUT,*) 'THICKNESS PARAMETERS FOR 4-DIGIT AIRFOIL' 
        WRITE(LUNOUT, '(A,F6.3)' ) '   t/c = ', TOC
      ELSEIF (PROFILE.EQ.'4-DIGITMOD') THEN
        WRITE(LUNOUT,*) 
     &    'THICKNESS PARAMETERS FOR 4-DIGIT MODIFIED AIRFOIL' 
        WRITE(LUNOUT, '(A,F9.6)' ) '   t/c         = ', TOC
        WRITE(LUNOUT, '(A,F9.6)' ) '   l.e. radius = ', RLE
        WRITE(LUNOUT, '(A,F9.6)' ) '   t.e. slope  = ', D1
      ELSE
        WRITE(LUNOUT,*) 'INVALID ENTRY FOR PROFILE'
        STOP 'Invalid entry for PROFILE'
      ENDIF  
      WRITE(LUNOUT,5000)   !--------------------------------------------
*
      WRITE(LUNOUT,*) 'CAMBER DEFINITION'
      IF (CAMBER.EQ.'NONE') THEN
        WRITE(LUNOUT,*) 'Uncambered airfoil'
      ELSEIF (CAMBER.EQ.'2-DIGIT') THEN
        WRITE(LUNOUT,*) '2-DIGIT CAMBER LINE PARAMETERS'
        WRITE(LUNOUT,'(A,F8.4)' ) '   max ordinate (y/c)  = ', CMB
        WRITE(LUNOUT,'(A,F8.4)' ) '   x/c of max ordinate = ', CM
      ELSEIF (CAMBER.EQ.'3-DIGIT' ) THEN
        WRITE(LUNOUT,*) '3-DIGIT CAMBER LINE PARAMETERS'
        WRITE(LUNOUT,'(A,F8.4)' ) '   k1 parameter          = ', CMB
        WRITE(LUNOUT,'(A,F8.4)' ) '   x/c of zero curvature = ', CM
      ELSEIF (CAMBER.EQ.'3-DIGITREF') THEN
        WRITE(LUNOUT,*) '3-DIGIT-REFLEX CAMBER LINE PARAMETERS'
        WRITE(LUNOUT,'(A,F8.4)' ) '   k1 parameter          = ', CMB
        WRITE(LUNOUT,'(A,F8.4)' ) '   x/c of zero curvature = ', CM
        WRITE(LUNOUT,'(A,F8.4)' ) '   aft/forward curvature = ', K2OK1
      ELSEIF (CAMBER.EQ.'6-SERIES' .OR.CAMBER.EQ.'6A-SERIES' ) THEN
        WRITE(LUNOUT,*) '6-SERIES CAMBER LINE PARAMETERS'
        WRITE(LUNOUT,*) '  #    CLI    A'
        WRITE(LUNOUT, '(I3,2F9.5)' )(i,CLI(i),A(i),i=1,ICKY)
      ELSE
        WRITE(6,*) 'INVALID ENTRY FOR CAMBER'
        STOP 'Invalid input for CAMBER'
      ENDIF
      WRITE(LUNOUT,5000)   !--------------------------------------------
*
*..... Compute the profile constants if this is a 4-digit modified section
      IF (PROFILE.NE.'4-DIGIT   '   ) THEN
         A0=SQRT(2.0*RLE)*0.2/TOC
         D0=0.002
         IF (D1.EQ.0.0) D1=.1*(2.24-5.42*XM+12.3*XM**2)/(1.-0.878*XM)
         D3=(3.*D1-0.588/(1.-XM))/(3.*(1.-XM)**2)
         D2=-1.5*(1.-XM)*D3-.5*D1/(1.-XM)
         A3=0.1/XM**3+(2.*D1*(1.-XM)-0.588)/
     &      (2.*XM*(1.-XM)**2)-3.*A0/(8.*XM**2.5)
         A2=-0.10/XM**2+.5*A0/XM**1.5-2.*XM*A3
         A1=-.5*A0/XM**.5-2.*XM*A2-3.*XM**2*A3
C     RC IS RADIUS OF CURVATURE AT X=XM   (printed but not used)
!!!         RC=((1.-XM)**2/(2.*D1*(1.-XM)-0.588))*.2/TOC ! radius of curv.
         rc=0.5/ABS(d2+3.0*d3*(1.0-xm))    ! changed by RLC 28Apr95
*
         WRITE(LUNOUT,*) 'DERIVED COEFFICIENTS FOR THE AIRFOIL'
         WRITE(LUNOUT,*) '        0           1          2          3'
         WRITE(LUNOUT,'(4X,A,4F10.6)') 'a',A0,A1,A2,A3, 'd',D0,D1,D2,D3
         !!!WRITE(LUNOUT,790) A0,A1,A2,A3,D0,D1,D2,D3,RC
         !!!790 FORMAT (10H A0,1,2,3=,4F13.6/10H D0,1,2,3=,4F13.6/4H RC=,F13.3//)
         WRITE(LUNOUT,'(A,F9.3)') 'Leading Edge Radius of Curvature=',RC
         WRITE(LUNOUT,5000)   !-----------------------------------------
      ENDIF
C
C     PROFILE, X LE XM
*
*..... Print the headers for the colimns of data
      IF (CAMBER .EQ. 'NONE' .OR. ABS(CMB).LT.EPS) THEN
        WRITE(LUNOUT,*) 
     &    ' DIMENSIONS OF UNCAMBERED AIRFOIL, CHORD LENGTH=1'
        WRITE(LUNOUT, '(//A,6X,A,5X,A,8X,A)' ) 
     &    '   X/C', 'Y/C', 'dY/dX', 'd2Y/dX2'
      ELSE
        WRITE(LUNOUT,*) 
     &    ' DIMENSIONS OF CAMBERED AIRFOIL, CHORD LENGTH=1'
        WRITE(LUNOUT,*) ' UNCAMBERED ---- UPPER SURFACE VALUES ---'
     &    //'    --- LOWER SURFACE VALUES ----'
        WRITE(LUNOUT, '(A,8X,A,7X,A,5X,A,6X,A,7X,A,5X,A)' ) 
     &    '   X/C', 'XU/C', 'YU/C', 'dYU/dXU', 
     &             'XL/C', 'YL/C', 'dYL/dXL'
      ENDIF

!!!  800 FORMAT (//,2X,'X/C',5X,'Y/C',8X,'DY/DX',6X,'D2Y/DX2')
  !!!810 FORMAT (//,' UNCAMBERED', 15X, 'UPPER SURFACE VALUES', 12X,
!!!     & 'LOWER SURFACE VALUES'/ 5X, 'X/C', 9X, 'XU/C', 6X, 'YU/C',5X,
!!!     & 'DYU/DXU', 8X, 'XL/C', 4X,'YL/C', 5X, 'DYL/DXL')

*
      X=0.0
      Y=0.0
      XC=0.0
      YC=0.0
      XU(1)=0.0
      YU(1)=0.0
      XL(1)=0.0
      YL(1)=0.0
      XUC=0.0
      YUC=0.0
      XLC=0.0
      YLC=0.0
      XAU(1)=0.0
      YAU(1)=0.0
      XAL(1)=0.0
      YAL(1)=0.0
!!!!      K=2      ! in original program; never used
      IF (CAMBER.EQ.'NONE') GOTO 110
      IF (CAMBER.EQ.'2-DIGIT') GO TO 120
      IF (CAMBER.EQ.'3-DIGIT') GO TO 130
      IF (CAMBER.EQ.'3-DIGITREF') GO TO 140
      IF (CAMBER.EQ.'6-SERIES') GO TO 150
      IF (CAMBER.EQ.'6A-SERIES') GO TO 150
*
*..... Uncambered airfoil ..............................................
  110 TANTH0(1)=EPS
      YP=BIG
      YPP=BIG
      YUP=-1.0/TANTH0(1)
      YLP=-1.0/TANTH0(1)
      GOTO 230
*
*..... Two digit camber line (leading edge) ............................
  120 TANTH0(1)=2.*CMB/CM
      IF (CMB .LT. EPS) TANTH0(1)=EPS
      YP=BIG
      YPP=BIG
      YUP=-1/TANTH0(1)
      YLP=-1/TANTH0(1)
      GO TO 230
*
*..... Three digit camber line (leading edge) ..........................
  130 TANTH0(1)=CMB*CM**2*(3.0-CM)/6.0
      IF (CMB .LT. EPS) TANTH0(1)=EPS
      YP=BIG
      YPP=BIG
      YUP=-1/TANTH0(1)
      YLP=-1/TANTH0(1)
      GO TO 230
*
*..... Three digit reflex (leading edge) ...............................
  140 TANTH0(1)=CMB*(3.*CM**2-K2OK1*(1-CM)**3-CM**3)/6
      IF (CMB .LT. EPS) TANTH0(1)=EPS
      YP=BIG
      YPP=BIG
      YUP=-1/TANTH0(1)
      YLP=-1/TANTH0(1)
      GO TO 230
*
*..... Six series camber line (leading edge) ...........................
  150 L=0
      CLIS=CLI(1)
      AS=A(1)
  160 L=L+1   ! keep coming back here until L=ICKY
      A(1)=A(L)
      CLI(1)=CLI(L)
!!!      K=2   ! never used
      U=0.005                   ! don't compute at x=0
      V=-(A(1)-U)/ABS(A(1)-U)
      OMXL=(1.-U)*ALOG(1.-U)
      AMXL=(A(1)-U)*ALOG(ABS(A(1)-U))
      OMXL1=-ALOG(1.-U)-1.
      AMXL1=-ALOG(ABS(A(1)-U))+V
      OMXL2=1./(1.-U)
      AMXL2=-V/ABS(A(1)-U)
      IF (A(1) .LT. EPS  .OR.  ABS(1.-A(1)) .LT. EPS) GO TO 170
      G=-(A(1)**2*(.5*ALOG(A(1))-0.25)+0.25)/(1.-A(1))
      Q=1.0
      H=(0.5*(1.-A(1))**2*ALOG(1.-A(1))-0.25*(1.-A(1))**2)/(1.-A(1))+G
      Z=.5*(A(1)-U)*AMXL-.5*(1.-U)*OMXL-.25*(A(1)-U)**2+.25*(1.-U)**2
      Z1=.5*((A(1)-U)*AMXL1-AMXL-(1.-U)*OMXL1+OMXL+(A(1)-U)-(1.-U))
      Z2=.5*(A(1)-U)*AMXL2-AMXL1-.5*(1.-U)*OMXL2+OMXL1
      GOTO 200              ! added by RLC 25Apr95.  Bug in original???
  170 CONTINUE
      IF (A(1) .LT. EPS) GO TO 180
      IF (ABS(A(1)-1.) .LT. EPS) GO TO 190
  180 H=-.5                                        ! special case of a=0
      Q=1.0
      Z1=U*ALOG(U)-.5*U-.5*(1.-U)*OMXL1+.5*OMXL-.5
      GO TO 200
  190 H=0.0                                        ! special case of a=1
      Q=H
      Z1=-OMXL1
      GO TO 200
  200 TANTH0(L)=CLI(1)*(Z1/(1.-Q*A(1))-1.-ALOG(U)-H)/PI/(A(1)+1.)/2.0
      IF (ICKY.GT.1.AND.L.LT.ICKY) GO TO 160  ! go back if more loadings
*
      DO 210 J=2,ICKY                ! old test for ICKY=1 removed   RLC
  210 TANTH0(1)=TANTH0(1)+TANTH0(J)
*
      IF (ABS(CMB) .LT. EPS) TANTH0(1)=EPS       ! don't understand this
      YP=1.E10
      YPP=1.E10
      YUP=-1/TANTH0(1)
      YLP=-1/TANTH0(1)
  230 CONTINUE
*
*..... Print the leading edge point ....................................
      I=1
      IF (CAMBER .EQ. 'NONE' .OR. ABS(CMB).LT.EPS) THEN
         WRITE(LUNOUT,1006) X,Y,YP,YPP
         xScaledU(I)=XC
         yScaledU(I)=YC
      ELSE
         WRITE(LUNOUT,1005) X,XU(I),YU(I),YUP, XL(I),YL(I),YLP
         xScaledU(I)=XUC
         yScaledU(I)=YUC
         xScaledL(I)=XLC
         yScaledL(I)=YLC
      ENDIF
*
*..... All of the logic of the previous 100 lines or so is just to do
*.     the leading edge point.
*..... Now proceed with the computation of the remainder of the airfoil.
*.
      X=.00025
  240 CONTINUE
*     upstream of the maximum thickness point
      IF (PROFILE.EQ.'4-DIGIT') GO TO 250
      IF (PROFILE.EQ.'4-DIGITMOD') GO TO 260
  250 Y=0.29690*SQRT(X)-0.12600*X-0.35160*X**2+0.28430*X**3-0.1015*X**4
      YP=.5*.2969/SQRT(X)-.126-2.*.3516*X+3.*.2843*X*X-4.*.1015*X**3
      YPP=-.5*.5*.2969/SQRT(X**3)-2.*.3516+2.*3.*.2843*X-3.*4.*.1015*X*X
      GO TO 270
  260 Y=A0*X**.5+A1*X+A2*X**2+A3*X**3
      YP=.5*A0/X**.5+A1+2.*A2*X+3.*A3*X**2
      YPP=-.25*A0/X**1.5+2.*A2+6.*A3*X
  270 CONTINUE
      Y=Y*TOC/.2
      YP=YP*TOC/.2
      YPP=YPP*TOC/.2
      XC=X*CHD
      YC=Y*CHD
      IF(ABS(CMB) .LT. EPS) CM=0.5
      IF (CAMBER.EQ.'NONE') GOTO 275
      IF (CAMBER.EQ.'2-DIGIT') GO TO 280
      IF (CAMBER.EQ.'3-DIGIT') GO TO 290
      IF (CAMBER.EQ.'3-DIGITREF') GO TO 300
      IF (CAMBER.EQ.'6-SERIES') GO TO 310
      IF (CAMBER.EQ.'6A-SERIES') GO TO 310
*
*..... Uncambered airfoil ..............................................
  275 YCMB(1)=0.0
      TANTH(1)=0.0
      F=1.0
      THP=0.0
      GOTO 480
*
*..... Two digit camber line ...........................................
  280 YCMB(1)=CMB*(2.0*CM*X-X*X)/CM**2                         ! was 240
      TANTH(1)=2.*CMB*(1.-X/CM)/CM
      IF (X.GT.CM) YCMB(1)=CMB*(1.-2.*CM+2.*CM*X-X*X)/(1.-CM)**2
      IF (X.GT.CM) TANTH(1)=(2.*CM-2.*X)*CMB/(1.-CM)**2
      F=SQRT(1.+TANTH(1)**2)
      THP=-2.*CMB/CM**2/F**2
      IF (X.GT.CM) THP=-2.*CMB/(1.-CM)**2/F**2
      GO TO 480
*
*..... Three digit camber line .........................................
  290 YCMB(1)=CMB*(X**3-3.*CM*X**2+CM**2*(3.-CM)*X)/6.         ! was 250
      TANTH(1)=CMB*(3.*X**2-6.*CM*X+CM**2*(3.-CM))/6.
      IF (X.GT.CM) YCMB(1)=CMB*CM**3*(1.-X)/6.
      IF (X.GT.CM) TANTH(1)=-CMB*CM**3/6.
      F=SQRT(1.+TANTH(1)**2)
      THP=CMB*(X-CM)/F**2
      IF (X.GT.CM) THP=0.0
      GO TO 480
*
*..... Three digit reflexed camber line ................................
  300 YCMB(1)=CMB*((X-CM)**3-K2OK1*(1-CM)**3*X-CM**3*X+CM**3)/6
      TANTH(1)=CMB*(3.*(X-CM)**2-K2OK1*(1-CM)**3-CM**3)/6.
      IF (X.GT.CM) YCMB(1)=CMB*(K2OK1*(X-CM)**3-K2OK1*(1-CM)**3*X-CM**3*
     &  X+CM**3)/6
      IF (X.GT.CM)
     &     TANTH(1)=CMB*(3*K2OK1*(X-CM)**2-K2OK1*(1-CM)**3-CM**3)/6
      F=SQRT(1.+TANTH(1)**2)
      THP=CMB*(X-CM)/F**2
      IF (X.GT.CM) THP=K2OK1*CMB*(X-CM)/F**2
      GO TO 480
*
*..... 6, 6A series camber line ........................................
  310 L=0              ! jump back here from #590 (aft of max thickness)
      A(1)=AS
      CLI(1)=CLIS
  320 L=L+1            ! jump back here from statement 450
      CLI(1)=CLI(L)
      XC=X*CHD
      YC=Y*CHD
      XLL=X*ALOG(X)
      Q=1.0
!..... Equations that follow could divide by zero or have non-positive
!      arguments for the logarithm function. Treat each case............
      IF (ABS(1.-A(1)) .LT. EPS  .AND.  ABS(1.-X).LT.EPS) GO TO 370
      IF (A(1) .LT. EPS  .AND.  (1.-X) .LT. EPS) GO TO 380
      IF (ABS(A(1)-X) .LT. EPS) GO TO 330
      IF (ABS(1.-X) .LT. EPS) GO TO 350
      IF (ABS(A(1)-1.) .LT. EPS) GO TO 360
      V=-(A(1)-X)/ABS(A(1)-X)
      OMXL=(1.-X)*ALOG(1.-X)
      AMXL=(A(1)-X)*ALOG(ABS(A(1)-X))
      OMXL1=-ALOG(1.-X)-1.
      AMXL1=-ALOG(ABS(A(1)-X))+V
      OMXL2=1./(1.-X)
      AMXL2=1./(A(1)-X)
      Z=.5*(A(1)-X)*AMXL-.5*(1.-X)*OMXL-.25*(A(1)-X)**2+.25*(1.-X)**2
      Z1=.5*((A(1)-X)*AMXL1-AMXL-(1.-X)*OMXL1+OMXL+(A(1)-X)-(1.-X))
      Z2=.5*(A(1)-X)*AMXL2-AMXL1-.5*(1.-X)*OMXL2+OMXL1
      IF (A(1) .LE. EPS) GO TO 340   ! got lost at Ames
      G=-(A(1)*A(1)*(.5*ALOG(A(1))-0.25)+0.25)/(1.-A(1))
      H=(0.5*(1.-A(1))**2*ALOG(1.-A(1))-0.25*(1.-A(1))**2)/(1.-A(1))+G
      GO TO 390
*
  330 Z=-.5*(1.-X)**2*ALOG(1.-X)+0.25*(1.-X)**2                    ! x=a
      Z1=-.5*(1.-X)*(-ALOG(1.-X)-1.)+.5*(1.-X)*ALOG(1.-X)-.5*(1.-X)
      Z2=-ALOG(1.-X)-0.5
      G=-(A(1)**2*(.5*ALOG(A(1))-0.25)+0.25)/(1.-A(1))
      H=(0.5*(1.-A(1))**2*ALOG(1.-A(1))-0.25*(1.-A(1))**2)/(1.-A(1))+G
      GO TO 390
  340 G=-.25                                                       ! a=0
      GO TO 390
  350 CONTINUE                                                     ! x=1
      G=-(A(1)**2*(.5*ALOG(A(1))-0.25)+0.25)/(1.-A(1))
      H=(0.5*(1.-A(1))**2*ALOG(1.-A(1))-0.25*(1.-A(1))**2)/(1.-A(1))+G
      Z=.5*(A(1)-1.)**2*ALOG(ABS(A(1)-1.))-0.25*(A(1)-1.)**2
      Z1=-(A(1)-1.)*ALOG(ABS(A(1)-1.))
      Z2=-1.E10
      GO TO 390
  360 G=0.0                                             ! a=1  ! was 320
      H=G
      Q=G
      Z=-(1.-X)*ALOG(1.-X)
      Z1=ALOG(1.-X)+1.
      Z2=-1./(1.-X)
      GO TO 390
  370 Z=0.0                                                ! a=1 and x=1
      G=Z
      H=Z
      Q=Z
      Z1=-1.E10
      Z2=-1.E10
      GO TO 390
  380 G=-.25                                               ! a=0 and x=1
      H=-.5
      Q=1.0
      Z=-.25
      Z1=0.0
      Z2=-10.E10
      GO TO 390
!..... Now compute the camber ordinate. Eq. 4-27 in A & V.
!      Also compute slope (TANTH) at this point
  390 YCMB(L)=CLI(1)*(Z/(1.-Q*A(1))-XLL+G-H*X)/PI/(A(1)+1.)/2.
      XSV=X
      IF (X.LT.0.005) X=0.005
      TANTH(L)=CLI(1)*(Z1/(1.-Q*A(1))-1.-ALOG(X)-H)/PI/(A(1)+1.)/2.0
      X=XSV
      IF (IF6XA(L).EQ.1) TANTH(L)=-5.
      IF (X.GT.0.005) GO TO 400
      YCP2(L)=0.0                                         ! leading edge
      GO TO 420
  400 CONTINUE
      IF (ABS(1.-X) .GT. EPS) GO TO 410
      YCP2(L)=1./EPS                                     ! trailing edge
      GO TO 420
  410 PIA=PI*(A(1)+1.)*2.
      YCP2(L)=CLI(1)*(Z2/(1.-Q*A(1))-1./X)/PIA
  420 CONTINUE
C
C       MODIFIED CAMBERLINE OPTION
C
      IF (CAMBER.EQ.'6A-SERIES') GO TO 430
      GO TO 450
  430 YCMB(L)=YCMB(L)*0.97948
      TANTH(L)=TANTH(L)*0.97948
      YCP2(L)=YCP2(L)*0.97948
      IF (ABS(A(1)-0.8) .LT. EPS) GO TO 440
      WRITE(LUNOUT,840)
      IF (KON.EQ.3) KON=0
      GO TO 20                                       ! error if a <> 0.8
*
  440 CONTINUE
      IF (TANTH(L).LE.-.24521*CLI(1)) YCMB(L)=0.24521*CLI(1)*(1.-X)
      IF (TANTH(L).LE.-.24521*CLI(1)) YCP2(L)=0.0
      IF (TANTH(L).LE.-.24521*CLI(1)) TANTH(L)=-0.24521*CLI(1)
      IF (TANTH(L).LE.-.24521*CLI(1)) IF6XA(L)=1
*
  450 CONTINUE
      IF (ICKY.GT.1.AND.L.LT.ICKY) GO TO 320             ! more loadings
*
      DO 460 J=2,ICKY
      YCMB(1)=YCMB(1)+YCMB(J)
      TANTH(1)=TANTH(1)+TANTH(J)
      YCP2(1)=YCP2(1)+YCP2(J)
  460 CONTINUE
*
      F=SQRT(1.+TANTH(1)**2)
      THP=YCP2(1)/F**2
  480 CONTINUE
      IF(X.GT.XM) GO TO 600        ! if we got to #310 by jump from #590
      IF(ABS(X-XM).LT.EPS) GO TO 600   ! then go back there
  490 CONTINUE
      SINTH=TANTH(1)/F
      COSTH=1./F
      I=I+1
      XU(I)=X-Y*SINTH
      YU(I)=YCMB(1)+Y*COSTH
      XL(I)=X+Y*SINTH
      YL(I)=YCMB(1)-Y*COSTH
      XAU(I) = XU(I)
      YAU(I) = YU(I)
      XAL(I) = XL(I)
      YAL(I) = YL(I)
  500 CONTINUE
      XUC=XU(I)*CHD
      YUC=YU(I)*CHD
      XLC=XL(I)*CHD
      YLC=YL(I)*CHD
      IF (ABS(CMB) .LE. EPS) GO TO 510
      YUP=0.0
      YLP=YUP
      IF (ABS(TANTH(1)) .LT. EPS) GO TO 510
      YUP=(TANTH(1)*F+YP-TANTH(1)*Y*THP)/(F-YP*TANTH(1)-Y*THP)
      YLP=(TANTH(1)*F-YP+TANTH(1)*Y*THP)/(F+YP*TANTH(1)+Y*THP)
      YUPR(I)=YUP            ! never used??
      YLPR(I)=YLP            ! never used??
  510 CONTINUE                                          ! was 460
*
*..... Print for this value of x (ahead of max thickness) ..............
      IF (CAMBER .EQ. 'NONE' .OR. ABS(CMB).LT.EPS) THEN
         WRITE(LUNOUT,1006) X,Y,YP,YPP
         xScaledU(I)=XC
         yScaledU(I)=YC
      ELSE
         WRITE(LUNOUT,1005) X,XU(I),YU(I),YUP, XL(I),YL(I),YLP
         xScaledU(I)=XUC
         yScaledU(I)=YUC
         xScaledL(I)=XLC
         yScaledL(I)=YLC
      ENDIF
*
*..... Set up the next value of X ......................................
      IF (X.LE.0.0975) FRAC=0.25
      IF (X .LE. 0.0124) FRAC=0.025
C     IF (X.LE.0.00225) FRAC=0.025
      X=X+FRAC*DX
      FRAC=1.0
      IF(ABS(X-XM) .LT. EPS) GO TO 520
      IF(X.LT.XM) GO TO 240           ! if still fore of max thick, jump
C
C     PROFILE - X GE XM
C
      X=XM                         ! insures that XM is a computed point
*
*..... downstream of the maximum thickness point
  520 CONTINUE
      IF (PROFILE.EQ.'4-DIGIT   ') GO TO 530
      IF (PROFILE.EQ.'4-DIGITMOD') GO TO 540
  530 Y=0.2969*SQRT(X)-0.126*X-0.3516*X**2+0.2843*X**3-0.1015*X**4
      YP=.5*.2969/SQRT(X)-.126-2.*.3516*X+3.*.2843*X*X-4.*.1015*X**3
      YPP=-.5*.5*.2969/SQRT(X**3)-2.*.3516+2.*3.*.2843*X-3.*4.*.1015*X*X
      GO TO 550
  540 Y=D0+D1*(1.-X)+D2*(1.-X)**2+D3*(1.-X)**3
      YP=-D1-2.*D2*(1.-X)-3.*D3*(1.-X)**2
      YPP=2.0*D2 + 6.0*D3*(1.-X)
  550 CONTINUE
      Y=Y*TOC/.2
      YP=YP*TOC/.2
      YPP=YPP*TOC/.2
      XC=X*CHD
      YC=Y*CHD
      IF (CAMBER.EQ.'NONE      ') GO TO 555
      IF (CAMBER.EQ.'2-DIGIT   ') GO TO 560
      IF (CAMBER.EQ.'3-DIGIT   ') GO TO 570
      IF (CAMBER.EQ.'3-DIGITREF') GO TO 580
      IF (CAMBER.EQ.'6-SERIES  ') GO TO 590
      IF (CAMBER.EQ.'6A-SERIES ') GO TO 590
*..... Uncambered airfoil ..............................................
  555 YCMB(1)=0.0
      TANTH(1)=0.0
      F=1.0
      THP=0.0
      GOTO 610

*
*..... Two digit camber line ...........................................
  560 YCMB(1)=CMB*(2.0*CM*X-X*X)/CM**2
      TANTH(1)=2.*CMB*(1.-X/CM)/CM
      IF (X.GT.CM) YCMB(1)=CMB*(1.-2.*CM+2.*CM*X-X*X)/(1.-CM)**2
      IF (X.GT.CM) TANTH(1)=(2.*CM-2.*X)*CMB/(1.-CM)**2
      F=SQRT(1.+TANTH(1)**2)
      THP=-2.*CMB/CM**2/F**2
      IF (X.GT.CM) THP=-2.*CMB/(1.-CM)**2/F**2
      GO TO 610
*..... Three digit camber line .........................................
  570 YCMB(1)=CMB*(X**3-3.*CM*X**2+CM**2*(3.-CM)*X)/6.
      TANTH(1)=CMB*(3.*X**2-6.*CM*X+CM**2*(3.-CM))/6.
      IF (X.GT.CM) YCMB(1)=CMB*CM**3*(1.-X)/6.
      IF (X.GT.CM) TANTH(1)=-CMB*CM**3/6.
      F=SQRT(1.+TANTH(1)**2)
      THP=CMB*(X-CM)/F**2
      IF (X.GT.CM) THP=0.0
      GO TO 610
*..... Three digit reflexed camber line ................................
  580 YCMB(1)=CMB*((X-CM)**3-K2OK1*(1-CM)**3*X-CM**3*X+CM**3)/6
      TANTH(1)=CMB*(3.*(X-CM)**2-K2OK1*(1-CM)**3-CM**3)/6.
      IF (X.GT.CM) YCMB(1)=CMB*(K2OK1*(X-CM)**3-K2OK1*(1-CM)**3*X-CM**3*
     &   X+CM**3)/6
      IF (X.GT.CM) TANTH(1)=
     &    CMB*(3*K2OK1*(X-CM)**2-K2OK1*(1-CM)**3-CM**3)/6
      F=SQRT(1.+TANTH(1)**2)
      THP=CMB*(X-CM)/F**2
      IF (X.GT.CM) THP=K2OK1*CMB*(X-CM)/F**2
      GO TO 610
*..... Six series camber line ..........................................
  590 GO TO 310                        ! big jump to use previous coding
*
  600 CONTINUE
  610 CONTINUE
      SINTH=TANTH(1)/F
      COSTH=1./F
      I=I+1
      XU(I)=X-Y*SINTH
      YU(I)=YCMB(1)+Y*COSTH
      XL(I)=X+Y*SINTH
      YL(I)=YCMB(1)-Y*COSTH
*
  620 CONTINUE
      XUC=XU(I)*CHD
      YUC=YU(I)*CHD
      XLC=XL(I)*CHD
      YLC=YL(I)*CHD
      IF (CAMBER .EQ. 'NONE') GOTO 630
      IF (ABS(CMB) .LE. EPS) GO TO 630
      YUP=0.0
      YLP=YUP
      IF (ABS(TANTH(1)).LT.1.E-10) GO TO 630
      YUP=TANTH(1)*(F+YP/TANTH(1)-Y*THP)/(F-YP*TANTH(1)-Y*THP)
      YLP=TANTH(1)*(F-YP/TANTH(1)+Y*THP)/(F+YP*TANTH(1)+Y*THP)
      YUPR(I)=YUP    ! never used??
      YLPR(I)=YLP    ! never used??
*
  630 CONTINUE
*
*..... Print for this value of X (aft of max thickness) ................
      IF (CAMBER .EQ. 'NONE' .OR. ABS(CMB).LT.EPS) THEN
         WRITE(LUNOUT,1006) X,Y,YP,YPP
         xScaledU(I)=XC
         yScaledU(I)=YC
      ELSE
         WRITE(LUNOUT,1005)
     &      X,XU(I),YU(I),YUP, XL(I),YL(I),YLP
         xScaledU(I)=XUC
         yScaledU(I)=YUC
         xScaledL(I)=XLC
         yScaledL(I)=YLC
      ENDIF
*
      X=X+DX
      XAU(I)=XU(I)
      YAU(I)=YU(I)
      XAL(I)=XL(I)
      YAL(I)=YL(I)
      IF (X .LE. 1.0) GO TO 520    ! keep looping until complete
*
      WRITE(LUNOUT,690) NAME
  632 CONTINUE
C      DO 633 II=1, I, 4
C      I1=II + 3
C      WRITE (LUN,870) (XAU(K),YAU(K),K=II,I1)
  633 CONTINUE
C      DO 634 II=1, I, 4
C      I1=II + 3
C      WRITE (LUN,870) (XAL(K),YAL(K),K=II,I1)
  634 CONTINUE
*
      WRITE(LUNOUT,5000)   !--------------------------------------------
      IF (CHD .NE. 1.0) THEN
        WRITE(LUNOUT, '(4X,A,F9.4)') 
     &    'AIRFOIL SCALED TO CHORD LENGTH=', CHD
        IF (CAMBER .EQ. 'NONE' .OR. ABS(CMB).LT.EPS) THEN
          WRITE(LUNOUT,*) '      x          y'
          WRITE(LUNOUT, '(2F11.5)') (xScaledU(J),yScaledU(J), J=1,I)
        ELSE
          WRITE(LUNOUT,*) '    x-upper    y-upper    x-lower    y-lower'
          WRITE(LUNOUT, '(4F11.5)' )
     &         (xScaledU(J),yScaledU(J), xScaledL(J),yScaledL(J), J=1,I)
        ENDIF
        WRITE(LUNOUT,5000)   !------------------------------------------
      ENDIF
*
*..... Write the section data to the plotting file .....................
!!!      CALL PlotSection(LUNFIG, I, XU,YU, XL,YL)
      WRITE(LUNFIG, '(2F12.6)' ) (XL(J),YL(J),J=I,1,-1)       ! te to le
      WRITE(LUNFIG, '(2F12.6)' ) (XU(J),YU(J),J=2,I)        ! le+1 to te
*
!!!  635 GO TO 20   ! big jump back
C
  690 FORMAT (/10X,A,10X,A/)
  700 FORMAT (//,' PROFILE PARAMETERS FOR FOUR DIGIT PROFILE'/
     & ' T/C=', F10.5/' L.E.RADIUS=',F10.5/
     & ' BASIC X INTERVAL=', F10.5/' CHORD=',F10.5/)
  710 FORMAT (//,' PROFILE PARAMETERS FOR MODIFIED FOUR DIGIT PROFILE'/
     & ' T/C=',F10.5/' L.E.RADIUS=',F10.5/
     & ' BASIC X INTERVAL=', F10.5/' CHORD=',F10.5/
     & ' POSITION OF MAXIMUM THICKNESS, XM=', F10.5/
     & ' CONSTANT D1=',F10.5/)
  720 FORMAT (' CAMBER LINE PARAMETERS, 2-DIGIT'/
     & ' CAMBER(YCMAX) =', F10.5/' POSITION OF MAXIMUM CAMBER=',F10.5/)
  730 FORMAT (' CAMBER LINE PARAMETERS, 3-DIGIT'/
     & ' CAMBER PARAMETER K1=',F10.5/
     & ' POSITION OF ZERO CAMBER LINE CURVATURE=', F10.5/)
  740 FORMAT (' CAMBER LINE PARAMETERS, 3-DIGIT-REFLEX'/
     & ' CAMBER PARAMETER K1=', F10.5/
     & ' POSITION OF ZERO CAMBER LINE CURVATURE=', F10.5/
     & ' RATIO OF AFT TO FORWARD CAMBER LINE CURVATURE FACTOR, K2OK1=',
     & F10.5/)
  750 FORMAT (' CAMBER LINE PARAMETERS, 6-SERIES'/7X,3HCLI,9X,1HA)
  760 FORMAT (2F10.3)
  !!!770 FORMAT (3H J=,I2/6H IERR=,I2/4H D1=,E13.6)
  !!! FORMAT (14H NO REAL ROOTS)
  790 FORMAT (10H A0,1,2,3=,4F13.6/10H D0,1,2,3=,4F13.6/4H RC=,F13.3//)
  800 FORMAT (//,2X,'X/C',5X,'Y/C',8X,'DY/DX',6X,'D2Y/DX2')
  810 FORMAT (//,' UNCAMBERED', 15X, 'UPPER SURFACE VALUES', 12X,
     & 'LOWER SURFACE VALUES'/ 5X, 'X/C', 9X, 'XU/C', 6X, 'YU/C',5X,
     & 'DYU/DXU', 8X, 'XL/C', 4X,'YL/C', 5X, 'DYL/DXL')
  830 FORMAT (2F8.5,2F12.5,8X,2F12.5)   ! uncambered
  840 FORMAT ('  MODIFIED CAMBER LINE OPTION ONLY ALLOWED IF A=0.8')
  850 FORMAT (F10.5,10X,4F10.5,E11.2,6X,4F10.5,E11.2) ! cambered
  870 FORMAT (8F10.6)
!   use 1005 for cambered, 1006 for uncambered print statements
 1005 FORMAT(F8.5,2F11.5,F11.4,2F11.5,F11.4)
 1006 FORMAT(0PF8.5,F9.5,G12.4,1PE12.3)
 
*      CALL PRWRIT ( LUNSAV, 200, NAME, I, I, XU, XL, YU, YL, 1 )
      CLOSE(UNIT=LUNRD)
      CLOSE(UNIT=LUNOUT)
      CLOSE(UNIT=LUNFIG)
*
      WRITE (*, '(A)' ) FAREWELL
      STOP 'NACA4 has terminated successfully.'

  999 WRITE(*,*) 'Unable to find a valid namelist record.'
      WRITE(*,*) 'No output produced'
      STOP 'No valid data'
*
 5000 FORMAT('-----------------------------------------------------'///)
      END ! ---------------------------------- End of main program NACA4
 
