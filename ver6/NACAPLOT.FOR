*+
      PROGRAM NACAairfoilPlotter
*   --------------------------------------------------------------------
*     PURPOSE - Make a PostScript file and a PCL file from the output of 
*        either NACA4 or NACA6.
*
*     AUTHOR - Ralph L. Carmichael, Public Domain Aeronautical Software
*
*     REVISION HISTORY
*   DATE  VERS PERSON  STATEMENT OF CHANGES
* 14Aug95  1.0   RLC   Original coding (made a metafile)
* 16Nov95  1.1   RLC   Make Ps and Pcl files directly
*  8Dec95  1.2   RLC   Now entering HP-GL mode properly
* 31Dec96  2.0   RLC   Protect against inability to open files
*
      IMPLICIT NONE
************************************************************************
*     C O N S T A N T S                                                *
************************************************************************
      INTEGER PLT,PS,PCL ! units
      PARAMETER (PLT=1, PS=2, PCL=3)
      CHARACTER GREETING*63, VERSION*30, FAREWELL*63
      CHARACTER AUTHOR*63, MODIFIER*63
      PARAMETER (GREETING='nacaplot - A program to plot airfoil data')
      PARAMETER (AUTHOR=
     & 'Ralph L. Carmichael, Public Domain Aeronautical Software')
      PARAMETER (MODIFIER=' ')
      PARAMETER (VERSION='1.2 (8Dec95)' )
      PARAMETER (FAREWELL=
     & 'Files naca.ps and naca.pcl have been added to your directory.')
************************************************************************
*     V A R I A B L E S                                                *
************************************************************************
      CHARACTER fileName*80
      INTEGER code
*-----------------------------------------------------------------------
      WRITE(*,*) GREETING
      WRITE(*,*) AUTHOR
      WRITE(*,*) 'Modified by '//MODIFIER
      WRITE(*,*) 'Version '//VERSION

    5 WRITE(*,*) 'Enter the name of the input file: '
      READ(*, '(A)' ) fileName
      IF (fileName .EQ. ' ') STOP
      OPEN (UNIT=PLT, FILE=fileName, STATUS='OLD', IOSTAT=code)
      IF (code .NE. 0) THEN
        WRITE (*,*) 'Unable to open this file. Try again'
        GOTO 5
      END IF

      OPEN (UNIT=PS,  FILE='naca.ps',  STATUS='UNKNOWN', IOSTAT=code)
      IF (code .NE. 0) STOP 'Unable to open PostScript file'
      OPEN (UNIT=PCL, FILE='naca.pcl', STATUS='UNKNOWN', IOSTAT=code)  
      IF (code .NE. 0) STOP 'Unable to open PostScript file'

      CALL PostScriptHeader(PS)
      CALL PCLheader(PCL)
      CALL CopyData(PLT, PS, PCL)
      CALL PostScriptTrailer(PS)
      CALL PCLtrailer(PCL)

      CLOSE(UNIT=PLT)
      CLOSE(UNIT=PS)
      CLOSE(UNIT=PCL)
      WRITE (*,*) FAREWELL
      STOP 'nacaplot has terminated successfully.'
      END ! --------------------- End of main program NACAairfoilPlotter
*+
      SUBROUTINE PostScriptHeader(efu)
*   --------------------------------------------------------------------
      IMPLICIT NONE
************************************************************************
*     A R G U M E N T S                                                *
************************************************************************
      INTEGER efu   ! external file unit
*-----------------------------------------------------------------------
      WRITE(efu,*) '%!PS-Adobe'
      WRITE(efu,*) '612 0 translate 90 rotate'
      RETURN
      END ! ------------------------- End of subroutine PostScriptHeader
*+
      SUBROUTINE PclHeader(efu)
*   --------------------------------------------------------------------
      IMPLICIT NONE
************************************************************************
*     A R G U M E N T S                                                *
************************************************************************
      INTEGER efu   ! external file unit
*-----------------------------------------------------------------------
      WRITE(efu, '(3A1,A)' ) CHAR(27),'E',CHAR(27),'%0BINRO90IPSP1;PU'
      RETURN
      END ! -------------------------------- End of subroutine PCLheader
*+
      SUBROUTINE CopyData(infile,psfile,pclfile)
*   --------------------------------------------------------------------
*     PURPOSE - Copy data from infile and make appropriate entries in
*         the PostScript file and the PCL file
*
*     NOTES-
*
      IMPLICIT NONE
************************************************************************
*     A R G U M E N T S                                                *
************************************************************************
      INTEGER infile    ! unit # of file from naca?.exe
      INTEGER psfile    ! unit # of the PostScript file
      INTEGER pclfile   ! unit # of the PCL file
************************************************************************
*     L O C A L   C O N S T A N T S                                    *
************************************************************************
      REAL VMINX,VMAXX, VMINY,VMAXY, XLOW,XHIGH, YLOW,YHIGH
      REAL SCLX,SCLY
      REAL MARGIN,SCALEPS,SCALEPCL
      PARAMETER (VMINX=0.02, VMAXX=0.98, VMINY=0.02, VMAXY=0.7)
      PARAMETER (XLOW=0.0, XHIGH=1.0, YLOW=-0.3125, YHIGH=0.3125)
      PARAMETER (SCLX=(VMAXX-VMINX)/(XHIGH-XLOW))
      PARAMETER (SCLY=(VMAXY-VMINY)/(YHIGH-YLOW))
      PARAMETER (MARGIN=32.0, SCALEPS=720.0, SCALEPCL=10160.0)
************************************************************************
*     L O C A L   V A R I A B L E S                                    *
************************************************************************
      CHARACTER*132 buffer
      INTEGER code
      LOGICAL firstPoint
      REAL x,y
      REAL scx,scy
      REAL psx,psy             ! PostScript coordinates are REAL
      INTEGER pclx,pcly        ! PCL coordinates are INTEGER
*-----------------------------------------------------------------------
      firstPoint = .TRUE.
      DO 
        READ(infile, '(A)', IOSTAT=code) buffer
        IF (code .LT. 0) EXIT
        IF (buffer(1:1) .EQ. '#') CYCLE               ! gnuplot comments
        READ (buffer, *, IOSTAT=code) x,y
        IF (code .NE. 0) CYCLE
        scx=VMINX + (x-XLOW)*SCLX
        scy=VMINY + (y-YLOW)*SCLY
        psx=MARGIN + scx*SCALEPS
        psy=MARGIN + scy*SCALEPS
        pclx= INT(scx*SCALEPCL)
        pcly= INT(scy*SCALEPCL)
        
        IF (firstPoint) THEN
          firstPoint = .FALSE.
          WRITE(psfile,  '(2F9.2, A)' ) psx,psy, ' moveto'
          WRITE(pclfile, '(A,2I7,A)' ) 'PUPA', pclx,pcly, ' PD'
        ELSE
          WRITE(psfile, '(2F9.2, A)' ) psx,psy, ' lineto'
          WRITE(pclfile, '(2I7)' ) pclx,pcly
        END IF
      END DO
      RETURN
      END ! --------------------------------- End of subroutine CopyData
*+
      SUBROUTINE PostScriptTrailer(efu)
*   --------------------------------------------------------------------
      IMPLICIT NONE
************************************************************************
*     A R G U M E N T S                                                *
************************************************************************
      INTEGER efu   ! external file unit
*-----------------------------------------------------------------------
      WRITE(efu,*) 'closepath stroke showpage'
      RETURN
      END ! ------------------------ End of subroutine PostScriptTrailer

*+
      SUBROUTINE PclTrailer(efu)
*   --------------------------------------------------------------------
      IMPLICIT NONE
************************************************************************
*     A R G U M E N T S                                                *
************************************************************************
      INTEGER efu   ! external file unit
*-----------------------------------------------------------------------
      WRITE (efu, '(A1,A3,2A1)' ) CHAR(27), '%0A', CHAR(27), 'E'
      RETURN
      END ! ------------------------------- End of subroutine PCLtrailer

