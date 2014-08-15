
COMPUTATION OF ORDINATES OF STANDARD NACA AIRFOILS        \naca456\readme.txt

The files for this program are in the directory \naca456 on the CD-ROM:
  readme.txt      this file
  naca4.txt       instructions and sample cases for naca4
  naca6.txt       instructions and sample cases for naca6
  naca4.for       source code for naca4
  naca6.for       source code for naca6
  naca4.exe       executable for MS-DOS  (F77L)
  naca6.exe       executable for MS-DOS  (F77L)
  nacaplot.for    source code for nacaplot.for
  nacaplot.exe    program to make plots on HP or PostScript printer
  inpdata.zip     input data sets for ALL of the airfoils in Abbott...
  tm4741.pdf      NASA Tech. Memo 4741  (Dec 1996)
  rep452.pdf      NACA Report 452 by T.Theodorsen
  rep824.pdf      NACA Report 824 by Abbott, von Doenhoeff, and Stivers

To use this program, create a directory on your hard disk and copy these 
files to that directory.

During the 1930's the National Advisory Committee for Aeronautics(NACA)
developed several families of airfoils and camber lines. Many of these 
airfoil shapes have been used successfully over the years as wing and tail 
sections for general aviation and military aircraft, as well as propellers 
and helicopter rotors.

The ordinates for numerous specific airfoils of these series are published 
and available. However, when performing parametric studies on effects of
such variables as thickness, location of maximum thickness, leading-
edge radius, location of maximum camber and others, it is not always easy
to obtain the ordinates of the desired shapes. To remedy thisproblem the 
NASA Langley Research Center sponsored the development of computer programs 
for generation of ordinates of standard NACA airfoils.

Two separate programs were written by Charles Ladson and Cuyler Brooks
of the NASA Langley Research Center in 1974-1975. The first, which is 
referred to here as NACA4, was documented in NASA TM X-3284 and produces 
ordinates for NACA 4-digit, 4-digit modified, 5-digit, and 16-series 
airfoils. All of these airfoils use the 4-digit or 4-digit modified 
thickness distributions. These thickness families are defined by algebraic 
equations and they are combined with appropriate mean lines to produce the 
final thick cambered airfoil.

The second program, referred to here as NACA6, was documented in NASA
TM X-3069 and produces ordinates for NACA 6-series and 6A-series airfoils.
Unlike the airfoils of NACA4, these thickness distributions are not
defined by algebraic equations, but use complex variable mapping of a
circle into an airfoil shape.

I obtained both programs at NASA Ames, where the formatted input
had been replaced with a Namelist input scheme and some of the original
coding with multiple GOTOs was replaced with IF-THEN-ELSE logic of Fortran
77. In preparing the programs for this collection, I have made only
minor modifications. The output originally designed for printing on a page 
with 132 columns has been reformatted so that users can view it on a CRT 
screen with 80 columns. I added a number of additional comments to 
the source code. The comments of the original authors are in upper case
and mine are in lower case.

In addition to the formatted print file naca4.prt or naca6.prt, 
a second file called naca4.plt or naca6.plt may be used for plotting or 
for further analysis. This file has a table of (x,y) pairs of points on 
the surface of the airfoil. The points begin at the trailing edge of the 
lower surface and proceed forward along the lower surface around the
leading edge and back along the upper surface to the trailing edge. This
is sometimes referred to as the clockwise way of defining an airfoil.
Many of the airfoils have slightly blunt trailing edges, so the first
point will not necessarily coincide with the last.

Refer to the file naca4.txt for additional description of program
naca4 and to naca4.inp for instructions on preparing input.

Refer to the file naca6.txt for additional description of program
naca6 and to naca6.inp for instructions on preparing input.

Refer to the file namelist.txt if you plan to recompile the program.
You may need to modify the source code or the data, depending upon
which compiler you use.

A program called nacaplot is included.  It reads the plotting file 
naca4.plt or naca6.plt and produce two files called naca.ps and naca.pcl.
The file naca.ps may be sent to a PostScript printer for plotting and 
naca.pcl may be sent to a HP-compatible printer for plotting.

In December 1996, NASA released a new report, TM 4741 that effectively
supercedes the original reports. I have included a copy of this
report in Portable Document Format (PDF) for your convenience. You
will need a PDF viewer, such as Adobe Acrobat Reader to view or print
this file. If you need a free copy of Adobe Acrobat Reader, go to the
Adobe web site at   http://www.adobe.com. I have also included a copy
of the classic papers by Theodorsen, NACA Report 452, and Abbott,
von Doenhoeff and Stivers, NACA Report 824. These two reports contain
most of the material in the classic book by Abbott and von Doenhoeff.
Every student of airfoil theory simply MUST have a copy of "The Theory
of Airfoil Sections".



