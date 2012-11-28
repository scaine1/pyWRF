      SUBROUTINE MICROINIT
!
!-- ABSTRACT:
!     Initializes arrays for new cloud microphysics
!
!-- Program History Log:
!     02-02-08  B. Ferrier
!     04-11-19 H CHUANG - WRF VERSION
!
!-- Input argument list:
!     None
!
!-- Output argument list:
!     None
!
!-- Subprograms called:
!     Function FPVS
!
!-- Common blocks:
!     CMASSI
!     RMASS_TABLES
!     MAPOT
!     CRHgrd
!
!-- Attributes:
!     Language: FORTRAN 90
!     Machine : IBM SP
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      use cmassi_mod
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
      REAL, PARAMETER :: RHOL=1000.
      real ax,C_N0r0
      integer i
!
!------------------------ START EXECUTION ------------------------
!
!---  READ IN MASSI FROM LOOKUP TABLES 
!
      OPEN (UNIT=1,FILE="eta_micro_lookup.dat",FORM="UNFORMATTED")
      DO I=1,3
        READ(1)
      ENDDO
      READ(1) MASSR
      DO I=1,5
        READ(1)
      ENDDO
      READ(1) MASSI
      CLOSE(1)
      RQR_DRmin=N0r0*MASSR(MDRmin)    ! Rain content for mean drop diameter of .05 mm
      RQR_DRmax=N0r0*MASSR(MDRmax)    ! Rain content for mean drop diameter of .45 mm
      PI=ACOS(-1.)
      C_N0r0=PI*RHOL*N0r0
      CN0r0=1.E6/C_N0r0**.25
      CN0r_DMRmin=1./(PI*RHOL*DMRmin**4)
      CN0r_DMRmax=1./(PI*RHOL*DMRmax**4)
      print *,'MICROINIT: MDRmin, MASSR(MDRmin)=',MDRmin,MASSR(MDRmin)
      print *,'MICROINIT: MDRmax, MASSR(MDRmax)=',MDRmax,MASSR(MDRmax)
!      print *,  'ETA2P:MASSI(50)= ', MASSI(50)
!      print *,  'ETA2P:MASSI(450)= ', MASSI(450)
!      print *,  'ETA2P:MASSI(1000)= ', MASSI(1000)
!
!--- Initialize saturation vapor pressure lookup tables (functions FPVS, FPVS0)
!
!BSF: skip      CALL GPVS
!
!--- Initialize RHgrd, grid-scale RH for onset of condensation. 
!    See GSMCONST in Eta model for algorithm with grid-size dependence.
!
!      AX=111.*(DPHD**2+DLMD**2)**.5
!      AX=111.*(DYVAL/1000.**2+DXVAL/1000.**2)**.5
!      AX=MIN(100., MAX(5., AX) )
!      RHgrd=0.90+.08*((100.-AX)/95.)**.5
!BSF skip:      RHgrd=1.
!--- 
      RETURN
      END
