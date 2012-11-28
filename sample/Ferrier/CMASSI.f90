  module CMASSI_mod
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   implicit none
!
!-----------------------------------------------------------------------
      REAL, PARAMETER :: DMImin=.05e-3, DMImax=1.e-3,      &
     &  XMImin=1.e6*DMImin, XMImax=1.e6*DMImax
      INTEGER, PARAMETER :: MDImin=XMImin, MDImax=XMImax
!-----------------------------------------------------------------------
!--- Mean mass of precpitation ice particles as functions of their mean
!    size (in microns)
!
      REAL MASSI(MDImin:MDImax)
!
!--- Mean rain drop diameters vary from 50 microns to 450 microns
!
      REAL, PARAMETER :: DMRmin=.05E-3, DMRmax=.45E-3, DelDMR=1.E-6    &
     &, XMRmin=1.E6*DMRmin, XMRmax=1.E6*DMRmax, N0r0=8.E6, N0rmin=1.e4
      INTEGER, PARAMETER :: MDRmin=XMRmin, MDRmax=XMRmax
!
!--- Various rain lookup tables
!
      REAL MASSR(MDRmin:MDRmax),RQR_DRmin,RQR_DRmax,           &
           CN0r0,CN0r_DMRmin,CN0r_DMRmax
!
  end module  CMASSI_mod
