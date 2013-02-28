      SUBROUTINE CALMICT(P1D,T1D,Q1D,C1D,FI1D,FR1D,FS1D,CUREFL,     &
                        QW1,QI1,QR1,QS1,DBZ1,DBZR1,DBZI1,DBZC1,     &
                        NLICE1,NSDIM,WEDIM,MASSR,MASSI,im,jm)
!
!
!
!To compile for python:
!
!       f2py -c -m --f77flags="-fconvert=big-endian" --f90flags="-fconvert=big-endian" 
!       dbzcalc_ferrier_py dbzcalc_ferrier_py.f90 --fcompiler=gfortran
!
!If successfull will create a *.so file which python can call.
!
!
!
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    CALMIC      COMPUTES HYDROMETEORS 
!   PRGRMMR: JIN         ORG: W/NP2      DATE: 01-08-14       
!     
! ABSTRACT:  
!     THIS ROUTINE COMPUTES THE MIXING RATIOS OF CLOUD WATER,
!     CLOUD ICE, RAIN, AND SNOW.  THE CODE IS BASED ON SUBROUTINES
!     GSMDRIVE & GSMCOLUMN IN THE NMM MODEL. 
!     
! PROGRAM HISTORY LOG:
!   01-08-14  YI JIN 
!   02-02-11  Brad Ferrier - Minor changes for consistency w/ NMM model
!   04-11-10  Brad Ferrier - Removed cloud fraction algorithm
!   04-11-17  H CHUANG - WRF VERSION     
! USAGE:    CALL CALMICT(P1D,T1D,Q1D,C1D,FI1D,FR1D,FS1D,CUREFL
!     &,                QW1,QI1,QR1,QS1,DBZ1,DBZR1,DBZI1,DBZC1)
!   INPUT ARGUMENT LIST:
!     P1D     - PRESSURE (PA)
!     T1D     - TEMPERATURE (K)
!     Q1D     - SPECIFIC HUMIDITY (KG/KG)
!     C1D     - TOTAL CONDENSATE (CWM, KG/KG)
!     FI1D    - F_ice (fraction of condensate in form of ice)
!     FR1D    - F_rain (fraction of liquid water in form of rain)
!     FS1D    - F_RimeF ("Rime Factor", ratio of total ice growth 
!                       to deposition growth)
!     CUREFL  - Radar reflectivity contribution from convection (mm**6/m**3)
!               NOTE: when the grid size is small cumulus parameterisation is
!               usually turned off and so this will be zero
!
!   OUTPUT ARGUMENT LIST:
!     QW1   - CLOUD WATER MIXING RATIO (KG/KG)
!     QI1   - CLOUD ICE MIXING RATIO (KG/KG)
!     QR1   - RAIN MIXING RATIO (KG/KG)
!     QS1   - "SNOW" (precipitation ice) MIXING RATIO (KG/KG)
!     DBZ1  - Equivalent radar reflectivity factor in dBZ; i.e., 10*LOG10(Z)
!     DBZR  - Equivalent radar reflectivity factor from rain in dBZ
!     DBZI  - Equivalent radar reflectivity factor from ice (all forms) in dBZ
!     DBZC  - Equivalent radar reflectivity factor from parameterized convection in dBZ
!
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     FUNCTIONS:
!        FPVS
!     UTILITIES:
!     LIBRARY:
!       NONE
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : IBM SP
!$$$  
!
!BSF: skip      use params_mod

      !use cmassi_mod
! start CMASSI hardcoded
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
      !LW REAL MASSI(MDImin:MDImax)
!
!--- Mean rain drop diameters vary from 50 microns to 450 microns
!
      REAL, PARAMETER :: DMRmin=.05E-3, DMRmax=.45E-3, DelDMR=1.E-6    &
     &, XMRmin=1.E6*DMRmin, XMRmax=1.E6*DMRmax, N0r0=8.E6, N0rmin=1.e4
      INTEGER, PARAMETER :: MDRmin=XMRmin, MDRmax=XMRmax
!
!--- Various rain lookup tables
!
      REAL RQR_DRmin,RQR_DRmax,           &
           CN0r0,CN0r_DMRmin,CN0r_DMRmax

      !REAL MASSR(MDRmin:MDRmax),RQR_DRmin,RQR_DRmax,           &
      !     CN0r0,CN0r_DMRmin,CN0r_DMRmax

!end CMASSI hardcoded


! hardcode in from params.f
      real, parameter :: D608=0.608
      real, parameter :: TFRZ=273.15
      real, parameter :: EPSQ=1.E-12
      real, parameter :: T_ICE=-30.
      real, parameter :: PI=3.141592653589793

!BSF: need to define array limits for im, jm
!
!     jm = WEdim and im=NSdim --> in pyWRF miy = WEdim and mjx=NSdim
      integer, intent(in) :: nsdim, wedim 
      integer, intent(in) :: im,jm
      
      INTEGER INDEXS, INDEXR
      REAL, PARAMETER :: Cice=1.634e13, DBZmin=-20., RD=287.04,   &
                         ONEPS=1.-18.015/28.964

!BSF: adjust following parameters based on whats in phys/module_mp_etanew.F:
      real, parameter :: NLImin=1.E3, NLImax=5.E3


      !real,dimension(nsdim,IM),intent(in) :: P1D,T1D,Q1D,C1D,FI1D,FR1D,     &
      !
      real,dimension(IM,JM),intent(in) :: P1D,T1D,Q1D,C1D,FI1D,FR1D,     &
           FS1D,CUREFL
      !real,dimension(nsdim,IM),intent(inout) ::  QW1,QI1,QR1,QS1,DBZ1,DBZR1,&
      !
      real,dimension(IM,JM),intent(inout) ::  QW1,QI1,QR1,QS1,DBZ1,DBZR1,&
           DBZI1,DBZC1,NLICE1
      
      REAL N0r,Ztot,Zrain,Zice,Zconv,Zmin
      integer I,J,JSTA,JEND


      real TC, Frain,Fice,Flimass,Flarge,     &
           Fsmall,RimeF,Xsimass,Qice,Qsat,ESAT,WV,RHO,RRHO,RQR,          &
           DRmm,Qsigrd,WVQW,Dum,XLi,Qlice,WC,DLI,xlimass
      
! From MICROINIT
      REAL, PARAMETER :: RHOL=1000.
      real ax,C_N0r0
      REAL, DIMENSION(MDRmin:MDRmax),intent(in) :: MASSR
      REAL, DIMENSION(MDImin:MDImax),intent(in) :: MASSI

      !integer i
! End MICROINIT 

!f2py intent(in) P1D
!f2py intent(in) T1D
!f2py intent(in) Q1D
!f2py intent(in) C1D
!f2py intent(in) FI1D
!f2py intent(in) FR1D
!f2py intent(in) FS1D
!f2py intent(in) CUREFL
!f2py intent(in) QW1
!f2py intent(in) QI1
!f2py intent(in) QR1
!f2py intent(in) QS1
!f2py intent(in,out) DBZ1
!f2py intent(in) DBZR1
!f2py intent(in) DBZI1
!f2py intent(in) DBZC1
!f2py intent(in) NLICE1
!f2py intent(in) JM
!f2py intent(in) IM
!f2py intent(in) MASSR
!f2py intent(in) MASSI

!BSF: skip      real,external :: fpvs
!************************************************************************
!--- Determine composition of condensate in the form of cloud water, 
!    total ice (cloud ice & snow), & rain following GSMDRIVE in NMM model
!


! Start MICROINIT hardcoded
!------------------------ START EXECUTION ------------------------
!
!---  READ IN MASSI FROM LOOKUP TABLES 
!
!      OPEN (UNIT=1,FILE="ETAMPNEW_DATA",FORM="UNFORMATTED",       &
!     &status='old')
!      DO I=1,3
!        READ(1)
!      ENDDO
!      READ(1) MASSR
!      DO I=1,5
!        READ(1)
!      ENDDO
!      READ(1) MASSI
!      CLOSE(1)
      RQR_DRmin=N0r0*MASSR(MDRmin)    ! Rain content for mean drop diameter of .05 mm
      RQR_DRmax=N0r0*MASSR(MDRmax)    ! Rain content for mean drop diameter of .45 mm
!      PI=ACOS(-1.)
      C_N0r0=PI*RHOL*N0r0
      CN0r0=1.E6/C_N0r0**.25
      CN0r_DMRmin=1./(PI*RHOL*DMRmin**4)
      CN0r_DMRmax=1./(PI*RHOL*DMRmax**4)
!      print *,'MICROINIT: MDRmin, MASSR(MDRmin)=',MDRmin,MASSR(MDRmin)
!      print *,'MICROINIT: MDRmax, MASSR(MDRmax)=',MDRmax,MASSR(MDRmax)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Initialise empty arrays - set dBZ arrays = -20 (clear sky)
      JSTA=1
      JEND=JM

      Zmin=10.**(0.1*DBZmin)
      DO J=JSTA,JEND
        !DO I=1,IM
        DO I=1,IM
          !QW1(I,J)=0.  ! Do not initialise these - these are now passed in from the WRF outfile
          !QI1(I,J)=0.
          !QR1(I,J)=0.
          !QS1(I,J)=0.
          NLICE1(I,J)=0.
          DBZ1(I,J)=DBZmin
          DBZR1(I,J)=DBZmin
          DBZI1(I,J)=DBZmin
          DBZC1(I,J)=DBZmin
        ENDDO
      ENDDO
      DO J=JSTA,JEND
        !DO I=1,IM
        DO I=1,IM
          Zrain=0.            !--- Radar reflectivity from rain
          Zice=0.             !--- Radar reflectivity from ice
          Zconv=CUREFL(I,J)   !--- Radar reflectivity from convection
          IF (C1D(I,J) .LE. EPSQ) THEN
!
!--- Skip rest of calculatiions if no condensate is present
!
            GO TO 10
          ELSE
            WC=C1D(I,J)
          ENDIF
!
!    QI1 - total ice (cloud ice & snow) mixing ratio
!    QW1 - cloud water mixing ratio
!    QR1 - rain mixing ratio

!BSF:  If you're running the ARW and it has as output arrays qw, qr, qi, then
!      start here.
!-- array QW is the same as the array 'qc' from ARW output and 'QCLOUD' from wrfout
!-- array QR is the same as the array 'qr' from ARW output and 'QRAIN' from wrfout
!-- array QI is the same as the array 'qs' from ARW output and 'QSNOW' from wrfout
!    => Here QI includes the total ice, large ice (snow) & small cloud ice
!       This is separated out down below in the code starting with defining
!       the constant 'FLARGE'.
!
!BSF: skip          WV=Q1D(I,J)/(1.-Q1D(I,J))
!
!--- Saturation vapor pressure w/r/t water ( >=0C ) or ice ( <0C )
!
!BSF: skip          ESAT=1000.*FPVS(T1D(I,J))
!BSF: skip          QSAT=EPS*ESAT/(P1D(I,J)-ESAT)
      RHO=P1D(I,J)/(RD*T1D(I,J)*(1.+D608*Q1D(I,J)))
          RRHO=1./RHO
  !
  !--- Based on code from GSMCOLUMN in model to determine reflectivity from rain
  !
  !      
  !    Using a Marshall-Palmer distribution:
  !       
  !       Ze = integral of D*6 No e*(-lamda*D)
  !    where 
  !    D is the droplet diameter in m
  !    No is the intercept parameter
  !    lamda is the slope
  !
  !    Assuming that rain droplets are perfectly spherical
  !    and that D<<1:
  !
  !        Ze,r = GAMMA(7) No,r / lamda*7
  !    where
  !    GAMMA(n) is the Euler gamma function and GAMMA(7)=720 (m)
  !    
  !    In the Ferrier scheme D=lamda*-1
  !    and so:
  !
  !        Ze,r = 720 * No,r * lamda*7
  !
  !
 
          IF (QR1(I,J) .GT. EPSQ) THEN
            RQR=RHO*QR1(I,J)
            IF (RQR .LE. RQR_DRmin) THEN
              N0r=MAX(N0rmin, CN0r_DMRmin*RQR)
              INDEXR=MDRmin
            ELSE IF (RQR .GE. RQR_DRmax) THEN
              N0r=CN0r_DMRmax*RQR
              INDEXR=MDRmax
           ELSE
              N0r=N0r0
              INDEXR=MAX( XMRmin, MIN(CN0r0*RQR**.25, XMRmax) )
            ENDIF
  !
  !--- INDEXR is the mean drop size in microns, DR is 'drop radius'; 
  !    convert all to mm:
  !
  !    GAMMA(7) = 720 * 1E-3 = 0.72
  !    No,r     = N0r   (intercept)
  !    lamda    = DRmm  (slope)
  !    
    
            DRmm=1.e-3*REAL(INDEXR)
            Zrain=0.72*N0r*DRmm*DRmm*DRmm*DRmm*DRmm*DRmm*DRmm
          ENDIF        !--- End IF (QR1(I,J) .GT. EPSQ) block
!
!--- Based on code from GSMCOLUMN in model to determine partition of 
!    total ice into cloud ice & snow (precipitation ice)
!
          IF (QI1(I,J) .GT. EPSQ) THEN
            QICE=QI1(I,J)
            RHO=P1D(I,J)/(RD*T1D(I,J)*(1.+ONEPS*Q1D(I,J)))
            RRHO=1./RHO
!BSF: skip            QSIgrd=RHgrd*QSAT
!BSF: skip            WVQW=WV+QW1(I,J)
!
!  * FLARGE  - ratio of number of large ice to total (large & small) ice
!  * FSMALL  - ratio of number of small ice crystals to large ice particles
!  ->  Small ice particles are assumed to have a mean diameter of 50 microns.
!  * XSIMASS - used for calculating small ice mixing ratio
!  * XLIMASS - used for calculating large ice mixing ratio
!  * INDEXS  - mean size of snow to the nearest micron (units of microns)
!  * RimeF   - Rime Factor, which is the mass ratio of total (unrimed &
!              rimed) ice mass to the unrimed ice mass (>=1)
!  * FLIMASS - mass fraction of large ice
!  * QTICE   - time-averaged mixing ratio of total ice
!  * QLICE   - time-averaged mixing ratio of large ice
!  * NLICE1  - time-averaged number concentration of large ice
       
!
!
!
            IF (TC.GE.0.) THEN    !-- BSF: skip   .OR. WVQW.LT.QSIgrd) THEN
              FLARGE=1.
            ELSE
              FLARGE=.1
              IF (TC.GE.-8. .AND. TC.LE.-3.) FLARGE=.5*FLARGE
            ENDIF
            FSMALL=(1.-FLARGE)/FLARGE
            XSIMASS=RRHO*MASSI(MDImin)*FSMALL
            DUM=XMImax*EXP(.0536*TC)
            INDEXS=MIN(MDImax, MAX(MDImin, INT(DUM) ) )
            RimeF=AMAX1(1., FS1D(I,J) )
            XLIMASS=RRHO*RimeF*MASSI(INDEXS)
            FLIMASS=XLIMASS/(XLIMASS+XSIMASS)
            QLICE=FLIMASS*QICE
            NLICE1(I,J)=QLICE/XLIMASS
            IF (NLICE1(I,J).LT.NLImin .OR. NLICE1(I,J).GT.NLImax) THEN
!
!--- Force NLICE1 to be between NLImin and NLImax
!
              DUM=MAX(NLImin, MIN(NLImax, NLICE1(I,J)) )
              XLI=RHO*(QICE/DUM-XSIMASS)/RimeF
              IF (XLI .LE. MASSI(MDImin) ) THEN
                INDEXS=MDImin
              ELSE IF (XLI .LE. MASSI(450) ) THEN
                DLI=9.5885E5*XLI**.42066         ! DLI in microns
                INDEXS=MIN(MDImax, MAX(MDImin, INT(DLI) ) )
              ELSE IF (XLI .LE. MASSI(MDImax) ) THEN
                DLI=3.9751E6*XLI**.49870         ! DLI in microns
                INDEXS=MIN(MDImax, MAX(MDImin, INT(DLI) ) )
              ELSE 
                INDEXS=MDImax
!
!--- 8/22/01: Increase density of large ice if maximum limits
!    are reached for number concentration (NLImax) and mean size
!    (MDImax).  Done to increase fall out of ice.
!
                IF (DUM .GE. NLImax)                              &
                  RimeF=RHO*(QICE/NLImax-XSIMASS)/MASSI(INDEXS)
              ENDIF             ! End IF (XLI .LE. MASSI(MDImin) )
              XLIMASS=RRHO*RimeF*MASSI(INDEXS)
              FLIMASS=XLIMASS/(XLIMASS+XSIMASS)
              QLICE=FLIMASS*QICE
              NLICE1(I,J)=QLICE/XLIMASS
            ENDIF               ! End IF (NLICE.LT.NLImin ...
!LW            QS1(I,J)=AMIN1(QI1(I,J), QLICE)
!LW            QI1(I,J)=AMAX1(0., QI1(I,J)-QS1(I,J))


!       Equation (C.8) in Ferrier (1994, JAS, p. 272)
!       Radar reflectivity of dry ice particles as a function of mass and number concentrations
!       
!       Zi = C'x(RHO*Qi)**2 / ni
!       where
!       C'x     = 0.224 * 10**12 ( 6/PI*RHO )**2 * GAMMA(7-ax) * ( GAMMA(1-ax) / GAMMA(4-ax)**2) 
!               = 1.63410E13 for exponential ice distributions (ax = 0)
!       ni      = NLICE number concetration in cm**-3
!       RHO     = Density in g/cm**3
!       Qi      = QLICE
!

   !--- Equation (C.8) in Ferrier (1994, JAS, p. 272), which when
   !    converted from cgs units to mks units results in the same
   !    value for Cice, which is equal to the {} term below:
   !
   !    Zi={.224*720*(10**18)/[(PI*RHOL)**2]}*(RHO*QLICE)**2/NLICE1(I,J),
   !    where RHOL=1000 kg/m**3 is the density of liquid water
   !
   !--- Valid only for exponential ice distributions
   !
            Zice=Cice*RHO*RHO*QLICE*QLICE/NLICE1(I,J) 
          ENDIF                 ! End IF (QI1(I,J) .GT. 0.) THEN
!
!---  Calculate total (convective + grid-scale) radar reflectivity
10        Ztot=Zrain+Zice+Zconv
          IF (Ztot .GT. Zmin)  DBZ1(I,J)= 10.*ALOG10(Ztot)
          IF (Zrain .GT. Zmin) DBZR1(I,J)=10.*ALOG10(Zrain)
          IF (Zice .GT. Zmin)  DBZI1(I,J)=10.*ALOG10(Zice)
          IF (Zconv .GT. Zmin) DBZC1(I,J)=10.*ALOG10(Zconv)
        ENDDO
      ENDDO
!
      RETURN
      END
