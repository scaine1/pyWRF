      subroutine dbzcalc(qvp,qra,qsn,qgr,tmk,prs,dbz,miy,mjx,mkzh,
     &   in0r,in0s,in0g,iliqskin)

c    this routine expects to be passed the following.
c    qvp = water vapor mixing ratio
c    qra = rain mixing ratio
c    qsn = snow mixing ratio
c    qgr = graupel mixing ratio
c    tmk = temperature in kelvin (I think)
c    prs = pressure
c    dbz = input dbz array to store your dbz data (zeros)
c    miy = horizontal dimenions (not sure which one)
c    mjx = horizontal dimenion (not sure which one)
c    mkzh= vertical dimension
c    in0r = can be 0 or 1, 1 means use variable intercepts, 0 means use the constants defined below as rn0_x
c    in0s = 0 or 1, dont use 1 for the Lin scheme
c    in0g = 0 or 1, even though it mentions thompson below, this scheme is not consistant with his assumptions! DO NOT USE WITH THOMPSON
c    iliqskin = 1, simulate bright-band effect for frozen particles above 0 degrees (assume they have liquid water on the outside)



c !set up dimensions for our arrays
      dimension qra(miy,mjx,mkzh),qvp(miy,mjx,mkzh),            
     &  qsn(miy,mjx,mkzh),qgr(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &  prs(miy,mjx,mkzh),dbz(miy,mjx,mkzh)



c The line below is required to get python to use the fortran code 

Cf2py intent(in,out) dbz(miy,mjx,mkzh)  

      character meth*1



c      include 'comconst'
c
c   Constant intercepts
c

c     in the lin scheme the particle-size drop distrbutions are modelled by a marshall-palmer distribution of the
c     form N(x) = rn0_x exp (-lambda_x D)   where D is the drop diameter, lambda is the slope of the distrbujtion and
c     lambda is defined as ((pi rn0_x rho_x) / (rho_air* q_x))**(1/4)  where q_x is the mixing ratio

c     rn0_r= intercept parameter for rain 
c     rn0_s= intercept parameter for snow 
c     rn0_g= intercept parameter for graupel

      rn0_r = 8.e6    ! m^-4
c      rn0_s = 2.e7    ! m^-4   simon changed for old lin
      rn0_s = 3.e6
      rn0_g = 4.e6    ! m^-4
c
c   Constants used to calculate variable intercepts
c
      rhowat=1000.                              !density of water
      pi=4.*atan(1.)                            !pi
      rgas=287.04                               !gas constant
      r1=1.e-15                                 !minimum rain threshold value (guessing) 
      ron=8.e6                                          
      ron2=1.e10
      son=2.e7
      gon=5.e7
      ron_min = 8.e6
      ron_qr0 = 0.00010
      ron_delqr0 = 0.25*ron_qr0
      ron_const1r = (ron2-ron_min)*0.5
      ron_const2r = (ron2+ron_min)*0.5
c
c   Other constants:
c

c    For a particle-drop size distrbution represented by a marshall-palmer distribution
c    simulated radar reflectivty can be represented as the 6th moment of the drop size
c    distribution which works out to be (Z_ex = equivelent reflectivity factor for x)
c    Z_ex = gamma(7) rn0_x  lambda_x**-7
c    where lambda_x = ((pi rn0_x rho_x) / (rho_air* q_x))**(1/4) 

c    substituing lambda_x into Z_e gives
c    Z_ex = (720(rho_air q_x)**(7/4))/ ( nr0_x**(3/4) ( pi rho_x) **(7/4))



      gamma_seven = 720.                !gamma(7)                 !
      rho_r = rhowat ! density of water  1000. kg m^-3
      rho_s = 100.   ! density of snow in kg m^-3
      rho_g = 400.   ! density of graupel in kg m^-3
      alpha = 0.224  ! this factor is used to take into account that ice scatters the radar beam differently to water

      factor_r = gamma_seven * 1.e18 * (1./(pi*rho_r))**1.75   ! 1e18 is to get it to convert to standard reflectivity units.
      factor_s = gamma_seven * 1.e18 * (1./(pi*rho_s))**1.75   ! these factors will bs used later, they are just calculated now
     &    * (rho_s/rhowat)**2 * alpha                          ! this is only half the relfectivty equation
      factor_g = gamma_seven * 1.e18 * (1./(pi*rho_g))**1.75
     &    * (rho_g/rhowat)**2 * alpha
c



      do k=1,mkzh
c      do j=1,mjx-1
c      do i=1,miy-1
      do j=1,mjx
      do i=1,miy

c
c         rhoair=prs(i,j,k)*100./
c     &      (rgas*virtual(tmk(i,j,k),qvp(i,j,k)))      ! air density

c########################SIMON CHANCE RHOAIR#############################
       
      rhoair=prs(i,j,k)*100.0/(rgas*(tmk(i,j,k)*                !calculate rho air, this should actually be dry air (fix)
     &    (0.622+qvp(i,j,k))/(0.622*(1.+qvp(i,j,k)))))

c
c      Adjust factor for brightband, where snow or graupel particle
c      scatters like liquid water (alpha=1.0) because it is assumed to
c      have a liquid skin.
c
         if (iliqskin.eq.1.and.tmk(i,j,k).gt.celkel) then
            factorb_s=factor_s/alpha
            factorb_g=factor_g/alpha
         else
            factorb_s=factor_s
            factorb_g=factor_g
         endif


c
c      Calculate variable intercept parameters if want ed  c
c     this is stupid, you can calculate the variable as in Thompson et al, but he uses totally different assumptions so you cannot actually apply this algorithm to his data!
 
c  SKIP TO LINE 172 and below         
     
         if (in0s.eq.1) then  ! N_0s as in Thompson et al. 
            temp_c = amin1(-0.001, tmk(i,j,k)-celkel)
            sonv = amin1(2.0e8, 2.0e6*exp(-0.12*temp_c))
         else
            sonv = rn0_s
         endif



c
         if (in0g.eq.1) then  ! N_0g as in Thompson et al.
            gonv = gon
            if (qgr(i,j,k).gt.r1) then
               gonv = 2.38*(pi*rho_g/
     +            (rhoair*qgr(i,j,k)))**0.92
               gonv = max(1.e4, min(gonv,gon))
            endif
         else
            gonv = rn0_g
         endif


c
         if (in0r.eq.1) then  ! N_0r as in Thompson et al.
            ronv = ron2
            if (qra(i,j,k).gt. r1) then
                ronv = ron_const1r*tanh((ron_qr0-qra(i,j,k))
     +             /ron_delqr0) + ron_const2r
            endif
c            if (i.eq.97.and.j.eq.100) then
c               print*,'k,ron2,qra(i,j,k),r1,ronv,rn0_r='
c               print*,'    ',k,ron2,qra(i,j,k),r1,ronv,rn0_r
c            endif
         else
            ronv = rn0_r
         endif


c
c      Total equivalent reflectivity factor (z_e, in mm^6 m^-3) is
c      the sum of z_e for each hydrometeor species:
c
         z_e =   factor_r  * (rhoair*qra(i,j,k))**1.75 /     ! rain    !this is the second half of the relfectivity equation, using factor_x defined above
     &           ronv**.75                                   ! ronv = rn0_r = intercept parameter
     &         + factorb_s * (rhoair*qsn(i,j,k))**1.75 /     ! snow
     &           sonv**.75
     &         + factorb_g * (rhoair*qgr(i,j,k))**1.75 /     ! graupel
     &           gonv**.75

c          print * , factor_r
c          print * , (rhoair*qra(i,j,k))**1.75      ! rain
c          print * , ronv**0.75




c
c      Adjust small values of Z_e so that dBZ is no lower than -20
c
         z_e = max(z_e,.01)


c
c      Convert to dBZ
c
         dbz(i,j,k) = 10. * log10(z_e)    !convert from dB to dBZ


c
      enddo
      enddo
      enddo
c
      return 
      end
