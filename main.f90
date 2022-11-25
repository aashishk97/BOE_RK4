module module_check 

implicit none

!qbar = dimensionless
!E0 = dimensionless
!eps0 = dimensionless
!We = dimensionless
!rs = dimensionless
!dq = dimensionless(hbaromd/Ef)
!Wph = not dimensionless(eV)
!D_{gamma} = dimensionless 
!hbaromd =  not dimensionless(eV)
!domega is ((wf-wi)/no of omega steps)
!mM = dimensionless( m/M)
!wla=9.875e13(1/sec), wta=5.0834e13(1/sec)
!omd=9.875 = not dimensionless(om[(LA] is 9.875e13, 10 power is calculated 
!                             and result is in the coefficient of K_ph_ph)
!hbar=1.0545= not dimensionless(again 10 power is calculated and is in the coefficient of K_ph_ph) 
!K_{gamma} = dimensionless
!D_{gamma} = dimensionless
!v_{gamma} = not dimensionless(m/s)
!E2, E3, lambda, mu are elastic coefficcients
!KT = not dimensionless(eV)
!Temp = dimensionless


integer :: i, j, k, t, n, m, p, ip1, ip2, ie, iw, iep, ix, code,nit
integer :: IERR, rank, size, npp, root, istart, iend, lpstart, lpend, lpp, &
        tpstart, tpend, tpp
integer, parameter :: Np=86 
real :: start, finish, par_cal_eph_LA, par_cal_eph_TA1, par_cal_phe_LA, par_cal_phe_TA1, &
        par_cal_eph_TA2, par_cal_phe_TA2, x, y, z, xp, yp, zp, omega_LApp, kph_ph, de ,zval
real, parameter :: qbar=1.173, E0=0.005558, eps0=0.02565, We=0.004275, dt=0.01, &
                   rs = 2.07,  pi = 4*atan(1.0), dq = 0.005558, Wph = 20e-4, domega_LA = 0.0125, &
                   hbaromd = 0.065, Kla = 2.036e-5, Kta1 = 5.39e-6, Dla = 0.66667, Dta1 = 0.246, &
                   mM = 2.0333e-5, Ef = 11.695, B = 0.01, domega_TA1 = 0.012554, &
                   Kta2 = 5.39e-6, Dta2 = 0.246, domega_TA2 = 0.012554, vla = 6471, vta = 3331.11, &
                   E2 = -133.83, E3 = -119.63, lambda = 53.16, mu = 29.96, sigma = 0.013, &
                   q0=1.4835, v0=6656.381, hbar=1.0545,Temp=0.00589397, KT = 0.026, bpar=119.63, &
                   omd=9.875
real, allocatable, dimension(:,:,:) :: Cee, Cph_php1, Cph_phm1, Cph_phm2, Cph_php2, Cph_php3, Cph_phm3, &
                   recv4
real, allocatable, dimension(:,:) :: f, phonon_LA, phonon_TA1, phonon_TA2, phonon_TA1p, phonon_TA2p
real, allocatable, dimension(:) :: eps, el_initial, omega_LA, omega_TA1, ph_initial_LA, ph_initial_TA1,&
                   omega_TA2, ph_initial_TA2, eq_QQ1, eq_QQ2, eq_QQ3, omega_TA1p, omega_TA2p, &
                   dfermi, dphonon_LA, dphonon_TA1, dphonon_TA2, dphonon_LA_temp, &
                   dphonon_TA1_temp, recv, recv2, recv3, f_old, pLA_old, pTA1_old, pTA2_old

contains 

!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine in_cond_el(eps, el_initial) 

implicit none

real, intent(in), dimension(n) :: eps
real, intent(out), dimension(n) :: el_initial

el_initial = 1.0/(exp((eps-1.0)/Temp)+1.0) + 0.1*exp(-((eps-1.0-eps0)/(2.0*We))**2) &
          - 0.1*exp(-((eps-1.0+eps0)/(2.0*We))**2) 

!open(unit = 1, file = 'in_el_dist.dat')
!do ie = 1, n
!  write(1,10) (eps(ie)-1.0)/E0, el_initial(ie)
!end do
!close(1)
!
!10 format(f15.8, 4X, f15.8)

end subroutine in_cond_el
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine in_cond_ph_LA(omega_LA, ph_initial_LA)

implicit none

real, intent(in), dimension(m) :: omega_LA
real, intent(out),dimension(m) :: ph_initial_LA

ph_initial_LA = (1.0/(exp(omega_LA*hbaromd/KT)-1.0)) + B*exp(-(((omega_LA-1)*hbaromd)/(2.0*Wph))**2)

!open(unit = 2, file = 'in_ph_dist.dat')
!do iw = 1, m
!  write(2,10) omega_LA(iw), ph_initial_LA(iw)
!end do
!close(2)
!
!10 format(f15.8, 4X, f15.8)

end subroutine in_cond_ph_LA
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine in_cond_ph_TA1(omega_TA1, ph_initial_TA1)

implicit none

real, dimension(p) :: omega_TA1
real,dimension(p) :: ph_initial_TA1

ph_initial_TA1 = (1.0/(exp(omega_TA1*hbaromd/KT)-1.0))

end subroutine in_cond_ph_TA1
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine in_cond_ph_TA2(omega_TA2, ph_initial_TA2)

implicit none

real, dimension(p) :: omega_TA2
real, dimension(p) :: ph_initial_TA2

ph_initial_TA2 = (1.0/(exp(omega_TA2*hbaromd/KT)-1.0))

end subroutine in_cond_ph_TA2
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
real function coupling_el(smax,smin)

implicit none

real :: smax,smin

coupling_el =((qbar*smax/(qbar**2+smax**2))+atan(smax/qbar)-&
              (qbar*smin/(qbar**2+smin**2))-atan(smin/qbar))

return
end function
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
subroutine deltaf(it)

use mpi
implicit none

integer :: ixp,phflag
integer, intent(in) :: it
real :: r1,xip,fie,fiep,fix,fixp,r2,r3,eixp,rph1, r4, rph2, rph3, r5
real :: r31,r32
real, allocatable :: fixpp(:)

allocate(fixpp(-n:2*n))

do ie=-n,0
  fixpp(ie)=1.0
end do
do ie=1,n
  fixpp(ie)=f(ie,it)
end do
do ie=n+1,2*n
  fixpp(ie)=0.0
end do

phflag=0
recv = 0

r2=de*(1.0/(2.0*(qbar**3)))*4.556
rph1=domega_LA*par_cal_eph_LA
rph2=domega_TA1*par_cal_eph_TA1
rph3=domega_TA2*par_cal_eph_TA2
!write(*,*) rph2, rph3
!stop

do ie= istart, iend
  r1=0.0
  r3=r2/sqrt(eps(ie))
  fie=f(ie,it)
  do iep=1,n
     fiep=f(iep,it)
     do ix=1,n
       fix=f(ix,it)
       ixp=ie+iep-ix
       r1=r1+r3*Cee(ie,iep,ix)*((-fie)*fiep*(1.0-fix)*(1.0-fixpp(ixp))&
                                +(1.0-fie)*(1.0-fiep)*fix*fixpp(ixp))
     end do
   end do

  r3=rph1/sqrt(eps(ie))
  r31=0.0
  do iw=1,m
     xip=eps(ie)-omega_LA(iw)*dq
     ixp=int((xip-eps(1))/dE)
     eixp=1-10.0*E0+float(ixp)*dE
     fixp = fixpp(ixp)+((xip-eixp)*(fixpp(ixp+1)-fixpp(ixp)))/dE
     r1=r1+r3*((omega_LA(iw))**2)* &
          ((fixp-fie)*phonon_LA(iw,it)-fie*(1.0-fixp))

     xip=eps(ie)+omega_LA(iw)*dq
     ixp=int((xip-eps(1))/dE)
     eixp=1-10.0*E0+float(ixp)*dE
     fixp = fixpp(ixp)+((xip-eixp)*(fixpp(ixp+1)-fixpp(ixp)))/dE
     r31=r31+r3*((omega_LA(iw))**2)* &
          ((fixp-fie)*phonon_LA(iw,it)+fixp*(1.0-fie))
  end do
!  write(*,*) 'r31=',r31
  r1=r1+r31

  r4=rph2/sqrt(eps(ie))
  r32=0.0
  do ip1=1,p
     xip=eps(ie)-omega_TA1(ip1)*dq
     ixp=int((xip-eps(1))/dE)
     eixp=1-10.0*E0+float(ixp)*dE
     fixp = fixpp(ixp)+((xip-eixp)*(fixpp(ixp+1)-fixpp(ixp)))/dE
     r32=r32+r4*((omega_TA1(ip1))**2)* &
          ((fixp-fie)*phonon_TA1(ip1,it)-fie*(1.0-fixp))
     
     xip=eps(ie)+omega_TA1(ip1)*dq
     ixp=int((xip-eps(1))/dE)
     eixp=1-10.0*E0+float(ixp)*dE
     fixp = fixpp(ixp)+((xip-eixp)*(fixpp(ixp+1)-fixpp(ixp)))/dE
     r32=r32+r3*((omega_TA1(ip1))**2)* &
          ((fixp-fie)*phonon_TA1(ip1,it)+fixp*(1.0-fie))
  end do
!  write(*,*) 'r32=',r32
  r1=r1+r32

  r5=rph3/sqrt(eps(ie))
  do ip2=1,p
     xip=eps(ie)-omega_TA2(ip2)*dq
     ixp=int((xip-eps(1))/dE)
     eixp=1-10.0*E0+float(ixp)*dE
     fixp = fixpp(ixp)+((xip-eixp)*(fixpp(ixp+1)-fixpp(ixp)))/dE
     r1=r1+r5*((omega_TA2(ip2))**2)* &
          ((fixp-fie)*phonon_TA2(ip2,it)-fie*(1.0-fixp))

     xip=eps(ie)+omega_TA2(ip2)*dq
     ixp=int((xip-eps(1))/dE)
     eixp=1-10.0*E0+float(ixp)*dE
     fixp = fixpp(ixp)+((xip-eixp)*(fixpp(ixp+1)-fixpp(ixp)))/dE
     r1=r1+r5*((omega_TA2(ip2))**2)* &
          ((fixp-fie)*phonon_TA2(ip2,it)+fixp*(1.0-fie))
  end do
  recv(ie)=r1
  end do
  !end if
  
  dfermi = 0

  call MPI_BARRIER(MPI_COMM_WORLD, IERR)

  call MPI_ALLREDUCE(recv, dfermi, n, MPI_REAL, &
          MPI_SUM, MPI_COMM_WORLD, IERR)

 end subroutine deltaf
!----------------------------------------------------------------------------
!----------------------------------------------------------------------------
subroutine deltap_LA(it)

use mpi
implicit none

integer, intent(in) :: it
real :: r2,xi,fix, omega_dp, phonon, rph1, fie, eix
integer :: ip1, iomega_dp, phflag
real, allocatable :: fixp(:)

allocate(fixp(-n:2*n))

do ie=-n,0
  fixp(ie)=1.0
end do
do ie=1,n
  fixp(ie)=f(ie,it)
end do
do ie=n+1,2*n
  fixp(ie)=0.0
end do

phflag=0
rph1= domega_LA*par_cal_phe_LA
recv2 = 0.0

do iw = lpstart, lpend
!do iw = 1,m
  r2 = 0.0
  do ie=1,n
    fie=f(ie,it)
    xi=eps(ie)+omega_LA(iw)*dq
    ix=int((xi-eps(1))/dE)
    eix=1-10.0*E0+float(ix)*dE
    fix=fixp(ix)+((xi-eix)*(fixp(ix+1)-fixp(ix)))/dE
    r2=r2+(-1)*rph1*fie*phonon_LA(iw,it)*(1.0-fix)
    
    xi=eps(ie)-omega_LA(iw)*dq
    ix=int((xi-eps(1))/dE)
    eix=1-10.0*E0+float(ix)*dE
    fix=fixp(ix)+((xi-eix)*(fixp(ix+1)-fixp(ix)))/dE

    r2=r2+rph1*fie*(1.0+phonon_LA(iw,it))*(1.0-fix)
   end do
   
   if(phflag==0) then
   do ip1=1,min(p,int(iw*domega_LA/domega_TA1))
    omega_dp=omega_LA(iw)-omega_TA1(ip1)
    phonon=1.0
    iomega_dp=max(int((omega_dp-omega_LA(1))/domega_LA),1)
    phonon=phonon_LA(iomega_dp, it)+((omega_dp-omega_LA(iomega_dp))* &
              (phonon_LA(iomega_dp+1, it)-phonon_LA(iomega_dp, it)))/domega_LA
    r2=r2+domega_LA*(kph_ph)*(1.0/2.0)*cph_phm1(iw,ip1,1)*(((1.0+phonon_LA(iw,it))* &
            phonon_TA1(ip1,it)*phonon)-(phonon_LA(iw,it)*(1.0+phonon_TA1(ip1,it))*(1.0+phonon)))
    
    omega_dp=omega_LA(iw)+omega_TA1(ip1)
    iomega_dp=min(int((omega_dp-omega_LA(1))/domega_LA),m-1)
    phonon=phonon_LA(iomega_dp,it)+((omega_dp-omega_LA(iomega_dp))* &
                (phonon_LA(iomega_dp+1,it)-phonon_LA(iomega_dp,it)))/domega_LA
    !  else 
    !    phonon=phonon_LA(iomega_dp, it)+((omega_dp-omega_LA(iomega_dp))* &
    !          (-phonon_LA(iomega_dp,it)))/domega_LA
    !  end if
    r2=r2+domega_LA*(kph_ph)*cph_php1(iw,ip1,1)*(((1.0+phonon_LA(iw,it))* &
            (1.0+phonon_TA1(ip1,it))*phonon)-(phonon_LA(iw,it)* &
            phonon_TA1(ip1,it)*(1.0+phonon)))
  end do     
  end if
  recv2(iw) = r2
  !dphonon_LA(iw)=r2
end do

dphonon_LA = 0
call MPI_BARRIER(MPI_COMM_WORLD, IERR)

call MPI_ALLREDUCE(recv2, dphonon_LA, m, MPI_REAL, &
          MPI_SUM, MPI_COMM_WORLD, IERR)

end subroutine deltap_LA
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine deltap_TA1(it)

use mpi
implicit none

integer, intent(in) :: it
real :: r3,xi,fix, omega_dp, phonon, eix, fie, rph1
integer :: ip1, iomega_dp, ip1p, phflag
real, allocatable :: fixp(:)

allocate(fixp(-n:2*n))

do ie=-n,0
  fixp(ie)=1.0
end do
do ie=1,n
  fixp(ie)=f(ie,it)
end do
do ie=n+1,2*n
  fixp(ie)=0.0
end do

phflag=0
rph1=domega_TA1*par_cal_phe_TA1
recv3 = 0.0

do ip1 = tpstart, tpend
  r3 = 0.0
  do ie=1,n
    fie=f(ie,it)
    xi=eps(ie)+omega_TA1(ip1)*dq
    ix=int((xi-eps(1))/dE)
    eix=1-10.0*E0+float(ix)*dE
    fix=fixp(ix)+((xi-eix)*(fixp(ix+1)-fixp(ix)))/dE
    r3=r3+rph1*fie*(-1)*phonon_TA1(ip1,it)*(1.0-fix)
    
    xi=eps(ie)-omega_TA1(ip1)*dq
    ix=int((xi-eps(1))/dE)
    eix=1-10.0*E0+float(ix)*dE
    fix=fixp(ix)+((xi-eix)*(fixp(ix+1)-fixp(ix)))/dE
    r3=r3+rph1*fie*(1.0+phonon_TA1(ip1,it))*(1.0-fix)
  end do

  if(phflag==1) then
  do ip1p=1,p
    omega_dp=omega_TA1(ip1)-omega_TA1p(ip1p)
    if(omega_dp .lt. omega_TA1(1)) then
      phonon=1.0
    else
      iomega_dp=int((omega_dp-omega_TA1(1))/domega_TA1)
      phonon=phonon_TA1(iomega_dp,it)+((omega_dp-omega_TA1(iomega_dp))* &
              (phonon_TA1(iomega_dp+1,it)-phonon_TA1(iomega_dp,it)))/domega_TA1
    end if
    r3=r3+domega_TA1*(kph_ph)*(1.0/2.0)*cph_phm2(ip1,ip1p,2)*(((1.0+phonon_TA1(ip1,it))* &
            phonon_TA1p(ip1p,it)*phonon)-(phonon_TA1(ip1,it)*(1.0+phonon_TA1p(ip1p,it))*(1.0+phonon)))
    
    omega_dp=omega_TA1(ip1)+omega_TA1p(ip1p)
    if(omega_dp .gt. omega_TA1(p)) then
      phonon=0.0
    else
      iomega_dp=int((omega_dp-omega_TA1(1))/domega_TA1)
      if(iomega_dp .lt. p) then
        phonon=phonon_TA1(iomega_dp,it)+((omega_dp-omega_TA1(iomega_dp))* &
                (phonon_TA1(iomega_dp+1,it)-phonon_TA1(iomega_dp,it)))/domega_TA1
      else
        phonon=phonon_TA1(iomega_dp, it)+((omega_dp-omega_TA1(iomega_dp))* &
              (-phonon_TA1(iomega_dp,it)))/domega_TA1
      end if
    end if
    r3=r3+domega_TA1*(kph_ph)*cph_php2(ip1,ip1p,2)*(((1.0+phonon_TA1(ip1,it))* &
            (1.0+phonon_TA1p(ip1p,it))*phonon)-(phonon_TA1(ip1,it)*phonon_TA1p(ip1p,it)*(1.0+phonon)))
  end do
  end if ! phflag
  recv3(ip1) = r3
  !dphonon_TA1(ip1)=r3
end do
   dphonon_TA1 = 0.0
   call MPI_BARRIER(MPI_COMM_WORLD, IERR)

   call MPI_ALLREDUCE(recv3, dphonon_TA1, p, MPI_REAL, &
          MPI_SUM, MPI_COMM_WORLD, IERR)

end subroutine deltap_TA1
!-----------------------------------------------------------------------------
!-----------------------------------------------------------------------------
subroutine deltap_TA2(it)

implicit none

integer, intent(in) :: it
real :: r4,xi,fix,omega_dp,phonon,rph1,fie,eix
integer :: ip2, ip2p, iomega_dp,phflag
real, allocatable :: fixp(:)

allocate(fixp(-n:2*n))

do ie=-n,0
  fixp(ie)=1.0
end do
do ie=1,n
  fixp(ie)=f(ie,it)
end do
do ie=n+1,2*n
  fixp(ie)=0.0
end do

phflag=0
rph1=domega_TA2*par_cal_phe_TA2

do ip2 = 1,p
  r4 = 0.0
  do ie=1,n
    fie=f(ie,it)
    xi=eps(ie)+omega_TA2(ip2)*dq
    ix=int((xi-eps(1))/dE)
    eix=1-10.0*E0+float(ix)*dE
    fix=fixp(ix)+((xi-eix)*(fixp(ix+1)-fixp(ix)))/dE
    r4=r4+rph1*fie*(-1.0)*phonon_TA2(ip2,it)*(1.0-fix)

    xi=eps(ie)-omega_TA2(ip2)*dq
    ix=int((xi-eps(1))/dE)
    eix=1-10.0*E0+float(ix)*dE
    fix=fixp(ix)+((xi-eix)*(fixp(ix+1)-fixp(ix)))/dE
    r4=r4+rph1*fie*(1.0+phonon_TA2(ip2,it))*(1.0-fix)
  end do

  if(phflag==1) then
  do ip2p=1,p
    omega_dp=omega_TA2(ip2)-omega_TA2p(ip2p)
    if(omega_dp .lt. omega_TA2(1)) then
      phonon=1.0
    else
      iomega_dp=int((omega_dp-omega_TA2(1))/domega_TA2)
      phonon=phonon_TA2(iomega_dp,it)+((omega_dp-omega_TA2(iomega_dp))* &
              (phonon_TA2(iomega_dp+1,it)-phonon_TA2(iomega_dp,it)))/domega_TA2
    end if
    r4=r4+domega_TA2*(kph_ph)*(1.0/2.0)*cph_phm3(ip2,ip2p,3)*(((1.0+phonon_TA2(ip2,it))* &
            phonon_TA2p(ip2p,it)*phonon)-(phonon_TA2(ip2,it)*(1.0+phonon_TA2p(ip2p,it))*(1.0+phonon)))

    omega_dp=omega_TA2(ip2)+omega_TA2p(ip2p)
    if(omega_dp .gt. omega_TA2(p)) then
      phonon=0.0
    else
      iomega_dp=int((omega_dp-omega_TA2(1))/domega_TA2)
      if(iomega_dp .lt. p) then
        phonon=phonon_TA2(iomega_dp,it)+((omega_dp-omega_TA2(iomega_dp))* &
                (phonon_TA2(iomega_dp+1,it)-phonon_TA2(iomega_dp,it)))/domega_TA2
      else
        phonon=phonon_TA2(iomega_dp, it)+((omega_dp-omega_TA2(iomega_dp))* &
              (-phonon_TA2(iomega_dp,it)))/domega_TA2
      end if
    end if
    r4=r4+domega_TA2*(kph_ph)*cph_php3(ip2,ip2p,3)*(((1.0+phonon_TA2(ip2,it))* &
            (1.0+phonon_TA2p(ip2p,it))*phonon)-(phonon_TA2(ip2,it)*phonon_TA2p(ip2p,it)*(1.0+phonon)))
  end do
  end if !phflag
  dphonon_TA2(ip2)=r4
end do

end subroutine deltap_TA2

end module module_check
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
program f_eps

use mpi
use module_check

implicit none


real :: smax,smin, Aph_part, coeff, wla, wta1, wta2, costheta, sintheta, Cph_ph1, diracdf, &
        qdd, qpqp, qw, qwp, vtap, wta1p, wta2p, Cph_ph2, Cph_ph3, Eie, Eiep, Eix, Eixp, T1, T2, wdd
integer :: ipm, s, ik, il,dbf1,dbf2,it, unit1, unit2, unit3
real, dimension(Np) :: xx, yy, zz
real, allocatable, dimension(:) :: ka1, ka2, ka3, ka4, l1, l2, l3, l4, ta1, ta2, ta3, ta4, sa1, sa2, sa3, sa4
real, allocatable, dimension(:) :: ka, l, ta, sa

call MPI_INIT(IERR)
call MPI_COMM_SIZE(MPI_COMM_WORLD, SIZE, IERR)
call MPI_COMM_RANK(MPI_COMM_WORLD, RANK, IERR)

n=200
m=80
p=41
nit=5280000
root=0

npp = n/size
lpp = m/size
tpp = p/size

!dbf1: Only electron part will evolve; Phonons will be stationary
!dbf2: Electrons are stationary and decoupled from phonons; Phonons should not evolve and that can be tested

allocate(eps(n))
allocate(f(n,nit))
allocate(f_old(n))
allocate(pLA_old(n))
allocate(pTA1_old(n))
allocate(pTA2_old(n))
allocate(el_initial(n))
allocate(dfermi(n))
allocate(Cee(n,n,n))
allocate(omega_LA(m))
allocate(omega_TA1(p))
allocate(omega_TA2(p))
allocate(omega_TA1p(p))
allocate(omega_TA2p(p))
allocate(ph_initial_LA(m))
allocate(ph_initial_TA1(p))
allocate(ph_initial_TA2(p))
allocate(phonon_LA(m,nit))
allocate(phonon_TA1(p,nit))
allocate(phonon_TA1p(p,nit))
allocate(phonon_TA2(p,nit))
allocate(phonon_TA2p(p,nit))
allocate(dphonon_LA(m))
allocate(dphonon_TA1(p))
allocate(dphonon_TA2(p))
allocate(Cph_php1(m,m,3))
allocate(Cph_php2(m,m,3))
allocate(Cph_php3(m,m,3))
allocate(Cph_phm1(m,m,3))
allocate(Cph_phm2(m,m,3))
allocate(Cph_phm3(m,m,3))
allocate(recv(n))
allocate(recv2(m))
allocate(recv3(p))
allocate(recv4(n,n,n))


allocate(ka1(n))
allocate(ka2(n))
allocate(ka3(n))
allocate(ka4(n))
allocate(ka(n))
allocate(l1(n))
allocate(l2(n))
allocate(l3(n))
allocate(l4(n))
allocate(l(n))
allocate(ta1(n))
allocate(ta2(n))
allocate(ta3(n))
allocate(ta4(n))
allocate(ta(n))
allocate(sa1(n))
allocate(sa2(n))
allocate(sa3(n))
allocate(sa4(n))
allocate(sa(n))

de=(20.0*E0/float(n))

do ie = 1, n
  eps(ie) = 1.0 -10.0*E0 + float(ie)*de
end do

do iw = 1, m
  omega_LA(iw) = 0.0 + float(iw)*domega_LA 
end do

do ip1 = 1, p
  omega_TA1(ip1) = 0.0 + float(ip1)*domega_TA1
end do

do ip2 = 1, p
  omega_TA2(ip2) = 0.0 + float(ip2)*domega_TA2
end do

do ip1 = 1, p
  omega_TA1p(ip1) = 0.0 + float(ip1)*domega_TA1
end do

do ip2 = 1, p
  omega_TA2p(ip2) = 0.0 + float(ip2)*domega_TA2
end do

!open(unit = 1, file = 'in_el_dist.dat', status='unknown')
!write(1,*) '#Energy steps=',n,'Time steps=',nit,'full electron distribution when ph-ph is 0' 
!                    
!open(unit = 2, file = 'in_ph_dist.dat', status='unknown')
!write(2,*) '#Energy steps=',n,'Time steps=',nit,'full electron distribution when ph-ph is 0' 
!
!open(unit = 4, file = 'neq_el_dist.dat', status='unknown')
!write(4,*) '#Energy steps=',n,'Time steps=',nit,'full electron distribution when ph-ph is 0' 
!
!open(unit = 5, file = 'neq_ph_dist.dat', status='unknown')
!write(5,*) '#Energy steps=',n,'Time steps=',nit,'full electron distribution when ph-e ph-ph is 0' 


call in_cond_el(eps, el_initial)
call in_cond_ph_LA(omega_LA, ph_initial_LA)
call in_cond_ph_TA1(omega_TA1, ph_initial_TA1)
call in_cond_ph_TA2(omega_TA2, ph_initial_TA2)

zval=3.0
par_cal_eph_LA = (3.0*pi*zval*(Dla**2)*(hbaromd**3)*mM)/(16.0*(Kla**2)*(Ef**3))
par_cal_eph_TA1 = (3.0*pi*zval*(Dta1**2)*(hbaromd**3)*mM)/(16.0*(Kta1**2)*(Ef**3))
par_cal_eph_TA2 = par_cal_eph_TA1
par_cal_phe_LA = (3.0*sqrt(2.0)*pi*zval*(Dla**2)*mM)/(8.0*sqrt(Kla))
par_cal_phe_TA1 = (3.0*sqrt(2.0)*pi*zval*(Dta1**2)*mM)/(8.0*sqrt(Kta1))
par_cal_phe_TA2 = par_cal_phe_TA1

!write(*,*) 'eph',par_cal_eph_LA, par_cal_eph_TA1, par_cal_eph_TA2
!write(*,*) 'phe',par_cal_phe_LA, par_cal_phe_TA1, par_cal_phe_TA2

kph_ph= 0.0*(((2.0e9*pi*dq)/((4.0*pi)**4))*hbar*(omd**4))/(bpar*(v0**3))

f = 0.0
f(:,1)=el_initial

phonon_LA = 0.0
phonon_LA(:,1) = ph_initial_LA

phonon_TA1 = 0.0
phonon_TA1(:,1) = ph_initial_TA1

phonon_TA1p = 0.0
phonon_TA1p(:,1) = ph_initial_TA1

phonon_TA2 = 0.0
phonon_TA2(:,1) = ph_initial_TA2

phonon_TA2p = 0.0
phonon_TA2p(:,1) = ph_initial_TA2

istart = (rank*npp)+1
iend = (rank*npp)+npp
if(rank==size-1) iend=n

lpstart = (rank*lpp)+1
lpend = (rank*lpp)+lpp
if(rank==size-1) lpend=m

tpstart = (rank*tpp)+1
tpend = (rank*tpp)+tpp
if(rank==size-1) tpend=p

!write(6,*) "rank: ",rank,"istart: ", istart,"iend: ",iend
!write(*,*) ""
!write(*,*) ""
!write(*,*) ""
!write(*,*) ""
!write(6,*) "rank: ",rank,"lpstart: ", lpstart,"lpend: ",lpend
!write(*,*) ""
!write(*,*) ""
!write(*,*) ""
!write(*,*) ""
!write(6,*) "rank: ",rank,"tpstart: ", tpstart,"tpend: ",tpend
call MPI_BARRIER(MPI_COMM_WORLD, ierr)

!if(rank == 0) then
!  write(*,*) "Parameters reading is done and Cee is started"
!end if

call cpu_time(T1)
Cee(:,:,:)=0.0
do ie = istart, iend
!do ie = 1, n
  Eie=sqrt(eps(ie))
  do iep = 1, n
    Eiep=sqrt(eps(iep))
    do ix = 1, n
      Eix=sqrt(eps(ix))
      Eixp=sqrt(Eie**2+Eiep**2-Eix**2)
      smax = min((Eix+Eie),(Eixp+Eiep)) 
      smin = max((Eix-Eie),(Eixp-Eiep))
      !Cee(ie,iep,ix)=coupling_el(smax,smin)
      recv4(ie,iep,ix)=coupling_el(smax,smin)
    end do
  end do
end do
call cpu_time(T2)

!if (rank == root) then
!  write(*,*) T1, T2, "The time taken by Cee loop is: ", T2-T1 
!end if
!Cee = 0.0

call MPI_BARRIER(MPI_COMM_WORLD, IERR)

call MPI_ALLREDUCE(recv4, Cee, n*n*n, MPI_REAL, &
                    MPI_SUM, MPI_COMM_WORLD, IERR)

if(rank == root) then
  open(unit=3,file='inp.dat',status='unknown')
  do i=1,Np
    read(3,*) xx(i),yy(i),zz(i)
  end do
  close(3)
end if
call MPI_BCAST(xx, Np, MPI_REAL, root, MPI_COMM_WORLD, IERR)
call MPI_BCAST(yy, Np, MPI_REAL, root, MPI_COMM_WORLD, IERR)
call MPI_BCAST(zz, Np, MPI_REAL, root, MPI_COMM_WORLD, IERR)

!if(rank == root) then
!  write(*,*) "Input file reading is done"
!end if

cph_ph1=0.0
cph_php1=0.0
cph_phm1=0.0
cph_ph2=0.0
cph_php2=0.0
cph_phm2=0.0
cph_ph3=0.0
cph_php3=0.0
cph_phm3=0.0

!write(*,*) 'Cph-ph initialized'
dbf1=1
if(dbf1.eq.0) then
!write(*,*) 'Cph-ph started'
do s=1,3
if(s==1) then
  do i=1,m
    do j=1,p
      wla=omega_LA(i)*omd
      wta1=omega_TA1(j)*omd
      do ipm=-1,1,2
        wdd=wla+ipm*wta1
        qw=wla/vla
        qwp=wta1/vta
        qdd=wdd/vla
        vtap=vta/v0
        coeff=(1.0/(vtap**3))*(omega_TA1(j)/(omega_LA(i)*omega_LA(i)))
        do ik=1,Np
          x=xx(ik); y=yy(ik); z=zz(ik)
          do il=1,Np
            xp=xx(il); yp=yy(il); zp=zz(il)
            qpqp=(qw*x+qwp*xp)**2
            qpqp=qpqp+(qw*y+qwp*yp)**2
            qpqp=sqrt(qpqp+(qw*z+qwp*zp)**2)
            costheta=min(x*xp+y*yp+z*zp,1.0)
            sintheta=sqrt(1.0-(costheta)**2)
            !write(6,*) wdd
            Aph_part=(qw*qwp/qdd)*(2.0*E2+4.0*E3+lambda+3.0*mu)*(qw**2-qwp**2)*(costheta*sintheta)/ &
              (bpar*(q0**3))
            !write(6,*) Aph_part
            !stop
            diracdf=(1.0/(sigma*sqrt(2.0*pi)))*exp(-(wdd/omd - (vla/omd)*qpqp)**2/(2.0*sigma**2))
            cph_ph1=cph_ph1+coeff*(Aph_part**2)*diracdf
          end do !il
        end do !ik
        if(ipm==-1) then
            cph_phm1(i,j,s)=(4.0*pi/float(Np))*cph_ph1
        else
            cph_php1(i,j,s)=(4.0*pi/float(Np))*cph_ph1
        end if
      end do !ipm
    end do ! j=1,p
  end do !i=1,m
!  if(rank==root) then
!  write(6,*) sum(cph_php1)
!  write(6,*) sum(cph_phm1)
!  end if

!else if (s==2) then
!  do i=1,p
!   do j=1,i-1
!     do k=1,m
!       wta1=omega_TA1(i)*omd
!       wta1p=omega_TA1p(j)*omd
!       do ipm=-1,1,2
!         wdd=wta1+ipm*wta1p
!         qw=wta1/vta
!         qwp=wta1p/vta
!         qdd=wdd/vla
!         vtap=vta/v0
!         coeff=(1.0/(vtap**3))*(omega_TA1p(j)/(omega_TA1(i)*omega_LA(k)))
!         do ik=1,Np
!           x=xx(ik); y=yy(ik); z=zz(ik)
!           do il=1,Np
!             xp=xx(il); yp=yy(il); zp=zz(il)
!             qpqp=(qw*x+qwp*xp)**2
!             qpqp=qpqp+(qw*y+qwp*yp)**2
!             qpqp=sqrt(qpqp+(qw*z+qwp*zp)**2)
!             costheta=min(x*xp+y*yp+z*zp,1.0)
!             sintheta=sqrt(1.0-(costheta)**2)
!             Aph_part=(qw*qwp/qdd)*(((E2+lambda+2.0*(E3+mu))*(qwp-((qw*costheta))**2))- &
!                     ((E2+2.0*E3+mu)*(qw**2)*sintheta))
!
!             !write(6,*) Aph_part
!             diracdf=(1.0/(sigma*sqrt(2.0*pi)))*exp(-((wdd/omd - (vla/omd)*qpqp)**2)/(2.0*sigma**2))
!             cph_ph2=cph_ph2+coeff*(Aph_part**2)*diracdf
!           end do
!         end do
!         if(ipm==-1) then
!             cph_phm2(i,j,s)=(4.0*pi/float(Np))*cph_ph2
!         else
!             cph_php2(i,j,s)=(4.0*pi/float(Np))*cph_ph2
!         end if
!       end do
!     end do ! k=1,m
!    end do !j=1,p
! end do !i=1,p
!  if(rank==root) then
!  write(6,*) sum(cph_php2)
!  write(6,*) sum(cph_phm2)
!  end if

else
  do i=1,p
   do j=1,p
     do k=1,m
       wta2=omega_TA2(i)*omd
       wta2p=omega_TA2p(j)*omd
       do ipm=-1,1,2
         wdd=wta2+ipm*wta2p
         qw=wta2/vta
         qwp=wta2p/vta
         qdd=wdd/vla
         vtap=vta/v0
         coeff=(1.0/(vtap**3))*(omega_TA2p(j)/(omega_TA2(i)*omega_LA(k)))
         do ik=1,Np
           x=xx(ik); y=yy(ik); z=zz(ik)
           do il=1,Np
             xp=xx(il); yp=yy(il); zp=zz(il)
             qpqp=(qw*x+qwp*xp)**2
             qpqp=qpqp+(qw*y+qwp*yp)**2
             qpqp=sqrt(qpqp+(qw*z+qwp*zp)**2)
             costheta=min(x*xp+y*yp+z*zp,1.0)
             sintheta=sqrt(1.0-(costheta)**2)
             Aph_part=(qw*qwp)*(((E2+lambda)*((qw*costheta)-qwp))+(2.0*(E3+mu)* &
                     (qw-qwp*costheta)*costheta))
             !write(6,*) Aph_part
             diracdf=(1.0/(sigma*sqrt(2.0*pi)))*exp(-((wdd/omd - (vla/omd)*qpqp)**2)/(2.0*sigma**2))
             cph_ph3=cph_ph3+coeff*(Aph_part**2)*diracdf
           end do
         end do
         if(ipm==-1) then
             cph_phm3(i,j,s)=(4.0*pi/float(Np))*cph_ph3
         else
             cph_php3(i,j,s)=(4.0*pi/float(Np))*cph_ph3
         end if
       end do
     end do ! k=1,m
    end do !j=1,p
 end do !i=1,p
!  if(rank==root) then
!  write(6,*) sum(cph_php3)
!  write(6,*) sum(cph_phm3)
!  end if
end if   
end do
!else
        !READ from database - Set up
end if



!if (rank == root) then
!  write(*,*) "Time loop has started"
!end if  

unit1=1000
unit2=4000
unit3=7000

do it = 1, nit-1
  !write(*,*) "Time", it
  !if (rank == root) write(*,*) "Time", it
  f_old=f(:,it)
  pLA_old=phonon_LA(:,it)
  pTA1_old=phonon_TA1(:,it)
  pTA2_old=phonon_TA2(:,it)
  call deltaf(it)
  call deltap_LA(it)
  call deltap_TA1(it)
  call deltap_TA2(it)

  ka1 = dfermi
  l1 = dphonon_LA
  ta1= dphonon_TA1
  sa1 = dphonon_TA2

  f(:,it)=f(:,it)+dt*ka1/2.0
  phonon_LA(:,it)=phonon_LA(:,it)+dt*l1/2.0
  phonon_TA1(:,it)=phonon_TA1(:,it)+dt*ta1/2.0
  phonon_TA2(:,it)=phonon_TA2(:,it)+dt*sa1/2.0
  call deltaf(it)
  call deltap_LA(it)
  call deltap_TA1(it)
  call deltap_TA2(it)
  ka2 = dfermi
  l2 = dphonon_LA
  ta2= dphonon_TA1
  sa2 = dphonon_TA2


  f(:,it)=f(:,it)+dt*ka2/2.0
  phonon_LA(:,it)=phonon_LA(:,it)+dt*l2/2.0
  phonon_TA1(:,it)=phonon_TA1(:,it)+dt*ta2/2.0
  phonon_TA2(:,it)=phonon_TA2(:,it)+dt*sa2/2.0
  call deltaf(it)
  call deltap_LA(it)
  call deltap_TA1(it)
  call deltap_TA2(it)
  ka3 = dfermi
  l3 = dphonon_LA
  ta3= dphonon_TA1
  sa3 = dphonon_TA2

  f(:,it)=f(:,it)+dt*ka3
  phonon_LA(:,it)=phonon_LA(:,it)+dt*l3
  phonon_TA1(:,it)=phonon_TA1(:,it)+dt*ta3
  phonon_TA2(:,it)=phonon_TA2(:,it)+dt*sa3
  call deltaf(it)
  call deltap_LA(it)
  call deltap_TA1(it)
  call deltap_TA2(it)
  ka4 = dfermi
  l4 = dphonon_LA
  ta4= dphonon_TA1
  sa4 = dphonon_TA2

  !if(rank==root) then
  !        write(6,*) sum(ka1),sum(ka2),sum(ka3),sum(ka4)
  !end if
  f(:,it+1)=f_old+dt*(ka1+2.0*ka2+2.0*ka3+ka4)/6.0
  phonon_LA(:,it+1)=pLA_old+dt*(l1+2.0*l2+2.0*l3+l4)/6.0
  phonon_TA1(:,it+1)=pTA1_old+dt*(ta1+2.0*ta2+2.0*ta3+ta4)/6.0
  phonon_TA2(:,it+1)=pTA2_old+dt*(sa1+2.0*sa2+2.0*sa3+sa4)/6.0

 !if(rank == root) then
!  open(unit = 4, file = 'neq_el_dist.dat', status='unknown')
!  open(unit = 5, file = 'neq_ph_dist.dat', status='unknown')
  if(mod(it,1760)==0 .and. rank==0) then  
    do ie = 1, n
      write(unit1,1) (eps(ie)-1.0)/E0, f(ie, it)
    end do
    close(unit1)  
    unit1 = unit1+1  
    
    do iw = 1, m
      write(unit2, 1) omega_LA(iw), phonon_LA(iw, it)
    end do
    close(unit2)
    unit2 = unit2+1  

    do ip1 = 1, p
      write(unit3, 1) omega_TA1(ip1), phonon_TA1(ip1, it)
    end do
    close(unit3)
    unit3 = unit3+1  
  end if

1 format(f15.8, 4X, f15.8)
end do


deallocate(eps)
deallocate(f)
deallocate(f_old)
deallocate(pLA_old)
deallocate(pTA1_old)
deallocate(pTA2_old)
deallocate(el_initial)
deallocate(dfermi)
deallocate(Cee)
deallocate(omega_LA)
deallocate(omega_TA1)
deallocate(omega_TA2)
deallocate(omega_TA1p)
deallocate(omega_TA2p)
deallocate(ph_initial_LA)
deallocate(ph_initial_TA1)
deallocate(ph_initial_TA2)
deallocate(phonon_LA)
deallocate(phonon_TA1)
deallocate(phonon_TA1p)
deallocate(phonon_TA2)
deallocate(phonon_TA2p)
deallocate(dphonon_LA)
deallocate(dphonon_TA1)
deallocate(dphonon_TA2)
deallocate(Cph_php1)
deallocate(Cph_php2)
deallocate(Cph_php3)
deallocate(Cph_phm1)
deallocate(Cph_phm2)
deallocate(Cph_phm3)
deallocate(recv)
deallocate(recv2)
deallocate(recv3)
deallocate(recv4)


deallocate(ka1)
deallocate(ka2)
deallocate(ka3)
deallocate(ka4)
deallocate(ka)
deallocate(l1)
deallocate(l2)
deallocate(l3)
deallocate(l4)
deallocate(l)
deallocate(ta1)
deallocate(ta2)
deallocate(ta3)
deallocate(ta4)
deallocate(ta)
deallocate(sa1)
deallocate(sa2)
deallocate(sa3)
deallocate(sa4)
deallocate(sa)

call MPI_FINALIZE(IERR)

end program f_eps
