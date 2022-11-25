Program energy_el_ph

integer, parameter :: np=80
integer :: i
real :: T0, wla
real, dimension(np)::ener, ph, phbe

T0=0.025/0.065

open(unit=1,file='dist_in.dat',status='unknown')
do i=1,np
  read(1,*) ener(i),ph(i)
end do
de=ener(2)-ener(1)

do i=1,np
    phbe(i)=1.0/(exp(ener(i)/T0)- 1.0)
end do

toten=(12657.83184/400)*de*sum((ener)**3*(ph-phbe))
write(6,*) 'Tot energy=',toten
open(unit=11,file='toten.dat',status='unknown')
write(11,*) toten
close(11)


open(unit=11,file='dist_out.dat',status='unknown')
do i=1,np
write(11,*) ener(i),ph(i),phbe(i)
end do
close(11)


end program energy_el_ph





