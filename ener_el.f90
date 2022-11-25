Program energy_el_ph

integer, parameter :: ne=200
integer :: i
real :: E0, T0, Ef
real, dimension(ne)::ener, f, fd

E0=0.005558
T0=0.025/11.695
Ef=1.0

open(unit=1,file='dist_in.dat',status='unknown')
do i=1,ne
  read(1,*) ener(i),f(i)
  ener(i)=ener(i)*E0 + 1.0
end do
de=ener(2)-ener(1)

do i=1,ne
  if(ener(i).le.1.0) then
    fd(i)=1.0/(exp((ener(i)-Ef)/T0)+ 1.0)
  else
    fd(i)=exp(-(ener(i)-Ef)/T0)/(exp(-(ener(i)-Ef)/T0)+ 1.0)
  end if
end do

toten=(2056.46*de*sum(sqrt(ener)*ener*(f-fd)))+6
!toten=sum((f-fd))
write(6,*) 'Tot energy=',toten
open(unit=11,file='toten.dat',status='unknown')
write(11,*) toten
close(11)


open(unit=11,file='dist_out.dat',status='unknown')
do i=1,ne
write(11,*) ener(i),f(i),fd(i)
end do
close(11)


end program energy_el_ph





