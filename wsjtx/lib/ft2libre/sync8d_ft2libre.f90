subroutine sync8d_ft2libre(cd0,i0,ctwk,itwk,sync)

! Compute sync power for a complex, downsampled FT2Libre signal.
! Uses four 4-symbol Costas arrays with h=1.0 modulation index.

  include 'ft2libre_params.f90'
  parameter(NP2=NMAX/NDOWN)
  parameter(NSPSD=NSPS/NDOWN)           !32 samples per symbol downsampled
  complex cd0(0:NP2-1)
  complex csync(0:3,NSPSD)
  complex csync2(NSPSD)
  complex ctwk(NSPSD)
  complex z1
  logical first
  integer icos_a(0:3),icos_b(0:3),icos_c(0:3),icos_d(0:3)
  data icos_a/0,1,3,2/
  data icos_b/1,0,2,3/
  data icos_c/2,3,1,0/
  data icos_d/3,2,0,1/
  data first/.true./
  save first,twopi,csync

  p(z1)=real(z1)**2 + aimag(z1)**2          !Statement function for power

! Set some constants and compute the csync array.
  if( first ) then
    twopi=8.0*atan(1.0)
    do i=0,3
      phi=0.0
      dphi=twopi*i*1.0/real(NSPSD)          !h=1.0 modulation index
      do j=1,NSPSD
        csync(i,j)=cmplx(cos(phi),sin(phi))
        phi=mod(phi+dphi,twopi)
      enddo
    enddo
    first=.false.
  endif

  sync=0
! Costas A at position 0
  do i=0,3
     i1=i0+i*NSPSD
     csync2=csync(icos_a(i),1:NSPSD)
     if(itwk.eq.1) csync2=ctwk*csync2
     z1=0.
     if(i1.ge.0 .and. i1+NSPSD-1.le.NP2-1) z1=sum(cd0(i1:i1+NSPSD-1)*conjg(csync2))
     sync = sync + p(z1)
  enddo
! Costas B at position 33
  do i=0,3
     i1=i0+(33+i)*NSPSD
     csync2=csync(icos_b(i),1:NSPSD)
     if(itwk.eq.1) csync2=ctwk*csync2
     z1=0.
     if(i1.ge.0 .and. i1+NSPSD-1.le.NP2-1) z1=sum(cd0(i1:i1+NSPSD-1)*conjg(csync2))
     sync = sync + p(z1)
  enddo
! Costas C at position 66
  do i=0,3
     i1=i0+(66+i)*NSPSD
     csync2=csync(icos_c(i),1:NSPSD)
     if(itwk.eq.1) csync2=ctwk*csync2
     z1=0.
     if(i1.ge.0 .and. i1+NSPSD-1.le.NP2-1) z1=sum(cd0(i1:i1+NSPSD-1)*conjg(csync2))
     sync = sync + p(z1)
  enddo
! Costas D at position 99
  do i=0,3
     i1=i0+(99+i)*NSPSD
     csync2=csync(icos_d(i),1:NSPSD)
     if(itwk.eq.1) csync2=ctwk*csync2
     z1=0.
     if(i1.ge.0 .and. i1+NSPSD-1.le.NP2-1) z1=sum(cd0(i1:i1+NSPSD-1)*conjg(csync2))
     sync = sync + p(z1)
  enddo

  return
end subroutine sync8d_ft2libre
