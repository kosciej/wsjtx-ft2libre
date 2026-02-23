subroutine sync8d_ft2libre(cd0,i0,ctwk,itwk,sync)

! Compute sync power for a complex, downsampled FT2Libre signal.
! Uses four 4-symbol Costas arrays with h=0.8 modulation index.
! Matches decodium3 sync2d approach: continuous-phase coherent correlation
! across all 4 symbols of each Costas array, stride-2 access, abs (not squared).

  include 'ft2libre_params.f90'
  parameter(NP2=NMAX/NDOWN)
  parameter(NSS=NSPS/NDOWN)            !32 samples per symbol downsampled
  parameter(NSYNC=2*NSS)               !64 samples per reference array
  complex cd0(0:NP2-1)
  complex csynca(NSYNC),csyncb(NSYNC),csyncc(NSYNC),csyncd(NSYNC)
  complex csync2(NSYNC)
  complex ctwk(NSYNC)
  complex z1,z2,z3,z4
  real fac
  logical first
  integer icos_a(0:3),icos_b(0:3),icos_c(0:3),icos_d(0:3)
  data icos_a/0,1,3,2/
  data icos_b/1,0,2,3/
  data icos_c/2,3,1,0/
  data icos_d/3,2,0,1/
  data first/.true./
  save first,twopi,csynca,csyncb,csyncc,csyncd

! Build continuous-phase reference arrays (computed once, saved).
! Each array spans 4 symbols but uses NSS/2 = 16 samples per symbol,
! yielding 2*NSS = 64 total reference samples.
! Phase accumulates continuously from symbol to symbol.
  if( first ) then
    twopi=8.0*atan(1.0)

    ! Costas A
    k=1; phi=0.0
    do i=0,3
      dphi=2.0*twopi*0.8*icos_a(i)/real(NSS)
      do j=1,NSS/2
        csynca(k)=cmplx(cos(phi),sin(phi))
        phi=mod(phi+dphi,twopi)
        k=k+1
      enddo
    enddo

    ! Costas B
    k=1; phi=0.0
    do i=0,3
      dphi=2.0*twopi*0.8*icos_b(i)/real(NSS)
      do j=1,NSS/2
        csyncb(k)=cmplx(cos(phi),sin(phi))
        phi=mod(phi+dphi,twopi)
        k=k+1
      enddo
    enddo

    ! Costas C
    k=1; phi=0.0
    do i=0,3
      dphi=2.0*twopi*0.8*icos_c(i)/real(NSS)
      do j=1,NSS/2
        csyncc(k)=cmplx(cos(phi),sin(phi))
        phi=mod(phi+dphi,twopi)
        k=k+1
      enddo
    enddo

    ! Costas D
    k=1; phi=0.0
    do i=0,3
      dphi=2.0*twopi*0.8*icos_d(i)/real(NSS)
      do j=1,NSS/2
        csyncd(k)=cmplx(cos(phi),sin(phi))
        phi=mod(phi+dphi,twopi)
        k=k+1
      enddo
    enddo

    first=.false.
  endif

  fac=1.0/(2.0*NSS)

! Correlate with stride-2 access into downsampled signal.
! Each correlation spans 4*NSS samples but takes every 2nd sample.
  i1=i0                      ! Costas A at position 0
  i2=i0+33*NSS               ! Costas B at position 33
  i3=i0+66*NSS               ! Costas C at position 66
  i4=i0+99*NSS               ! Costas D at position 99

  z1=cmplx(0.0,0.0)
  z2=cmplx(0.0,0.0)
  z3=cmplx(0.0,0.0)
  z4=cmplx(0.0,0.0)

  if(itwk.eq.1) then
    csync2=ctwk*csynca
  else
    csync2=csynca
  endif
  if(i1.ge.0 .and. i1+4*NSS-1.le.NP2-1) then
    z1=sum(cd0(i1:i1+4*NSS-1:2)*conjg(csync2))
  endif

  if(itwk.eq.1) then
    csync2=ctwk*csyncb
  else
    csync2=csyncb
  endif
  if(i2.ge.0 .and. i2+4*NSS-1.le.NP2-1) then
    z2=sum(cd0(i2:i2+4*NSS-1:2)*conjg(csync2))
  endif

  if(itwk.eq.1) then
    csync2=ctwk*csyncc
  else
    csync2=csyncc
  endif
  if(i3.ge.0 .and. i3+4*NSS-1.le.NP2-1) then
    z3=sum(cd0(i3:i3+4*NSS-1:2)*conjg(csync2))
  endif

  if(itwk.eq.1) then
    csync2=ctwk*csyncd
  else
    csync2=csyncd
  endif
  if(i4.ge.0 .and. i4+4*NSS-1.le.NP2-1) then
    z4=sum(cd0(i4:i4+4*NSS-1:2)*conjg(csync2))
  endif

  sync=abs(z1*fac)+abs(z2*fac)+abs(z3*fac)+abs(z4*fac)

  return
end subroutine sync8d_ft2libre
