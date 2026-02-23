subroutine gen_ft2libre_wave(itone,nsym,nsps,bt,fsample,f0,cwave,wave,icmplx,nwave)
!
! Generate FT2Libre waveform using Gaussian-filtered frequency pulses.
! 4-GFSK modulation with h=0.8.
! Input:  itone(nsym) = channel symbols (103 sync+data)
! Output: wave/cwave of length nwave = (nsym+2)*nsps including ramp symbols
!
  use timer_module, only: timer
  parameter(MAX_SECONDS=20,NTAB=65536)
  real wave(nwave)
  complex cwave(nwave),ctab(0:NTAB-1)
  real pulse(23040)
  real dphi(0:(nsym+2)*nsps-1)
  integer itone(nsym)
  data fchk0/0.0/
  save pulse,twopi,dt,hmod,fchk0,ctab

  ibt=nint(10*1.0)
  fchk=nsym+nsps+1.0+fsample
  if(fchk.ne.fchk0) then
     twopi=8.0*atan(1.0)
     dt=1.0/fsample
     hmod=0.8
! Compute the frequency-smoothing pulse
     do i=1,3*nsps
        tt=(i-1.5*nsps)/real(nsps)
        pulse(i)=gfsk_pulse(1.0,tt)
     enddo
     do i=0,NTAB-1
        phi=i*twopi/NTAB
        ctab(i)=cmplx(cos(phi),sin(phi))
     enddo
     fchk0=fchk
  endif

! Compute the smoothed frequency waveform.
! Length = (nsym+2)*nsps samples, first and last symbols extended
  dphi_peak=twopi*hmod/real(nsps)
  dphi=0.0
  do j=1,nsym
     ib=(j-1)*nsps
     ie=ib+3*nsps-1
     dphi(ib:ie) = dphi(ib:ie) + dphi_peak*pulse(1:3*nsps)*itone(j)
  enddo
! No dummy symbols in ramp regions — ramps carry carrier only (matching decodium3)

! Calculate and insert the audio waveform
! Output includes ramp-up and ramp-down symbols: nwave = (nsym+2)*nsps
  phi=0.0
  dphi = dphi + twopi*f0*dt                      !Shift frequency up by f0
  if(icmplx .eq. 0) wave=0.
  if(icmplx .ne. 0) cwave=0.

  k=0
  do j=0,nwave-1                                 !Include ramp symbols
     k=k+1
     if(icmplx.eq.0) then
        wave(k)=sin(phi)
     else
        i=phi*float(NTAB)/twopi
        cwave(k)=ctab(i)
     endif
     phi=mod(phi+dphi(j),twopi)
  enddo

! Apply cosine taper to ramp-up (first nsps samples) and ramp-down (last nsps samples)
  nramp=nsps
  if(icmplx.eq.0) then
     wave(1:nramp)=wave(1:nramp) *                                          &
          (1.0-cos(twopi*(/(i,i=0,nramp-1)/)/(2.0*nramp)))/2.0
     k1=nwave-nramp+1
     wave(k1:nwave)=wave(k1:nwave) *                                        &
          (1.0+cos(twopi*(/(i,i=0,nramp-1)/)/(2.0*nramp)))/2.0
  else
     cwave(1:nramp)=cwave(1:nramp) *                                        &
          (1.0-cos(twopi*(/(i,i=0,nramp-1)/)/(2.0*nramp)))/2.0
     k1=nwave-nramp+1
     cwave(k1:nwave)=cwave(k1:nwave) *                                      &
          (1.0+cos(twopi*(/(i,i=0,nramp-1)/)/(2.0*nramp)))/2.0
  endif

  return
end subroutine gen_ft2libre_wave
