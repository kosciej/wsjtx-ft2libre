subroutine subtractft2libre(dd0,itone,f0,dt)

! Subtract an FT2Libre signal
!
! Measured signal  : dd(t)    = a(t)cos(2*pi*f0*t+theta(t))
! Reference signal : cref(t)  = exp( j*(2*pi*f0*t+phi(t)) )
! Complex amp      : cfilt(t) = LPF[ dd(t)*CONJG(cref(t)) ]
! Subtract         : dd(t)    = dd(t) - 2*REAL{cref*cfilt}

  parameter (NMAX=45000,NFRAME=NMAX)
  parameter (NFFT=NMAX,NFILT=600)
  real dd(NMAX),dd0(NMAX)
  real window(-NFILT/2:NFILT/2)
  real x(NFFT+2)
  real endcorrection(NFILT/2+1)
  complex cx(0:NFFT/2)
  complex cref,camp,cfilt,cw,z
  integer itone(79)
  logical first
  data first/.true./
  common/heap2libre/cref(NFRAME),camp(NMAX),cfilt(NMAX),cw(NMAX)
  equivalence (x,cx)
  save first,/heap2libre/,endcorrection

  if(first) then                         ! Create and normalize the filter
     pi=4.0*atan(1.0)
     fac=1.0/float(nfft)
     sumw=0.0
     do j=-NFILT/2,NFILT/2
        window(j)=cos(pi*j/NFILT)**2
        sumw=sumw+window(j)
     enddo
     cw=0.
     cw(1:NFILT+1)=window/sumw
     cw=cshift(cw,NFILT/2+1)
     call four2a(cw,nfft,1,-1,1)
     cw=cw*fac
     first=.false.
     do j=1,NFILT/2+1
       endcorrection(j)=1.0/(1.0-sum(window(j-1:NFILT/2))/sumw)
     enddo
  endif

! Generate complex reference waveform cref (h=0.75 via gen_ft2libre_wave)
  call gen_ft2libre_wave(itone,79,360,2.0,12000.0,f0,cref,xjunk,1,NFRAME)

  nstart=dt*12000+1
  camp=0.
  dd=dd0
  do i=1,nframe
     j=nstart-1+i
     if(j.ge.1.and.j.le.NMAX) camp(i)=dd(j)*conjg(cref(i))
  enddo

  cfilt(1:nframe)=camp(1:nframe)
  cfilt(nframe+1:)=0.0
  call four2a(cfilt,nfft,1,-1,1)
  cfilt(1:nfft)=cfilt(1:nfft)*cw(1:nfft)
  call four2a(cfilt,nfft,1,1,1)
  cfilt(1:NFILT/2+1)=cfilt(1:NFILT/2+1)*endcorrection
  cfilt(nframe:nframe-NFILT/2:-1)=cfilt(nframe:nframe-NFILT/2:-1)*endcorrection

  do i=1,nframe
     j=nstart+i-1
     if(j.ge.1 .and. j.le.NMAX) then
        z=cfilt(i)*cref(i)
        dd(j)=dd(j)-2.0*real(z)      !Subtract the reconstructed signal
     endif
  enddo
  dd0=dd                               !Return dd0 with this signal subtracted

  return
end subroutine subtractft2libre
