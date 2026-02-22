program ft2libre_loopback

! Loopback test: encode FT2Libre message, generate waveform, add noise,
! then call full decode chain (sync_ft2libre + ft2libreb) and check result.
! Runs multiple trials per SNR to get statistical decode rate.

  use packjt77
  use timer_module, only: timer

  include 'ft2libre_params.f90'
  parameter (NWAVE=NN2*NSPS)
  parameter (MAXCAND=200)

  real dd(NMAX)
  real candidate(3,MAXCAND)
  real sbase(NH1)
  complex cwave(NWAVE)
  real xjunk(NWAVE)
  integer itone(NN)
  integer itone_out(NN)
  integer*1 msgbits(77)
  character*37 msg37,msgsent37,decoded37
  character*12 mycall12,hiscall12
  character arg*12
  logical nagain,newdat,lsubtract,lapon,lapcqonly,ldecoded

  nargs=iargc()
  if(nargs.lt.4) then
     print*,'Usage: ft2libre_loopback "message" f0 ntrials snr1 [snr2 snr3...]'
     print*,'Example: ft2libre_loopback "K1ABC W9XYZ EN37" 1500.0 20 -14 -15 -16'
     go to 999
  endif
  call getarg(1,msg37)
  call getarg(2,arg)
  read(arg,*) f0
  call getarg(3,arg)
  read(arg,*) ntrials

  nsnr_points=nargs-3
  ndepth=3

  ! Encode (once)
  i3=-1
  n3=-1
  call genft2libre(msg37,i3,n3,msgsent37,msgbits,itone)
  write(*,'(a,a37)') 'Message:  ',msgsent37

  ! Generate complex waveform (once)
  call gen_ft2libre_wave(itone,NN,NSPS,1.0,12000.0,f0,cwave,xjunk,1,NWAVE)

  twopi=8.0*atan(1.0)
  fs=12000.0
  dt=1.0/fs
  bandwidth_ratio=2500.0/(fs/2.0)

  write(*,'(a,f8.1,a,i4)') 'f0=',f0,' ntrials=',ntrials
  write(*,'(a)') '  SNR    decoded/trials   rate'
  write(*,'(a)') '  ---    --------------   ----'

  call sgran()

  do isnr=1,nsnr_points
     call getarg(3+isnr,arg)
     read(arg,*) snrdb
     sig=sqrt(2*bandwidth_ratio) * 10.0**(0.05*snrdb)
     if(snrdb.gt.90.0) sig=1.0

     ngood=0
     do itrial=1,ntrials
        ! Place waveform in audio buffer at ~0.5s offset
        dd=0.0
        nstart=nint(0.5*fs)+1
        do i=1,NWAVE
           j=nstart+i-1
           if(j.ge.1.and.j.le.NMAX) dd(j)=sig*aimag(cwave(i))
        enddo

        ! Add noise
        if(snrdb.lt.90) then
           do i=1,NMAX
              dd(i)=dd(i)+gran()
           enddo
        endif

        ! Sync search
        nfa=200
        nfb=4000
        nfqso=nint(f0)
        syncmin=1.3
        maxc=MAXCAND
        call sync_ft2libre(dd,NMAX,nfa,nfb,syncmin,nfqso,maxc, &
             candidate,ncand,sbase)

        ! Try decoding each candidate
        mycall12='K1ABC       '
        hiscall12='W9XYZ       '
        nagain=.false.
        newdat=.true.
        lsubtract=.false.
        lapon=.false.
        lapcqonly=.false.
        ncontest=0
        nzhsym=NHSYM
        napwid=50

        ldecoded=.false.
        do ic=1,ncand
           if(ldecoded) exit
           f1=candidate(1,ic)
           xdt=candidate(2,ic)
           ibin=max(1,nint(f1/(12000.0/NFFT1)))
           ibin=min(ibin,NH1)
           xbase=10.0**(0.1*(sbase(ibin)-40.0))
           decoded37='                                     '
           iaptype=0

           do imetric=1,2
              call ft2libreb(dd,newdat,nfqso,ndepth,nzhsym,lsubtract, &
                   nagain,ncontest,imetric,f1,xdt,xbase,nharderrors,  &
                   dmin,nbadcrc,iaptype,decoded37,xsnr,itone_out,     &
                   mycall12,hiscall12,0,0,lapon,lapcqonly,napwid)
              if(nbadcrc.eq.0) then
                 ldecoded=.true.
                 exit
              endif
           enddo
        enddo
        if(ldecoded) ngood=ngood+1
     enddo
     rate=100.0*real(ngood)/real(ntrials)
     write(*,'(f6.1,i10,a,i4,f10.1,a)') snrdb,ngood,'/',ntrials,rate,'%'
  enddo

999 end program ft2libre_loopback
