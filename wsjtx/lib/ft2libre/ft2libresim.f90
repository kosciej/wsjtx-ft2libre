program ft2libresim

! Generate simulated FT2Libre .wav files for loopback testing.
! Usage: ft2libresim "message" f0 DT snr nfiles

  use wavhdr
  use packjt77
  include 'ft2libre_params.f90'
  parameter (NWAVE=NN2*NSPS)              !30240 samples including ramps
  type(hdr) h
  character arg*12,fname*17
  character msg37*37,msgsent37*37
  character c77*77
  complex cwave(NWAVE)
  real wave(NMAX)
  real xjunk(NWAVE)
  integer itone(NN)
  integer*1 msgbits(77)
  integer*2 iwave(NMAX)

  nargs=iargc()
  if(nargs.ne.5) then
     print*,'Usage:    ft2libresim "message" f0 DT nfiles snr'
     print*,'Example:  ft2libresim "K1ABC W9XYZ EN37" 1500.0 0.0 10 -15'
     go to 999
  endif
  call getarg(1,msg37)
  call getarg(2,arg)
  read(arg,*) f0
  call getarg(3,arg)
  read(arg,*) xdt
  call getarg(4,arg)
  read(arg,*) nfiles
  call getarg(5,arg)
  read(arg,*) snrdb

  twopi=8.0*atan(1.0)
  fs=12000.0
  dt=1.0/fs
  tt=NSPS*dt
  baud=1.0/tt
  bandwidth_ratio=2500.0/(fs/2.0)
  sig=sqrt(2*bandwidth_ratio) * 10.0**(0.05*snrdb)
  if(snrdb.gt.90.0) sig=1.0

  ! Encode message to channel symbols
  i3=-1
  n3=-1
  call pack77(msg37,i3,n3,c77)
  call genft2libre(msg37,i3,n3,msgsent37,msgbits,itone)

  write(*,'(a,a37)') 'Message:  ',msgsent37
  write(*,'(a,i1,a,i1)') 'i3=',i3,' n3=',n3
  write(*,'(a,f9.3,a,f6.2,a,f6.1)') 'f0:',f0,'  DT:',xdt,'  SNR:',snrdb
  write(*,'(a)') 'Channel symbols:'
  write(*,'(20i2)') itone
  write(*,*)

  ! Generate complex waveform (includes ramp symbols)
  call gen_ft2libre_wave(itone,NN,NSPS,1.0,fs,f0,cwave,xjunk,1,NWAVE)

  call sgran()

  do ifile=1,nfiles
     wave=0.0
     ! Place waveform at correct time offset
     ! Nominal start is at 0.5s into the file
     nstart=nint((xdt+0.5)*fs)+1
     do i=1,NWAVE
        j=nstart+i-1
        if(j.ge.1.and.j.le.NMAX) then
           wave(j)=sig*aimag(cwave(i))
        endif
     enddo

     if(snrdb.lt.90) then
        do i=1,NMAX
           wave(i)=wave(i)+gran()
        enddo
     endif

     gain=100.0
     if(snrdb.lt.90.0) then
        wave=gain*wave
     else
        datpk=maxval(abs(wave))
        fac=32766.9/datpk
        wave=fac*wave
     endif
     if(any(abs(wave).gt.32767.0)) print*,"Warning - data will be clipped."
     iwave=nint(wave)
     h=default_header(12000,NMAX)
     write(fname,1102) ifile
1102 format('000000_',i6.6,'.wav')
     open(10,file=fname,status='unknown',access='stream')
     write(10) h,iwave
     close(10)
     write(*,1110) ifile,xdt,f0,snrdb,fname
1110 format(i4,f7.2,f8.2,f7.1,2x,a17)
  enddo

999 end program ft2libresim
