subroutine ft2libreb(dd0,newdat,nfqso,ndepth,nzhsym,lsubtract,nagain, &
     ncontest,imetric,f1,xdt,xbase,nharderrors,dmin,nbadcrc,ipass,    &
     msg37,xsnr,itone)

  use crc
  use timer_module, only: timer
  use packjt77
  include 'ft2libre_params.f90'
  parameter(NP2=NMAX/NDOWN)
  parameter(NSPSD=NSPS/NDOWN)           !30 samples per symbol downsampled
  parameter(NFFT_SYM=40)                !FFT size for soft symbols (zero-padded)
  character*37 msg37
  character*77 c77
  real a(5)
  real s8(0:7,NN)
  real s2(0:511)
  real bmeta(174),bmetb(174),bmetc(174),bmetd(174),bmete(174)
  real llra(174),llrb(174),llrc(174),llrd(174),llre(174),llrz(174)
  real dd0(NMAX)
  real ss(9)
  real temp(3)
  integer*1 message77(77),message91(91),apmask(174),cw(174)
  integer itone(NN)
  integer icos7(0:6),ip(1)
  logical one(0:511,0:8)
  integer graymap(0:7)
  integer iloc(1)
  complex cd0(0:NP2-1)
  complex ctwk(NSPSD)
  complex csymb(NFFT_SYM)
  complex cs(0:7,NN)
  logical first,newdat,lsubtract,nagain,unpk77_success
  data icos7/3,1,4,0,6,5,2/
  data first/.true./
  data graymap/0,1,3,2,5,6,4,7/
  save one

  if(first) then
     one=.false.
     do i=0,511
       do j=0,8
         if(iand(i,2**j).ne.0) one(i,j)=.true.
       enddo
     enddo
     first=.false.
  endif

  max_iterations=30
  nharderrors=-1
  fs2=12000.0/NDOWN
  dt2=1.0/fs2
  twopi=8.0*atan(1.0)
  delfbest=0.
  ibest=0

  call timer('f2l_down',0)
  call ft2libre_downsample(dd0,newdat,f1,cd0)
  call timer('f2l_down',1)

  i0=nint((xdt+0.5)*fs2)
  smax=0.0
  do idt=i0-10,i0+10
     call sync8d_ft2libre(cd0,idt,ctwk,0,sync)
     if(sync.gt.smax) then
        smax=sync
        ibest=idt
     endif
  enddo

! Now peak up in frequency
  smax=0.0
  do ifr=-5,5
    delf=ifr*0.5
    dphi=twopi*delf*dt2
    phi=0.0
    do i=1,NSPSD
      ctwk(i)=cmplx(cos(phi),sin(phi))
      phi=mod(phi+dphi,twopi)
    enddo
    call sync8d_ft2libre(cd0,ibest,ctwk,1,sync)
    if( sync .gt. smax ) then
      smax=sync
      delfbest=delf
    endif
  enddo
  a=0.0
  a(1)=-delfbest
  call twkfreq1(cd0,NP2,fs2,a,cd0)
  f1=f1+delfbest

  call timer('f2l_down',0)
  call ft2libre_downsample(dd0,.false.,f1,cd0)
  call timer('f2l_down',1)

  smax=0.0
  do idt=-4,4
     call sync8d_ft2libre(cd0,ibest+idt,ctwk,0,sync)
     ss(idt+5)=sync
  enddo
  smax=maxval(ss)
  iloc=maxloc(ss)
  ibest=iloc(1)-5+ibest
  xdt=(ibest-1)*dt2
  sync=smax

! Extract soft symbols using zero-padded FFT
  do k=1,NN
    i1=ibest+(k-1)*NSPSD
    csymb=cmplx(0.0,0.0)
    if( i1.ge.0 .and. i1+NSPSD-1 .le. NP2-1 ) csymb(1:NSPSD)=cd0(i1:i1+NSPSD-1)
    call four2a(csymb,NFFT_SYM,1,-1,1)
    cs(0:7,k)=csymb(1:8)/1e3
    s8(0:7,k)=abs(csymb(1:8))
  enddo

! sync quality check
  is1=0
  is2=0
  is3=0
  do k=1,7
    ip=maxloc(s8(:,k))
    if(icos7(k-1).eq.(ip(1)-1)) is1=is1+1
    ip=maxloc(s8(:,k+36))
    if(icos7(k-1).eq.(ip(1)-1)) is2=is2+1
    ip=maxloc(s8(:,k+72))
    if(icos7(k-1).eq.(ip(1)-1)) is3=is3+1
  enddo
  nsync=is1+is2+is3
  syncmin=6
  if(imetric.eq.2) syncmin=7
  if(ndepth.le.2) syncmin=8
  if(nsync.le.syncmin) then
    nbadcrc=1
    return
  endif

  do nsym=1,3
    nt=2**(3*nsym)
    do ihalf=1,2
      do k=1,29,nsym
        if(ihalf.eq.1) ks=k+7
        if(ihalf.eq.2) ks=k+43
        amax=-1.0
        do i=0,nt-1
          i1=i/64
          i2=iand(i,63)/8
          i3=iand(i,7)
          if(nsym.eq.1) then
            s2(i)=abs(cs(graymap(i3),ks))
          elseif(nsym.eq.2) then
            s2(i)=abs(cs(graymap(i2),ks)+cs(graymap(i3),ks+1))
          elseif(nsym.eq.3) then
            s2(i)=abs(cs(graymap(i1),ks)+cs(graymap(i2),ks+1)+  &
                 cs(graymap(i3),ks+2))
          else
            print*,"Error - nsym must be 1, 2, or 3."
          endif
        enddo
        if(imetric.eq.2) s2=s2**2
        i32=1+(k-1)*3+(ihalf-1)*87
        if(nsym.eq.1) ibmax=2
        if(nsym.eq.2) ibmax=5
        if(nsym.eq.3) ibmax=8
        do ib=0,ibmax
          bm=maxval(s2(0:nt-1),one(0:nt-1,ibmax-ib)) - &
             maxval(s2(0:nt-1),.not.one(0:nt-1,ibmax-ib))
          if(i32+ib .gt.174) cycle
          if(nsym.eq.1) then
            bmeta(i32+ib)=bm
            den=max(maxval(s2(0:nt-1),one(0:nt-1,ibmax-ib)), &
                    maxval(s2(0:nt-1),.not.one(0:nt-1,ibmax-ib)))
            if(den.gt.0.0) then
              cm=bm/den
            else
              cm=0.0
            endif
            bmetd(i32+ib)=cm
          elseif(nsym.eq.2) then
            bmetb(i32+ib)=bm
          elseif(nsym.eq.3) then
            bmetc(i32+ib)=bm
          endif
        enddo
      enddo
    enddo
  enddo
  do i=1,174
    temp(1)=bmeta(i)
    temp(2)=bmetb(i)
    temp(3)=bmetc(i)
    ip=maxloc(abs(temp))
    bmete(i)=temp(ip(1))
  enddo

  call normalizebmet(bmeta,174)
  call normalizebmet(bmetb,174)
  call normalizebmet(bmetc,174)
  call normalizebmet(bmetd,174)
  call normalizebmet(bmete,174)

  scalefac=2.83
  llra=scalefac*bmeta
  llrb=scalefac*bmetb
  llrc=scalefac*bmetc
  llrd=scalefac*bmetd
  llre=scalefac*bmete

! No AP decoding for FT2Libre (simplified)
  npasses=5

  iaptype=0
  do ipass=1,npasses
     llrz=llra
     if(ipass.eq.2) llrz=llrb
     if(ipass.eq.3) llrz=llrc
     if(ipass.eq.4) llrz=llrd
     if(ipass.eq.5) llrz=llre
     apmask=0

     cw=0
     dmin=0.0
     norder=2
     maxosd=2
     if(ndepth.eq.1) maxosd=-1
     call timer('dec174_91 ',0)
     Keff=91
     call decode174_91(llrz,Keff,maxosd,norder,apmask,message91,cw,  &
                       nharderrors,dmin)
     if(nharderrors.ge.0) message77=message91(1:77)
     call timer('dec174_91 ',1)

     msg37='                                     '
     nbadcrc=1
     if(nharderrors.lt.0 .or. nharderrors.gt.36) cycle
     if(count(cw.eq.0).eq.174) cycle
     write(c77,'(77i1)') message77
     read(c77(72:74),'(b3)') n3
     read(c77(75:77),'(b3)') i3
     if(i3.gt.5 .or. (i3.eq.0.and.n3.gt.6)) cycle
     if(i3.eq.0 .and. n3.eq.2) cycle
     call unpack77(c77,1,msg37,unpk77_success)
     if(.not.unpk77_success) cycle
     nbadcrc=0
     call get_ft8_tones_from_77bits(message77,itone)
     if(lsubtract) then
        call timer('sub_f2l ',0)
        call subtractft2libre(dd0,itone,f1,xdt)
        call timer('sub_f2l ',1)
     endif
     xsig=0.0
     xnoi=0.0
     do i=1,79
        xsig=xsig+s8(itone(i),i)**2
        ios=mod(itone(i)+4,7)
        xnoi=xnoi+s8(ios,i)**2
     enddo
     xsnr=0.001
     xsnr2=0.001
     arg=xsig/xnoi-1.0
     if(arg.gt.0.1) xsnr=arg
     arg=xsig/xbase/3.0e6-1.0
     if(arg.gt.0.1) xsnr2=arg
     xsnr=10.0*log10(xsnr)-27.0
     xsnr2=10.0*log10(xsnr2)-27.0
     if(.not.nagain) then
       xsnr=xsnr2
     endif
     if(nsync.le.10 .and. xsnr.lt.-25.0) then
       nbadcrc=1
       return
     endif
     if(xsnr .lt. -25.0) xsnr=-25.0
     return
  enddo
  return
end subroutine ft2libreb
