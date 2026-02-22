subroutine ft2libreb(dd0,newdat,nfqso,ndepth,nzhsym,lsubtract,nagain, &
     ncontest,imetric,f1,xdt,xbase,nharderrors,dmin,nbadcrc,iaptype,  &
     msg37,xsnr,itone,mycall12,hiscall12,nQSOProgress,nftx,lapon,     &
     lapcqonly,napwid)

  use crc
  use timer_module, only: timer
  use packjt77
  include 'ft2libre_params.f90'
  parameter(NP2=NMAX/NDOWN)
  parameter(NSS=NSPS/NDOWN)            !32 samples per symbol downsampled
  parameter(NSYNC=2*NSS)               !64 samples for sync reference arrays
  parameter(NFFT_SYM=32)               !FFT size for soft symbols (= NSS)
  character*37 msg37
  character*77 c77
  character*12 mycall12,hiscall12
  real a(5)
  real s4(0:3,NN)
  real s2(0:255)
  real bitmetrics(2*NN,3)              !206 bit positions x 3 modes
  real bmet(2*NN)                      !Temporary for normalization
  real llra(174),llrb(174),llrc(174),llrd(174),llre(174),llrz(174)
  real dd0(NMAX)
  real temp(3)
  real rvec_pm1(77)                    !rvec in +/-1.0 form for AP scrambling
  real apmag
  integer*1 message77(77),message91(91),apmask(174),cw(174)
  integer*1 hbits(2*NN)
  integer itone(NN)
  integer ip(1)
  logical one(0:255,0:7)
  integer graymap(0:3)
  integer iloc(1)
  complex cd0(0:NP2-1)
  complex ctwk(NSYNC)
  complex csymb(NFFT_SYM)
  complex cs(0:3,NN)
  logical first,newdat,lsubtract,nagain,unpk77_success
  logical lapon,lapcqonly,lcross
  integer nQSOProgress,nftx,napwid,ntype

  integer icos_a(0:3),icos_b(0:3),icos_c(0:3),icos_d(0:3)
  data icos_a/0,1,3,2/
  data icos_b/1,0,2,3/
  data icos_c/2,3,1,0/
  data icos_d/3,2,0,1/

! Gray-coded sync bit patterns for hard quality check
  integer*1 sync_a(8),sync_b(8),sync_c(8),sync_d(8)
  data sync_a/0,0,0,1,1,0,1,1/
  data sync_b/0,1,0,0,1,1,1,0/
  data sync_c/1,1,1,0,0,1,0,0/
  data sync_d/1,0,1,1,0,0,0,1/

  integer*1 rvec(77)
  data rvec/0,1,0,0,1,0,1,0,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0,1,1,0,1,1,0, &
            1,0,0,1,0,1,1,0,0,0,0,1,0,0,0,1,0,1,0,0,1,1,1,1,0,0,1,0,1, &
            0,1,0,1,0,1,1,0,1,1,1,1,1,0,0,0,1,0,1/

! AP message patterns (same as FT8, in {0,1} form; converted to +/-1 in first block)
  integer mcq(29),mcqru(29),mcqfd(29),mcqtest(29),mcqww(29)
  integer mrrr(19),m73(19),mrr73(19)
  data     mcq/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0/
  data   mcqru/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,0,0,1,1,0,0/
  data   mcqfd/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0/
  data mcqtest/0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1,0,1,1,1,1,1,1,0,0,1,0/
  data   mcqww/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,1,1,1,1,0/
  data    mrrr/0,1,1,1,1,1,1,0,1,0,0,1,0,0,1,0,0,0,1/
  data     m73/0,1,1,1,1,1,1,0,1,0,0,1,0,1,0,0,0,0,1/
  data   mrr73/0,1,1,1,1,1,1,0,0,1,1,1,0,1,0,1,0,0,1/

! AP pass configuration tables (same as FT8)
  integer nappasses(0:5)
  integer naptypes(0:5,4)
  integer apsym0(58),aph10(10)
  integer ncontest0

  data first/.true./
  data graymap/0,1,3,2/
  save one,mcq,mcqru,mcqfd,mcqtest,mcqww,mrrr,m73,mrr73,       &
       rvec_pm1,nappasses,naptypes,ncontest0

  if(first) then
     one=.false.
     do i=0,255
       do j=0,7
         if(iand(i,2**j).ne.0) one(i,j)=.true.
       enddo
     enddo

! Convert CQ/response patterns to +/-1 form (once only)
     mcq=2*mcq-1
     mcqru=2*mcqru-1
     mcqfd=2*mcqfd-1
     mcqtest=2*mcqtest-1
     mcqww=2*mcqww-1
     mrrr=2*mrrr-1
     m73=2*m73-1
     mrr73=2*mrr73-1

! Compute rvec in +/-1 form for AP scrambling
     do i=1,77
       rvec_pm1(i)=1.0-2.0*rvec(i)
     enddo

! AP pass configuration (same as FT8)
     nappasses(0)=2
     nappasses(1)=2
     nappasses(2)=2
     nappasses(3)=4
     nappasses(4)=4
     nappasses(5)=3

! iaptype
!------------------------
!   1        CQ     ???    ???           (29+3=32 ap bits)
!   2        MyCall ???    ???           (29+3=32 ap bits)
!   3        MyCall DxCall ???           (58+3=61 ap bits)
!   4        MyCall DxCall RRR           (77 ap bits)
!   5        MyCall DxCall 73            (77 ap bits)
!   6        MyCall DxCall RR73          (77 ap bits)

     naptypes(0,1:4)=(/1,2,0,0/)  ! Tx6 selected (CQ)
     naptypes(1,1:4)=(/2,3,0,0/)  ! Tx1
     naptypes(2,1:4)=(/2,3,0,0/)  ! Tx2
     naptypes(3,1:4)=(/3,4,5,6/)  ! Tx3
     naptypes(4,1:4)=(/3,4,5,6/)  ! Tx4
     naptypes(5,1:4)=(/3,1,2,0/)  ! Tx5

     ncontest0=ncontest
     first=.false.
  endif

! Initialize AP symbols for current mycall/hiscall
  dxcall13=hiscall12
  mycall13=mycall12
  call ft8apset(mycall12,hiscall12,ncontest,apsym0,aph10)

  max_iterations=40                    !Fix 7: was 30
  nharderrors=-1
  iaptype=0
  fs2=12000.0/NDOWN
  dt2=1.0/fs2
  twopi=8.0*atan(1.0)
  delfbest=0.
  ibest=0

! Fix 5+12: Downsample and normalize
  call timer('f2l_down',0)
  call ft2libre_downsample(dd0,newdat,f1,cd0)
  call timer('f2l_down',1)
  sum2=sum(abs(cd0)**2)/(real(NMAX)/real(NDOWN))
  if(sum2.gt.0.0) cd0=cd0/sqrt(sum2)
  write(73,'(A,F12.6)') 'FT2L: sum2_initial=', sum2
  flush(73)

! Fix 6: 3-segment DT search with coarse/fine passes
! Matches decodium3 search ranges and strategy
  smax_all=0.0
  ibest_all=0
  delfbest_all=0.0

  do iseg=1,3
    if(iseg.eq.1) then
      ibmin=216; ibmax=1120
    elseif(iseg.eq.2) then
      ibmin=1120; ibmax=2024
    else
      ibmin=-688; ibmax=216
    endif

    ! Early exit for segments 2,3 if segment 1 sync too low
    if(iseg.ge.2 .and. smax_all.lt.0.10) exit

    smax_coarse=0.0
    ibest_coarse=ibmin
    delf_coarse=0.0

    ! Coarse pass: freq -12 to +12 Hz step 3, time step 4
    do ifr=-4,4
      delf=ifr*3.0
      dphi=twopi*delf*2.0*dt2
      phi=0.0
      do i=1,NSYNC
        ctwk(i)=cmplx(cos(phi),sin(phi))
        phi=mod(phi+dphi,twopi)
      enddo
      do idt=ibmin,ibmax,4
        call sync8d_ft2libre(cd0,idt,ctwk,1,sync)
        if(sync.gt.smax_coarse) then
          smax_coarse=sync
          ibest_coarse=idt
          delf_coarse=delf
        endif
      enddo
    enddo

    ! Fine pass: freq +/-4 Hz step 1, time +/-5 step 1
    smax_seg=smax_coarse
    ibest_seg=ibest_coarse
    delf_seg=delf_coarse

    do ifr=-4,4
      delf=delf_coarse+ifr*1.0
      dphi=twopi*delf*2.0*dt2
      phi=0.0
      do i=1,NSYNC
        ctwk(i)=cmplx(cos(phi),sin(phi))
        phi=mod(phi+dphi,twopi)
      enddo
      do idt=max(ibmin,ibest_coarse-5),min(ibmax,ibest_coarse+5)
        call sync8d_ft2libre(cd0,idt,ctwk,1,sync)
        if(sync.gt.smax_seg) then
          smax_seg=sync
          ibest_seg=idt
          delf_seg=delf
        endif
      enddo
    enddo

    write(73,'(A,I2,A,F10.3,A,I6,A,F8.2)') 'FT2L: seg=',iseg, &
         ' sync=',smax_seg,' ibest=',ibest_seg,' delf=',delf_seg
    flush(73)
    if(smax_seg.gt.smax_all) then
      smax_all=smax_seg
      ibest_all=ibest_seg
      delfbest_all=delf_seg
    endif
  enddo

  ibest=ibest_all
  delfbest=delfbest_all
  sync=smax_all

! Apply frequency correction and re-downsample
  f1=f1+delfbest
  call timer('f2l_down',0)
  call ft2libre_downsample(dd0,.false.,f1,cd0)
  call timer('f2l_down',1)
  sum2=sum(abs(cd0)**2)/(real(NMAX)/real(NDOWN))
  if(sum2.gt.0.0) cd0=cd0/sqrt(sum2)

! Fine DT refinement after frequency correction
  smax=0.0
  do idt=ibest-5,ibest+5
    call sync8d_ft2libre(cd0,idt,ctwk,0,sync)
    if(sync.gt.smax) then
      smax=sync
      ibest=idt
    endif
  enddo
  sync=smax
  write(73,'(A,F10.3,A,I6,A,F8.2,A,F8.3)') 'FT2L: final sync=',sync, &
       ' ibest=',ibest,' delf=',delfbest,' xdt=',ibest*dt2-0.5
  flush(73)

! Fix 9: DT reporting (matching decodium3: ibest/1333.33 - 0.5)
  xdt=ibest*dt2-0.5

! Note: removed signal-region-only normalization — the full-array
! normalization above is sufficient. Region-only normalization
! reduces effective SNR by normalizing out the signal power.

! Extract soft symbols using FFT (NFFT_SYM = NSS = 32)
  do k=1,NN
    i1=ibest+(k-1)*NSS
    csymb=cmplx(0.0,0.0)
    if( i1.ge.0 .and. i1+NSS-1 .le. NP2-1 ) csymb(1:NSS)=cd0(i1:i1+NSS-1)
    call four2a(csymb,NFFT_SYM,1,-1,1)
    cs(0:3,k)=csymb(1:4)/1e3
    s4(0:3,k)=abs(csymb(1:4))
  enddo

! Sync quality check against all four Costas arrays (hard-decoded tones)
  is1=0
  is2=0
  is3=0
  is4=0
  do k=0,3
    ip=maxloc(s4(:,k+1))
    if(icos_a(k).eq.(ip(1)-1)) is1=is1+1
    ip=maxloc(s4(:,k+34))
    if(icos_b(k).eq.(ip(1)-1)) is2=is2+1
    ip=maxloc(s4(:,k+67))
    if(icos_c(k).eq.(ip(1)-1)) is3=is3+1
    ip=maxloc(s4(:,k+100))
    if(icos_d(k).eq.(ip(1)-1)) is4=is4+1
  enddo
  nhsync=is1+is2+is3+is4
  write(73,'(A,4I3,A,I4)') 'FT2L: is1/is2/is3/is4=',is1,is2,is3,is4,' nhsync=',nhsync
  flush(73)
  if(nhsync.lt.4) then
    nbadcrc=1
    return
  endif

! Compute bit metrics over ALL NN=103 symbols (including sync)
! Produces bitmetrics(2*NN, 3) = 206 bit positions x 3 modes
! For multi-symbol modes (nseq=2,3), skip groups that cross
! sync/data boundaries to avoid contaminating bit metrics.
  bitmetrics=0.0
  do nseq=1,3
    if(nseq.eq.1) nsym_actual=1
    if(nseq.eq.2) nsym_actual=2
    if(nseq.eq.3) nsym_actual=4
    nt=4**nsym_actual

    do ks=1,NN-nsym_actual+1,nsym_actual
! Skip multi-symbol groups crossing sync/data boundaries:
! Sync at 1-4, 34-37, 67-70, 100-103; Data at 5-33, 38-66, 71-99
      if(nsym_actual.gt.1) then
        ks_end=ks+nsym_actual-1
        lcross=.false.
        if(ks.le.4 .and. ks_end.ge.5) lcross=.true.
        if(ks.le.33 .and. ks_end.ge.34) lcross=.true.
        if(ks.le.37 .and. ks_end.ge.38) lcross=.true.
        if(ks.le.66 .and. ks_end.ge.67) lcross=.true.
        if(ks.le.70 .and. ks_end.ge.71) lcross=.true.
        if(ks.le.99 .and. ks_end.ge.100) lcross=.true.
        if(lcross) cycle
      endif
      do i=0,nt-1
        if(nsym_actual.eq.1) then
          i1=i
          s2(i)=abs(cs(graymap(i1),ks))
        elseif(nsym_actual.eq.2) then
          i1=i/4
          i2=iand(i,3)
          s2(i)=abs(cs(graymap(i1),ks)+cs(graymap(i2),ks+1))
        elseif(nsym_actual.eq.4) then
          i1=i/64
          i2=iand(i,63)/16
          i3=iand(i,15)/4
          i4=iand(i,3)
          s2(i)=abs(cs(graymap(i1),ks)+cs(graymap(i2),ks+1)+  &
               cs(graymap(i3),ks+2)+cs(graymap(i4),ks+3))
        endif
      enddo
      if(imetric.eq.2) s2=s2**2

      ipt=1+(ks-1)*2
      if(nsym_actual.eq.1) ibmax=1
      if(nsym_actual.eq.2) ibmax=3
      if(nsym_actual.eq.4) ibmax=7
      do ib=0,ibmax
        if(ipt+ib.gt.2*NN) cycle
        bm=maxval(s2(0:nt-1),one(0:nt-1,ibmax-ib)) - &
           maxval(s2(0:nt-1),.not.one(0:nt-1,ibmax-ib))
        bitmetrics(ipt+ib,nseq)=bm
      enddo
    enddo
  enddo

! Fall back to single-symbol metrics where multi-symbol groups were
! skipped (boundary crossings) or incomplete (end of symbol sequence)
  do ipt=1,2*NN
    if(bitmetrics(ipt,2).eq.0.0) bitmetrics(ipt,2)=bitmetrics(ipt,1)
    if(bitmetrics(ipt,3).eq.0.0) bitmetrics(ipt,3)=bitmetrics(ipt,2)
  enddo

! Sync quality check on raw mode-1 bitmetrics (before normalization)
  bmet=bitmetrics(:,1)
  hbits=0
  where(bmet.ge.0.0) hbits=1
  ns1=count(hbits(1:8).eq.sync_a)
  ns2=count(hbits(67:74).eq.sync_b)
  ns3=count(hbits(133:140).eq.sync_c)
  ns4=count(hbits(199:206).eq.sync_d)
  nsync_qual=ns1+ns2+ns3+ns4
  write(73,'(A,4I3,A,I4)') 'FT2L: ns1/ns2/ns3/ns4=',ns1,ns2,ns3,ns4,' nsync_qual=',nsync_qual
  flush(73)
  if(nsync_qual.lt.15) then
    nbadcrc=1
    return
  endif

! Extract data bits by skipping sync positions, then normalize
! Sync: bit positions 1-8, 67-74, 133-140, 199-206
! Data: bit positions 9-66, 75-132, 141-198
! Normalize over 174 data positions only (matching FT8 approach)
  scalefac=2.83

  llra(1:58)=bitmetrics(9:66,1)
  llra(59:116)=bitmetrics(75:132,1)
  llra(117:174)=bitmetrics(141:198,1)
  call normalizebmet(llra,174)
  write(73,'(A,3F10.4)') 'FT2L: llra min/max/rms=', &
       minval(llra),maxval(llra),sqrt(sum(llra**2)/174.0)
  flush(73)
  llra=scalefac*llra

  llrb(1:58)=bitmetrics(9:66,2)
  llrb(59:116)=bitmetrics(75:132,2)
  llrb(117:174)=bitmetrics(141:198,2)
  call normalizebmet(llrb,174)
  llrb=scalefac*llrb

  llrc(1:58)=bitmetrics(9:66,3)
  llrc(59:116)=bitmetrics(75:132,3)
  llrc(117:174)=bitmetrics(141:198,3)
  call normalizebmet(llrc,174)
  llrc=scalefac*llrc

! llrd = best-of (max absolute value among a,b,c)
! llre = average of a,b,c
  do i=1,174
    temp(1)=llra(i)
    temp(2)=llrb(i)
    temp(3)=llrc(i)
    ip=maxloc(abs(temp))
    llrd(i)=temp(ip(1))
    llre(i)=(temp(1)+temp(2)+temp(3))/3.0
  enddo

! Determine number of decode passes (5 base + AP passes)
  npasses=5
  if(lapon.or.ncontest.eq.7) then
     if(.not.lapcqonly) then
        npasses=5+2*nappasses(nQSOProgress)
     else
        npasses=7
     endif
  endif
  if(nzhsym.lt.50) npasses=5

! Decode passes (Fix 7: maxosd=3, norder=3)
  do ipass=1,npasses
     llrz=llra
     if(ipass.eq.2) llrz=llrb
     if(ipass.eq.3) llrz=llrc
     if(ipass.eq.4) llrz=llrd
     if(ipass.eq.5) llrz=llre
     if(ipass.le.5) then
        apmask=0
        iaptype=0
     endif

! AP passes (6+): apply a priori constraints with rvec scrambling
     if(ipass.gt.5) then
        if(mod(ipass-5,2).eq.1) llrz=llra
        if(mod(ipass-5,2).eq.0) llrz=llrc
        apmag=maxval(abs(llrz))*1.1
        if(.not.lapcqonly) then
           iaptype=naptypes(nQSOProgress,(ipass-4)/2)
        else
           iaptype=1
        endif

! Conditions that cause us to bail out of AP decoding
        if(ncontest.le.5 .and. iaptype.ge.3 .and.                     &
           (abs(f1-nfqso).gt.napwid .and.                              &
            abs(f1-nftx).gt.napwid)) cycle
        if(ncontest.eq.6) cycle                     ! No AP for Foxes
        if(ncontest.eq.7.and.f1.gt.950.0) cycle     ! Hounds: AP below 950 Hz only
        if(iaptype.ge.2 .and. apsym0(1).gt.1) cycle ! No/nonstandard mycall
        if(ncontest.eq.7 .and. iaptype.ge.2 .and. aph10(1).gt.1) cycle
        if(iaptype.ge.3 .and. apsym0(30).gt.1) cycle ! No/nonstandard dxcall

! FT2 AP: all known message-bit constraints must be rvec-scrambled.
! rvec_pm1(i) = 1-2*rvec(i) converts rvec from {0,1} to {+1,-1}.
! For a known message bit v (in +/-1) at position p, the LDPC info bit
! constraint is v * rvec_pm1(p).

        if(iaptype.eq.1) then ! CQ or CQ RU or CQ TEST or CQ FD
           apmask=0
           apmask(1:29)=1
           if(ncontest.eq.0) llrz(1:29)=apmag*mcq(1:29)*rvec_pm1(1:29)
           if(ncontest.eq.1) llrz(1:29)=apmag*mcqtest(1:29)*rvec_pm1(1:29)
           if(ncontest.eq.2) llrz(1:29)=apmag*mcqtest(1:29)*rvec_pm1(1:29)
           if(ncontest.eq.3) llrz(1:29)=apmag*mcqfd(1:29)*rvec_pm1(1:29)
           if(ncontest.eq.4) llrz(1:29)=apmag*mcqru(1:29)*rvec_pm1(1:29)
           if(ncontest.eq.5) llrz(1:29)=apmag*mcqww(1:29)*rvec_pm1(1:29)
           if(ncontest.eq.7) llrz(1:29)=apmag*mcq(1:29)*rvec_pm1(1:29)
           if(ncontest.eq.8) llrz(1:29)=apmag*mcqtest(1:29)*rvec_pm1(1:29)
           apmask(75:77)=1
           llrz(75)=apmag*(-1)*rvec_pm1(75)
           llrz(76)=apmag*(-1)*rvec_pm1(76)
           llrz(77)=apmag*(+1)*rvec_pm1(77)
        endif

        if(iaptype.eq.2) then ! MyCall,???,???
           apmask=0
           if(ncontest.eq.0.or.ncontest.eq.1.or.ncontest.eq.5              &
              .or.ncontest.eq.8) then
              apmask(1:29)=1
              llrz(1:29)=apmag*apsym0(1:29)*rvec_pm1(1:29)
              apmask(75:77)=1
              llrz(75)=apmag*(-1)*rvec_pm1(75)
              llrz(76)=apmag*(-1)*rvec_pm1(76)
              llrz(77)=apmag*(+1)*rvec_pm1(77)
           else if(ncontest.eq.2) then
              apmask(1:28)=1
              llrz(1:28)=apmag*apsym0(1:28)*rvec_pm1(1:28)
              apmask(72:74)=1
              llrz(72)=apmag*(-1)*rvec_pm1(72)
              llrz(73)=apmag*(+1)*rvec_pm1(73)
              llrz(74)=apmag*(-1)*rvec_pm1(74)
              apmask(75:77)=1
              llrz(75)=apmag*(-1)*rvec_pm1(75)
              llrz(76)=apmag*(-1)*rvec_pm1(76)
              llrz(77)=apmag*(-1)*rvec_pm1(77)
           else if(ncontest.eq.3) then
              apmask(1:28)=1
              llrz(1:28)=apmag*apsym0(1:28)*rvec_pm1(1:28)
              apmask(75:77)=1
              llrz(75)=apmag*(-1)*rvec_pm1(75)
              llrz(76)=apmag*(-1)*rvec_pm1(76)
              llrz(77)=apmag*(-1)*rvec_pm1(77)
           else if(ncontest.eq.4) then
              apmask(2:29)=1
              llrz(2:29)=apmag*apsym0(1:28)*rvec_pm1(2:29)
              apmask(75:77)=1
              llrz(75)=apmag*(-1)*rvec_pm1(75)
              llrz(76)=apmag*(+1)*rvec_pm1(76)
              llrz(77)=apmag*(+1)*rvec_pm1(77)
           else if(ncontest.eq.7) then
              apmask(29:56)=1
              llrz(29:56)=apmag*apsym0(1:28)*rvec_pm1(29:56)
              apmask(57:66)=1
              llrz(57:66)=apmag*aph10(1:10)*rvec_pm1(57:66)
              apmask(72:77)=1
              llrz(72)=apmag*(-1)*rvec_pm1(72)
              llrz(73)=apmag*(-1)*rvec_pm1(73)
              llrz(74)=apmag*(+1)*rvec_pm1(74)
              llrz(75)=apmag*(-1)*rvec_pm1(75)
              llrz(76)=apmag*(-1)*rvec_pm1(76)
              llrz(77)=apmag*(-1)*rvec_pm1(77)
           endif
        endif

        if(iaptype.eq.3) then ! MyCall,DxCall,???
           apmask=0
           if(ncontest.eq.0.or.ncontest.eq.1.or.ncontest.eq.2              &
              .or.ncontest.eq.5.or.ncontest.eq.7.or.ncontest.eq.8) then
              apmask(1:58)=1
              llrz(1:58)=apmag*apsym0(1:58)*rvec_pm1(1:58)
              apmask(75:77)=1
              llrz(75)=apmag*(-1)*rvec_pm1(75)
              llrz(76)=apmag*(-1)*rvec_pm1(76)
              llrz(77)=apmag*(+1)*rvec_pm1(77)
           else if(ncontest.eq.3) then
              apmask(1:56)=1
              llrz(1:28)=apmag*apsym0(1:28)*rvec_pm1(1:28)
              llrz(29:56)=apmag*apsym0(30:57)*rvec_pm1(29:56)
              apmask(72:74)=1
              apmask(75:77)=1
              llrz(75)=apmag*(-1)*rvec_pm1(75)
              llrz(76)=apmag*(-1)*rvec_pm1(76)
              llrz(77)=apmag*(-1)*rvec_pm1(77)
           else if(ncontest.eq.4) then
              apmask(2:57)=1
              llrz(2:29)=apmag*apsym0(1:28)*rvec_pm1(2:29)
              llrz(30:57)=apmag*apsym0(30:57)*rvec_pm1(30:57)
              apmask(75:77)=1
              llrz(75)=apmag*(-1)*rvec_pm1(75)
              llrz(76)=apmag*(+1)*rvec_pm1(76)
              llrz(77)=apmag*(+1)*rvec_pm1(77)
           endif
        endif

        if(iaptype.eq.5.and.ncontest.eq.7) cycle !Hound
        if(iaptype.eq.4 .or. iaptype.eq.5 .or. iaptype.eq.6) then
           apmask=0
           if(ncontest.le.5 .or. (ncontest.eq.7.and.iaptype.eq.6)          &
              .or. ncontest.eq.8) then
              apmask(1:77)=1
              llrz(1:58)=apmag*apsym0(1:58)*rvec_pm1(1:58)
              if(iaptype.eq.4) llrz(59:77)=apmag*mrrr(1:19)*rvec_pm1(59:77)
              if(iaptype.eq.5) llrz(59:77)=apmag*m73(1:19)*rvec_pm1(59:77)
              if(iaptype.eq.6) llrz(59:77)=apmag*mrr73(1:19)*rvec_pm1(59:77)
           else if(ncontest.eq.7.and.iaptype.eq.4) then
              apmask(1:28)=1
              llrz(1:28)=apmag*apsym0(1:28)*rvec_pm1(1:28)
              apmask(57:66)=1
              llrz(57:66)=apmag*aph10(1:10)*rvec_pm1(57:66)
              apmask(72:77)=1
              llrz(72)=apmag*(-1)*rvec_pm1(72)
              llrz(73)=apmag*(-1)*rvec_pm1(73)
              llrz(74)=apmag*(+1)*rvec_pm1(74)
              llrz(75)=apmag*(-1)*rvec_pm1(75)
              llrz(76)=apmag*(-1)*rvec_pm1(76)
              llrz(77)=apmag*(-1)*rvec_pm1(77)
           endif
        endif
     endif

     cw=0
     dmin=0.0
     norder=3                          !Fix 7: was 2
     maxosd=3                          !Fix 7: was 2
     if(ndepth.eq.1) maxosd=-1
     call timer('dec174_91 ',0)
     Keff=91
     call decode174_91(llrz,Keff,maxosd,norder,apmask,message91,cw,  &
                       ntype,nharderrors,dmin)
     write(73,'(A,I3,A,I4,A,F8.2,A,I3)') 'FT2L: pass=',ipass, &
          ' nharderr=',nharderrors,' dmin=',dmin,' iaptype=',iaptype
     flush(73)
     if(nharderrors.ge.0) then
! Descramble message bits with rvec
        do i=1,77
           message77(i)=ieor(message91(i),rvec(i))
        enddo
     endif
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
     call get_ft2libre_tones_from_77bits(message77,itone)
     if(lsubtract) then
        call timer('sub_f2l ',0)
        call subtractft2libre(dd0,itone,f1,xdt)
        call timer('sub_f2l ',1)
     endif
     xsig=0.0
     xnoi=0.0
     do i=1,NN
        xsig=xsig+s4(itone(i),i)**2
        ios=mod(itone(i)+2,4)
        xnoi=xnoi+s4(ios,i)**2
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
     if(nhsync.le.8 .and. xsnr.lt.-25.0) then
       nbadcrc=1
       return
     endif
     if(xsnr .lt. -25.0) xsnr=-25.0
     return
  enddo
  return
end subroutine ft2libreb
