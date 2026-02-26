module ft2l_decode

   type :: ft2l_decoder
      procedure(ft2l_decode_callback), pointer :: callback
   contains
      procedure :: decode
   end type ft2l_decoder

   abstract interface
      subroutine ft2l_decode_callback (this,sync,snr,dt,freq,decoded,nap,qual)
         import ft2l_decoder
         implicit none
         class(ft2l_decoder), intent(inout) :: this
         real, intent(in) :: sync
         integer, intent(in) :: snr
         real, intent(in) :: dt
         real, intent(in) :: freq
         character(len=37), intent(in) :: decoded
         integer, intent(in) :: nap
         real, intent(in) :: qual
      end subroutine ft2l_decode_callback
   end interface

contains

   subroutine decode(this,callback,iwave,nQSOProgress,nfqso,    &
      nfa,nfb,ndepth,lapcqonly,ncontest,mycall,hiscall)
      use timer_module, only: timer
      use packjt77
      include 'ft2l/ft2l_params.f90'
      parameter (MAXCAND=200)
      class(ft2l_decoder), intent(inout) :: this
      procedure(ft2l_decode_callback) :: callback
      parameter (NSS=NSPS/NDOWN,NDMAX=NMAX/NDOWN)
      character message*37,msgsent*37
      character c77*77
      character*37 decodes(100)
      character*17 cdatetime0
      character*12 mycall,hiscall
      character*12 mycall0,hiscall0
      character*6 hhmmss

      complex cd2(0:NDMAX-1)                  !Complex waveform
      complex cb(0:NDMAX-1)
      complex cd(0:NN*NSS-1)                       !Complex waveform
      complex ctwk(2*NSS),ctwk2(2*NSS,-42:42)

      real a(5)
      real bitmetrics(2*NN,5)
      real dd(NMAX)
      real llr(2*ND),llra(2*ND),llrb(2*ND),llrc(2*ND),llrd(2*ND),llre(2*ND)
      real candidate(2,MAXCAND)
      real savg(NH1),sbase(NH1)

      integer apbits(2*ND)
      integer apmy_ru(28),aphis_fd(28)
      integer*2 iwave(NMAX)                 !Raw received data
      integer*1 message77(77),rvec(77),apmask(2*ND),cw(2*ND)
      integer*1 message91(91)
      integer*1 hbits(2*NN)
      integer i4tone(103)
      integer nappasses(0:5)    ! # of decoding passes for QSO States 0-5
      integer naptypes(0:5,4)   ! nQSOProgress, decoding pass
      integer mcq(29),mcqru(29),mcqfd(29),mcqtest(29),mcqww(29)
      integer mrrr(19),m73(19),mrr73(19)

      logical nohiscall,unpk77_success
      logical first, dobigfft
      logical dosubtract,doosd
      logical badsync
      logical, intent(in) :: lapcqonly

      data first/.true./
      data     mcq/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0/
      data   mcqru/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,1,1,1,0,0,1,1,0,0/
      data   mcqfd/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,0,1,0,0,0,1,0/
      data mcqtest/0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,1,0,1,0,1,1,1,1,1,1,0,0,1,0/
      data   mcqww/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,1,0,1,1,1,1,0/      
      data    mrrr/0,1,1,1,1,1,1,0,1,0,0,1,0,0,1,0,0,0,1/
      data     m73/0,1,1,1,1,1,1,0,1,0,0,1,0,1,0,0,0,0,1/
      data   mrr73/0,1,1,1,1,1,1,0,0,1,1,1,0,1,0,1,0,0,1/
      data rvec/0,1,0,0,1,0,1,0,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0,1,1,0,1,1,0, &
         1,0,0,1,0,1,1,0,0,0,0,1,0,0,0,1,0,1,0,0,1,1,1,1,0,0,1,0,1, &
         0,1,0,1,0,1,1,0,1,1,1,1,1,0,0,0,1,0,1/
      save fs,dt,tt,txt,twopi,h,first,apbits,nappasses,naptypes, &
         mycall0,hiscall0,ctwk2

      this%callback => callback
      hhmmss=cdatetime0(8:13)
      dxcall13=hiscall        ! initialize for use in packjt77
      mycall13=mycall

      smax1=0.
      nd1=0

      if(first) then
         fs=12000.0/NDOWN                !Sample rate after downsampling
         dt=1/fs                         !Sample interval after downsample (s)
         tt=NSPS*dt                      !Duration of "itone" symbols (s)
         txt=NZ*dt                       !Transmission length (s) without ramp up/down
         twopi=8.0*atan(1.0)
         h=1.0
         
         do idf=-42,42
            a=0.
            a(1)=real(idf)
            ctwk=1.
            call twkfreq1(ctwk,2*NSS,fs/2.0,a,ctwk2(:,idf))
         enddo

         nappasses(0)=2
         nappasses(1)=2
         nappasses(2)=2
         nappasses(3)=2
         nappasses(4)=2
         nappasses(5)=3

         naptypes(0,1:2)=(/1,2/)
         naptypes(1,1:2)=(/1,2/)
         naptypes(2,1:2)=(/1,2/)
         naptypes(3,1:2)=(/1,3/)
         naptypes(4,1:2)=(/1,3/)
         naptypes(5,1:3)=(/1,3,6/)

         mycall0=''
         hiscall0=''
         first=.false.
      endif

      l1=index(mycall,char(0))
      if(l1.ne.0) mycall(l1:)=" "
      l1=index(hiscall,char(0))
      if(l1.ne.0) hiscall(l1:)=" "

      if(first .or. mycall.ne.mycall0 .or. hiscall.ne.hiscall0) then
         apbits=0
         apbits(1)=99
         apbits(30)=99
         apmy_ru=0
         aphis_fd=0

         if(len(trim(mycall)) .lt. 3) go to 10

         nohiscall=.false.
         hiscall0=hiscall
         if(len(trim(hiscall0)).lt.3) then
            hiscall0=mycall  ! use mycall for dummy hiscall - mycall won't be hashed.
            nohiscall=.true.
         endif
         message=trim(mycall)//' '//trim(hiscall0)//' RR73'
         i3=-1
         n3=-1
         call pack77(message,i3,n3,c77)
         call unpack77(c77,1,msgsent,unpk77_success)
         if(i3.ne.1 .or. (message.ne.msgsent) .or. .not.unpk77_success) go to 10
         read(c77,'(77i1)') message77
         apmy_ru=2*mod(message77(1:28)+rvec(2:29),2)-1
         aphis_fd=2*mod(message77(30:57)+rvec(29:56),2)-1
         message77=mod(message77+rvec,2)
         call encode174_91(message77,cw)
         apbits=2*cw-1
         if(nohiscall) apbits(30)=99

10       continue
         mycall0=mycall
         hiscall0=hiscall
      endif

      dd=iwave/1000.0
      ndecodes=0
      decodes=''
      fa=nfa
      fb=nfb


! ndepth=3: 5 passes, including osd
! ndepth=2: 3 passes, bp only
! ndepth=1: 1 pass, no subtraction

      max_iterations=100
      syncmin=1.05
      dosubtract=.true.
      doosd=.true.
      nsp=3
      if(ndepth.eq.2) then
         doosd=.false.
      endif
      if(ndepth.eq.1) then
         nsp=1
         dosubtract=.false.
         doosd=.false.
      endif

      do isp = 1,nsp
         if(isp.eq.2) then
            if(ndecodes.eq.0) exit
            nd1=ndecodes
         elseif(isp.eq.3) then
            nd2=ndecodes-nd1
            if(nd2.eq.0) exit
         endif

         candidate=0.0
         ncand=0
         call timer('getcand2l',0)
         call getcandidates2l(dd,fa,fb,syncmin,nfqso,MAXCAND,savg,candidate,   &
            ncand,sbase)
         call timer('getcand2l',1)
         dobigfft=.true.
         do icand=1,ncand
            f0=candidate(1,icand)
            snr=candidate(2,icand)-1.0
            call timer('ft2l_down',0)
            call ft2l_downsample(dd,dobigfft,f0,cd2)  !Downsample to 32 Sam/Sym
            call timer('ft2l_down',1)
            if(dobigfft) dobigfft=.false.
            sum2=sum(cd2*conjg(cd2))/(real(NMAX)/real(NDOWN))
            if(sum2.gt.0.0) cd2=cd2/sqrt(sum2)
! Sample rate is now 12000/9 = 1333.33 samples/second
            do iseg=1,3                ! DT search is done over 3 segments
               do isync=1,2          
                  if(isync.eq.1) then
                     idfmin=-24
                     idfmax=24
                     idfstp=6
                     ibmin=-688
                     ibmax=2024
                     if(iseg.eq.1) then
                        ibmin=216
                        ibmax=1120
                     elseif(iseg.eq.2) then
                        smax1=smax
                        ibmin=1120
                        ibmax=2024
                     elseif(iseg.eq.3) then
                        ibmin=-688
                        ibmax=216
                     endif
                     ibstp=1
                  else
                     idfmin=idfbest-8
                     idfmax=idfbest+8
                     idfstp=1
                     ibmin=ibest-10
                     ibmax=ibest+10
                     ibstp=1
                  endif
                  ibest=-1
                  idfbest=0
                  smax=-99.
                  call timer('sync2ld ',0)
                  do idf=idfmin,idfmax,idfstp
                     do istart=ibmin,ibmax,ibstp
                        call sync2ld(cd2,istart,ctwk2(:,idf),1,sync)  !Find sync power
                        if(sync.gt.smax) then
                           maxit=1
                           smax=sync
                           ibest=istart
                           idfbest=idf
                        endif
                     enddo
                  enddo
                  call timer('sync2ld ',1)
               enddo
               if(iseg.eq.1) smax1=smax
               if(smax.lt.1.0) cycle
               if(iseg.gt.1 .and. smax.lt.smax1) cycle 
               f1=f0+real(idfbest)
               if( f1.le.10.0 .or. f1.ge.4990.0 ) cycle
               call timer('ft2ldown ',0)
               call ft2l_downsample(dd,dobigfft,f1,cb) !Final downsample, corrected f0
               call timer('ft2ldown ',1)
               sum2=sum(abs(cb)**2)/(real(NSS)*NN)
               if(sum2.gt.0.0) cb=cb/sqrt(sum2)
               cd=0.
               if(ibest.ge.0) then
                  it=min(NDMAX-1,ibest+NN*NSS-1)
                  np=it-ibest+1
                  cd(0:np-1)=cb(ibest:it)
               else
                  cd(-ibest:ibest+NN*NSS-1)=cb(0:NN*NSS+2*ibest-1)
               endif
               call timer('bitmet  ',0)
               call get_ft2l_bitmetrics(cd,bitmetrics,badsync)
               call timer('bitmet  ',1)
               if(badsync) then
                  cycle
               endif
               hbits=0
               where(bitmetrics(:,1).ge.0) hbits=1
               ns1=count(hbits(  1:  8).eq.(/0,0,0,1,1,0,1,1/))
               ns2=count(hbits( 67: 74).eq.(/0,1,0,0,1,1,1,0/))
               ns3=count(hbits(133:140).eq.(/1,1,1,0,0,1,0,0/))
               ns4=count(hbits(199:206).eq.(/1,0,1,1,0,0,0,1/))
               nsync_qual=ns1+ns2+ns3+ns4
               if(nsync_qual.lt. 15) cycle

               scalefac=2.83
               llra(  1: 58)=bitmetrics(  9: 66, 1)
               llra( 59:116)=bitmetrics( 75:132, 1)
               llra(117:174)=bitmetrics(141:198, 1)
               llra=scalefac*llra
               llrb(  1: 58)=bitmetrics(  9: 66, 2)
               llrb( 59:116)=bitmetrics( 75:132, 2)
               llrb(117:174)=bitmetrics(141:198, 2)
               llrb=scalefac*llrb
               llrc(  1: 58)=bitmetrics(  9: 66, 3)
               llrc( 59:116)=bitmetrics( 75:132, 3)
               llrc(117:174)=bitmetrics(141:198, 3)
               llrc=scalefac*llrc
               llrd(  1: 58)=bitmetrics(  9: 66, 4)
               llrd( 59:116)=bitmetrics( 75:132, 4)
               llrd(117:174)=bitmetrics(141:198, 4)
               llrd=scalefac*llrd
               llre(  1: 58)=bitmetrics(  9: 66, 5)
               llre( 59:116)=bitmetrics( 75:132, 5)
               llre(117:174)=bitmetrics(141:198, 5)
               llre=scalefac*llre

               apmag=maxval(abs(llra))*1.1
               npasses=5+nappasses(nQSOProgress)
               if(lapcqonly) npasses=6
               if(ndepth.eq.1) npasses=5
               if(ncontest.ge.6) npasses=5  ! Don't support Fox and Hound
               do ipass=1,npasses
                  if(ipass.eq.1) llr=llra
                  if(ipass.eq.2) llr=llrb
                  if(ipass.eq.3) llr=llrc
                  if(ipass.eq.4) llr=llrd
                  if(ipass.eq.5) llr=llre
                  if(ipass.le.5) then
                     apmask=0
                     iaptype=0
                  endif

                  if(ipass .gt. 5) then
                     llrd=llrc
                     iaptype=naptypes(nQSOProgress,ipass-5)
                     if(lapcqonly) iaptype=1

                     napwid=75
                     if(ncontest.le.5 .and. iaptype.ge.3 .and. (abs(f1-nfqso).gt.napwid) ) cycle
                     if(iaptype.ge.2 .and. apbits(1).gt.1) cycle  ! No, or nonstandard, mycall
                     if(iaptype.ge.3 .and. apbits(30).gt.1) cycle ! No, or nonstandard, dxcall

                     if(iaptype.eq.1) then  ! CQ or CQ TEST or CQ FD or CQ RU or CQ WW 
                        apmask=0
                        apmask(1:29)=1
                        if( ncontest.eq.0 ) llrd(1:29)=apmag*mcq(1:29)
                        if( ncontest.eq.1 ) llrd(1:29)=apmag*mcqtest(1:29)
                        if( ncontest.eq.2 ) llrd(1:29)=apmag*mcqtest(1:29)
                        if( ncontest.eq.3 ) llrd(1:29)=apmag*mcqfd(1:29)
                        if( ncontest.eq.4 ) llrd(1:29)=apmag*mcqru(1:29) 
                        if( ncontest.eq.5 ) llrd(1:29)=apmag*mcqww(1:29) 
                     endif

                     if(iaptype.eq.2) then ! MyCall,???,???
                        apmask=0
                        if(ncontest.eq.0.or.ncontest.eq.1.or.ncontest.eq.5) then
                           apmask(1:29)=1
                           llrd(1:29)=apmag*apbits(1:29)
                        else if(ncontest.eq.2) then
                           apmask(1:28)=1
                           llrd(1:28)=apmag*apbits(1:28)
                        else if(ncontest.eq.3) then
                           apmask(1:28)=1
                           llrd(1:28)=apmag*apbits(1:28)
                        else if(ncontest.eq.4) then
                           apmask(2:29)=1
                           llrd(2:29)=apmag*apmy_ru(1:28)
                        endif
                     endif

                     if(iaptype.eq.3) then ! MyCall,DxCall,???
                        apmask=0
                        if(ncontest.eq.0.or.ncontest.eq.1.or.ncontest.eq.2.or.ncontest.eq.5) then
                           apmask(1:58)=1
                           llrd(1:58)=apmag*apbits(1:58)
                        else if(ncontest.eq.3) then ! Field Day
                           apmask(1:56)=1
                           llrd(1:28)=apmag*apbits(1:28)
                           llrd(29:56)=apmag*aphis_fd(1:28)
                        else if(ncontest.eq.4) then 
                           apmask(2:57)=1
                           llrd(2:29)=apmag*apmy_ru(1:28)
                           llrd(30:57)=apmag*apbits(30:57)
                        endif
                     endif

                     if(iaptype.eq.4 .or. iaptype.eq.5 .or. iaptype.eq.6) then
                        apmask=0
                        if(ncontest.le.5) then
                           apmask(1:116)=1
                           llrd(1:116)=apmag*apbits(1:116)
                           apmask(117:135)=1
                           if(iaptype.eq.4) llrd(117:135)=apmag*mrrr(1:19)
                           if(iaptype.eq.5) llrd(117:135)=apmag*m73(1:19)
                           if(iaptype.eq.6) llrd(117:135)=apmag*mrr73(1:19)
                        endif
                     endif

                     llr=llrd
                  endif
                  message77=0
                  dmin=0.0

                  ndeep=3
                  maxosd=3  
                  if(abs(nfqso-f1).le.napwid) then
                     ndeep=3
                     maxosd=4
                  endif
                  if(.not.doosd) maxosd = -1
                  call timer('dec174_91 ',0)
                  Keff=91
                  call decode174_91(llr,Keff,maxosd,ndeep,apmask,message91,cw, &
                                    ntype,nharderror,dmin)
                  message77=message91(1:77)
                  call timer('dec174_91 ',1)

                  if(sum(message77).eq.0) cycle
                  if( nharderror.ge.0 ) then
                     message77=mod(message77+rvec,2) ! remove rvec scrambling
                     write(c77,'(77i1)') message77(1:77)
                     call unpack77(c77,1,message,unpk77_success)
                     if(.not.unpk77_success) exit
                     if(dosubtract) then
                        call get_ft2l_tones_from_77bits(message77,i4tone)
                        dt=real(ibest)/1333.33
                        call timer('subtract',0)
                        call subtractft2l(dd,i4tone,f1,dt)
                        call timer('subtract',1)
                     endif
                     idupe=0
                     do i=1,ndecodes
                        if(decodes(i).eq.message) idupe=1
                     enddo
                     if(idupe.eq.1) exit
                     ndecodes=ndecodes+1
                     decodes(ndecodes)=message
                     if(snr.gt.0.0) then
                        xsnr=10*log10(snr)-14.8
                     else
                        xsnr=-21.0
                     endif
                     nsnr=nint(max(-21.0,xsnr))
                     xdt=ibest/1333.33 - 0.5
                     qual=1.0-(nharderror+dmin)/60.0 
                     call this%callback(smax,nsnr,xdt,f1,message,iaptype,qual)
                     exit
                  endif
               enddo                      !Sequence estimation
               if(nharderror.ge.0) exit
            enddo                         !3 DT segments
         enddo                            !Candidate list
      enddo                               !Subtraction loop
      return
   end subroutine decode

end module ft2l_decode
