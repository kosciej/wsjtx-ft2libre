module ft2libre_decode

  type :: ft2libre_decoder
     procedure(ft2libre_decode_callback), pointer :: callback
   contains
     procedure :: decode
  end type ft2libre_decoder

  abstract interface
     subroutine ft2libre_decode_callback (this,sync,snr,dt,freq,decoded,nap,qual)
       import ft2libre_decoder
       implicit none
       class(ft2libre_decoder), intent(inout) :: this
       real, intent(in) :: sync
       integer, intent(in) :: snr
       real, intent(in) :: dt
       real, intent(in) :: freq
       character(len=37), intent(in) :: decoded
       integer, intent(in) :: nap
       real, intent(in) :: qual
     end subroutine ft2libre_decode_callback
  end interface

contains

  subroutine decode(this,callback,iwave,nfqso,newdat,nutc,nfa,nfb,   &
       nzhsym,ndepth,ncontest,nagain,mycall12,hiscall12,              &
       nQSOProgress,nftx,lapon,lapcqonly,napwid)

    use timer_module, only: timer

    include 'ft2libre/ft2libre_params.f90'

    class(ft2libre_decoder), intent(inout) :: this
    procedure(ft2libre_decode_callback) :: callback
    parameter (MAXCAND=200,MAX_DECODES=100,NPTS=NMAX)
    real sbase(NH1)
    real candidate(3,MAXCAND)
    real dd(NPTS)
    logical, intent(in) :: nagain
    logical newdat,lsubtract,lapon,lapcqonly
    character*12 mycall12,hiscall12
    integer*2 iwave(NPTS)
    integer nQSOProgress,nftx,napwid
    character*37 msg37
    character*37 allmessages(MAX_DECODES)
    integer allsnrs(MAX_DECODES)
    integer itone(NN)
    logical ldupe
    real rms_dd

    save ndecodes,allmessages,allsnrs

    this%callback => callback

    dd=iwave
    ndecodes=0
    allmessages='                                     '
    allsnrs=0

    open(73,file='/tmp/ft2l_debug.log',status='unknown',              &
         position='append',action='write')
    rms_dd=sqrt(sum(dd**2)/NPTS)
    write(73,'(A,I6.6,A,F10.1,A,I6,A,I6)') '=== FT2L decode UTC=',  &
         nutc,' rms=',rms_dd,' nfa=',nfa,' nfb=',nfb
    flush(73)

    ifa=nfa
    ifb=nfb
    if(nagain) then
       ifa=nfqso-20
       ifb=nfqso+20
    endif

    npass=3
    imetric=1
    if(ndepth.eq.1) npass=2
    do ipass=1,npass
      newdat=.true.
      syncmin=1.3
      if(ndepth.le.2) syncmin=2.1
      if(ipass.eq.1) then
        lsubtract=.true.
        imetric=1
      elseif(ipass.eq.2) then
        imetric=2
        lsubtract=.true.
      elseif(ipass.eq.3) then
        imetric=2
        if(ndecodes.eq.0) cycle
        lsubtract=.true.
      endif
      write(73,'(A,I2,A,I4)') 'FT2L: pass_start ipass=',ipass,        &
           ' ndecodes=',ndecodes
      flush(73)
      call timer('syn_f2l ',0)
      maxc=MAXCAND
      call sync_ft2libre(dd,NPTS,ifa,ifb,syncmin,nfqso,maxc,         &
           candidate,ncand,sbase)
      call timer('syn_f2l ',1)
      write(73,'(A,I2,A,I4)') 'FT2L: post_sync ipass=',ipass,        &
           ' ndecodes=',ndecodes
      write(73,'(A,I2,A,I2,A,I4,A,I6,A,I6)') 'FT2L: pass=',ipass,  &
           ' metric=',imetric,' ncand=',ncand,' ifa=',ifa,' ifb=',ifb
      do icand=1,min(ncand,10)
        write(73,'(A,I3,A,F8.1,A,F8.3,A,F8.3)') 'FT2L: cand',icand, &
             ' f1=',candidate(1,icand),' xdt=',candidate(2,icand),   &
             ' sync=',candidate(3,icand)
      enddo
      flush(73)
      do icand=1,ncand
        sync=candidate(3,icand)
        f1=candidate(1,icand)
        xdt=candidate(2,icand)
        xbase=10.0**(0.1*(sbase(nint(f1/(12000.0/NFFT1)))-40.0))
        msg37='                                     '
        iaptype=0
        call timer('f2libreb',0)
        call ft2libreb(dd,newdat,nfqso,ndepth,nzhsym,lsubtract,      &
             nagain,ncontest,imetric,f1,xdt,xbase,nharderrors,dmin,   &
             nbadcrc,iaptype,msg37,xsnr,itone,mycall12,hiscall12,    &
             nQSOProgress,nftx,lapon,lapcqonly,napwid)
        call timer('f2libreb',1)
        nsnr=nint(xsnr)
        if(nbadcrc.eq.0) then
           ldupe=.false.
           do id=1,ndecodes
              if(msg37.eq.allmessages(id)) ldupe=.true.
           enddo
           write(73,'(A,I2,A,I3,A,L1,A,I3,A,A37)') 'FT2L: DECODED pass=', &
                ipass,' cand=',icand,' dupe=',ldupe,' ap=',iaptype,        &
                ' msg=',msg37
           if(.not.ldupe .and. ndecodes.gt.0) then
              write(73,'(A,I4)') 'FT2L: DUPE_MISS ndecodes=',ndecodes
              write(73,'(A,20I4)') 'FT2L: msg   =',                       &
                   (ichar(msg37(k:k)),k=1,20)
              write(73,'(A,20I4)') 'FT2L: stored=',                        &
                   (ichar(allmessages(1)(k:k)),k=1,20)
              flush(73)
           endif
           if(.not.ldupe) then
              if(ndecodes.ge.MAX_DECODES) cycle
              ndecodes=ndecodes+1
              allmessages(ndecodes)=msg37
              allsnrs(ndecodes)=nsnr
           endif
           if(.not.ldupe .and. associated(this%callback)) then
              qual=1.0-(nharderrors+dmin)/60.0
              call this%callback(sync,nsnr,xdt,f1,msg37,iaptype,qual)
              write(73,'(A)') 'FT2L: callback fired'
           endif
           flush(73)
        endif
      enddo  ! icand
    enddo  ! ipass

    write(73,'(A,I4)') 'FT2L: total decodes=',ndecodes
    flush(73)
    close(73)
    return
  end subroutine decode

end module ft2libre_decode
