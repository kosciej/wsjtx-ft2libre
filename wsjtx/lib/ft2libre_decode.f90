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
       nzhsym,ndepth,ncontest,nagain,mycall12,hiscall12)

    use timer_module, only: timer

    include 'ft2libre/ft2libre_params.f90'

    class(ft2libre_decoder), intent(inout) :: this
    procedure(ft2libre_decode_callback) :: callback
    parameter (MAXCAND=200,MAX_DECODES=100,NPTS=NMAX)
    real sbase(NH1)
    real candidate(3,MAXCAND)
    real dd(NPTS)
    logical, intent(in) :: nagain
    logical newdat,lsubtract
    character*12 mycall12,hiscall12
    integer*2 iwave(NPTS)
    character*37 msg37
    character*37 allmessages(MAX_DECODES)
    integer allsnrs(MAX_DECODES)
    integer itone(NN)
    logical ldupe

    this%callback => callback

    dd=iwave
    ndecodes=0
    allmessages='                                     '
    allsnrs=0

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
      call timer('syn_f2l ',0)
      maxc=MAXCAND
      call sync_ft2libre(dd,NPTS,ifa,ifb,syncmin,nfqso,maxc,         &
           candidate,ncand,sbase)
      call timer('syn_f2l ',1)
      do icand=1,ncand
        sync=candidate(3,icand)
        f1=candidate(1,icand)
        xdt=candidate(2,icand)
        xbase=10.0**(0.1*(sbase(nint(f1/(12000.0/NFFT1)))-40.0))
        msg37='                                     '
        call timer('f2libreb',0)
        call ft2libreb(dd,newdat,nfqso,ndepth,nzhsym,lsubtract,      &
             nagain,ncontest,imetric,f1,xdt,xbase,nharderrors,dmin,   &
             nbadcrc,iappass,msg37,xsnr,itone)
        call timer('f2libreb',1)
        nsnr=nint(xsnr)
        xdt=xdt-0.5
        if(nbadcrc.eq.0) then
           ldupe=.false.
           do id=1,ndecodes
              if(msg37.eq.allmessages(id)) ldupe=.true.
           enddo
           if(.not.ldupe) then
              if(ndecodes.ge.MAX_DECODES) cycle
              ndecodes=ndecodes+1
              allmessages(ndecodes)=msg37
              allsnrs(ndecodes)=nsnr
           endif
           if(.not.ldupe .and. associated(this%callback)) then
              iaptype=0
              qual=1.0-(nharderrors+dmin)/60.0
              call this%callback(sync,nsnr,xdt,f1,msg37,iaptype,qual)
           endif
        endif
      enddo  ! icand
    enddo  ! ipass

    return
  end subroutine decode

end module ft2libre_decode
