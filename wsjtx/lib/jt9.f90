program jt9

! Decoder for JT9.  Can run stand-alone, reading data from *.wav files;
! or as the back end of wsjt-x, with data placed in a shared memory region.

  use options
  use prog_args
  use, intrinsic :: iso_c_binding
  use FFTW3
  use timer_module, only: timer
  use timer_impl, only: init_timer, fini_timer
  use readwav
  use ft8_mod1, only : dd8
  use jt65_mod6, only : dd

  include 'jt9com.f90'

  integer(C_INT) iret
  type(wav_header) wav
  real*8 TRperiod
  character c
  character(len=500) optarg, infile
  character wisfile*256

  integer :: arglen,stat,offset,remain,mode=0,flow=200,fsplit=2700,          &
       fhigh=4000,nrxfreq=1500,ndepth=1,nexp_decode=0,nQSOProg=0,ncycles=3,  &
       nft8rxfsens=3,nmt=0,nmtft8decsens=3,ndecoderstart=3
  logical :: read_files = .true., tx9 = .false., display_help = .false.,     &
       nexp_decode_set = .false., have_ntol = .false., multift8 = .false.,   &
       hidedupes = .false.
  type (option) :: long_options(42) = [                                      &
    option ('help', .false., 'h', 'Display this help message', ''),          &
    option ('shmem',.true.,'s','Use shared memory for sample data','KEY'),   &
    option ('tr-period', .true., 'p', 'Tx/Rx period, default SECONDS=60',    &
        'SECONDS'),                                                          &
    option ('executable-path', .true., 'e',                                  &
        'Location of subordinate executables (KVASD) default PATH="."',      &
        'PATH'),                                                             &
    option ('data-path', .true., 'a',                                        &
        'Location of writeable data files, default PATH="."', 'PATH'),       &
    option ('temp-path', .true., 't',                                        &
        'Temporary files path, default PATH="."', 'PATH'),                   &
    option ('lowest', .true., 'L',                                           &
        'Lowest frequency decoded (JT65), default HERTZ=200', 'HERTZ'),      &
    option ('highest', .true., 'H',                                          &
        'Highest frequency decoded, default HERTZ=4000', 'HERTZ'),           &
    option ('split', .true., 'S',                                            &
        'Lowest JT9 frequency decoded, default HERTZ=2700', 'HERTZ'),        &
    option ('rx-frequency', .true., 'f',                                     &
        'Receive frequency offset, default HERTZ=1500', 'HERTZ'),            &
    option ('freq-tolerance', .true., 'F',                                   &
        'Receive frequency tolerance, default HERTZ=20', 'HERTZ'),           &
    option ('patience', .true., 'w',                                         &
        'FFTW3 planing patience (0-4), default PATIENCE=1', 'PATIENCE'),     &
    option ('fft-threads', .true., 'm',                                      &
        'Number of threads to process large FFTs, default THREADS=1',        &
        'THREADS'),                                                          &
    option ('multithreadft8', .false., 'M', 'Use Multithread FT8 Decoder',   &
        ''),                                                                 &
    option ('MTft8-cycles', .true., 'C',                                     &
        'Number of cycles for MTft8 decoder (1-3), default CYCLES=3',        &
        'CYCLES'),                                                           &
    option ('MTft8-hidedupes', .false., 'U', 'Hide MTft8 Dupes',             &
        ''),                                                                 &
    option ('MTft8-RxFreqSens', .true., 'R',                                 &
        'RxFreqSensitivity for MTft8 decoder (1-3), default RXFREQSENS=3',   &
        'RXFREQSENS'),                                                       &
    option ('MTft8-NumThreads', .true., 'N',                                 &
        'Number of threads for MTft8 decoder (0-12), default NUMTHREADS=0',  &
        'NUMTHREADS'),                                                       &
    option ('MTft8-DecSens', .true., 'E',                                    &
        'MTft8 decoder sensitivity (1-3), default SENSITIVITY=3',            &
        'SENSITIVITY'),                                                      &
    option ('MTft8-StartTime', .true., 'D',                                  &
        'MTft8 decoder start (0-4), default START_TIME=3 (Normal)',          &
        'START_TIME'),                                                       &
    option ('Skip-MTft8-WideDxCallSearch', .false., 'Z',                     & 
        'SKIP MTft8 Wideband DX Call Search', ''),                           &        
    option ('q65', .false., '3', 'Q65 mode', ''),                            &
    option ('jt4', .false., '4', 'JT4 mode', ''),                            &
    option ('ft2l', .false., '2', 'FT2L mode', ''),                          &
    option ('ft4', .false., '5', 'FT4 mode', ''),                            &
    option ('jt65', .false.,'6', 'JT65 mode', ''),                           &
    option ('fst4', .false., '7', 'FST4 mode', ''),                          &
    option ('fst4w', .false., 'W', 'FST4W mode', ''),                        &
    option ('fst4w', .false., 'Y', 'FST4W mode, print hash22 values', ''),   &
    option ('ft8', .false., '8', 'FT8 mode', ''),                            &
    option ('jt9', .false., '9', 'JT9 mode', ''),                            &
    option ('quiet', .false., 'q', 'Quiet mode, no <DecodeFinished>', ''),   &
    option ('msk144', .false., 'k', 'MSK144 mode', ''),                      &
    option ('QSOprog', .true., 'Q', 'QSO progress (0-5), default PROGRESS=1',&
        'QSOprogress'),                                                      &
    option ('sub-mode', .true., 'b', 'Sub mode, default SUBMODE=A', 'A'),    &
    option ('depth', .true., 'd',                                            &
        'Decoding depth (1-3), default DEPTH=1', 'DEPTH'),                   &
    option ('tx-jt9', .false., 'T', 'Tx mode is JT9', ''),                   &
    option ('my-call', .true., 'c', 'my callsign', 'CALL'),                  &
    option ('my-grid', .true., 'G', 'my grid locator', 'GRID'),              &
    option ('his-call', .true., 'x', 'his callsign', 'CALL'),                &
    option ('his-grid', .true., 'g', 'his grid locator', 'GRID'),            &
    option ('experience-decode', .true., 'X',                                &
        'experience based decoding flags (1..n), default FLAGS=0',           &
        'FLAGS') ]

  type(dec_data), allocatable :: shared_data
  character(len=12) :: mycall='K1ABC', hiscall='W9XYZ'
  character(len=6) :: mygrid='', hisgrid='EN37'
  common/patience/npatience,nthreads
  common/decstats/ntry65a,ntry65b,n65a,n65b,num9,numfano
  data npatience/1/,nthreads/1/,wisfile/' '/

  nsubmode = 0
  ntol = 20
  TRperiod=60.d0

  do
     call getopt('hs:e:a:b:r:m:p:d:f:F:w:t:98765432WYqkTMUSZL:S:H:c:G:x:g:X:Q:C:R:N:E:D:',     &
          long_options,c,optarg,arglen,stat,offset,remain,.true.)
     if (stat .ne. 0) then
        exit
     end if
     select case (c)
        case ('h')
           display_help = .true.
        case ('s')
           read_files = .false.
           shm_key = optarg(:arglen)
        case ('e')
           exe_dir = optarg(:arglen)
        case ('a')
           data_dir = optarg(:arglen)
        case ('b')
           nsubmode = ichar (optarg(:1)) - ichar ('A')
        case ('t')
           temp_dir = optarg(:arglen)
        case ('m')
           read (optarg(:arglen), *) nthreads
        case ('p')
           read (optarg(:arglen), *) TRperiod
        case ('d')
           read (optarg(:arglen), *) ndepth
        case ('f')
           read (optarg(:arglen), *) nrxfreq
        case ('F')
           read (optarg(:arglen), *) ntol
           have_ntol = .true.
        case ('L')
           read (optarg(:arglen), *) flow
        case ('S')
           read (optarg(:arglen), *) fsplit
        case ('H')
           read (optarg(:arglen), *) fhigh
        case ('M')
           multift8 = .true.
        case ('C')
           read (optarg(:arglen), *) ncycles
        case ('U')
           hidedupes = .true.
        case ('R')
           read (optarg(:arglen), *) nft8rxfsens
        case ('N')
           read (optarg(:arglen), *) nmt
        case ('E')
           read (optarg(:arglen), *) nmtft8decsens
        case ('D')
           read (optarg(:arglen), *) ndecoderstart
        case ('Z')
           lwidedxcsearch = .false.
        case ('q')
           lquiet=.true.
        case ('k')
           mode = 144
        case ('Q')
           read (optarg(:arglen), *) nQSOProg
        case ('3')
           mode = 66
        case ('4')
           mode = 4
        case ('2')
           mode = 6
        case ('5')
           mode = 5
        case ('6')
           if (mode.lt.65) mode = mode + 65
        case ('7')
           mode = 240
        case ('8')
           mode = 8
        case ('9')
           if (mode.lt.9.or.mode.eq.65) mode = mode + 9
        case ('T')
           tx9 = .true.
        case ('w')
           read (optarg(:arglen), *) npatience
        case ('W')
           mode = 241
        case ('Y')
           mode = 242
        case ('c')
           read (optarg(:arglen), *) mycall
        case ('G')
           read (optarg(:arglen), *) mygrid
        case ('x')
           read (optarg(:arglen), *) hiscall
        case ('g')
           read (optarg(:arglen), *) hisgrid
        case ('X')
           read (optarg(:arglen), *) nexp_decode
           nexp_decode_set = .true.
     end select
  end do
  
  if (display_help .or. stat .lt. 0                      &
       .or. (.not. read_files .and. remain .gt. 0)       &
       .or. (read_files .and. remain .lt. 1)) then

     print *, 'Usage: jt9 [OPTIONS] file1 [file2 ...]'
     print *, '       Reads data from *.wav files.'
     print *, ''
     print *, '       jt9 -s <key> [-w patience] [-m threads] [-e path] [-a path] [-t path]'
     print *, '       Gets data from shared memory region with key==<key>'
     print *, ''
     print *, 'OPTIONS:'
     print *, ''
     do i = 1, size (long_options)
       call long_options(i) % print (6)
     end do
     go to 999
  endif

  iret=fftwf_init_threads()            !Initialize FFTW threading 

! Default to 1 thread, but use nthreads for the big ones
  call fftwf_plan_with_nthreads(1)

! Import FFTW wisdom, if available
  wisfile=trim(data_dir)//'/jt9_wisdom.dat'// C_NULL_CHAR
  iret=fftwf_import_wisdom_from_filename(wisfile)

  ntry65a=0
  ntry65b=0
  n65a=0
  n65b=0
  num9=0
  numfano=0

  if (.not. read_files) then
     call jt9a()          !We're running under control of WSJT-X
     go to 999
  endif

  if(mycall.eq.'b') mycall='            '
  if(hiscall.eq.'b') then
     hiscall='            '
     hisgrid='      '
  endif

  if (mode .eq. 241 .or. mode .eq. 242) then
     ntol = min (ntol, 100)
  else if (mode .eq. 65 + 9 .and. .not. have_ntol) then
     ntol = 20
  else if (mode .eq. 66 .and. .not. have_ntol) then
     ntol = 10
  else
     ntol = min (ntol, 1000)
  end if

  allocate(shared_data)
  nflatten=0
  do iarg = offset + 1, offset + remain
     call get_command_argument (iarg, optarg, arglen)
     infile = optarg(:arglen)
     call wav%read (infile)
     nfsample=wav%audio_format%sample_rate
     i1=index(infile,'.wav')
     if(i1.lt.1) i1=index(infile,'.WAV')
     if(infile(i1-5:i1-5).eq.'_') then
        read(infile(i1-4:i1-1),*,err=1) nutc
     else
        read(infile(i1-6:i1-1),*,err=1) nutc
     endif
     go to 2
1    nutc=0
2    nsps=6912
     npts=nint(TRperiod*12000.0)
     kstep=nsps/2
     
     do while (.true.)
        k=0
        nhsym=0
        if(iarg .eq. offset + 1 .and. nutc.eq.0) then
           call init_timer (trim(data_dir)//'/timer.out')
           call timer('jt9     ',0)
        endif
        shared_data%id2=0
        
        do iblk=1,npts/kstep
           k=iblk*kstep
           read(unit=wav%lun,end=3) shared_data%id2(k-kstep+1:k)
           go to 4
3          if(iblk.eq.1 .and. k.eq.0) goto 5 ! EOF at start of block
           shared_data%id2(k+1:)=0
           exit
4          nhsym=(k-2048)/kstep
        enddo

        shared_data%params%nutc=nutc
        shared_data%params%ndiskdat=.true.
        shared_data%params%ntr=TRperiod
        shared_data%params%nfqso=nrxfreq
        shared_data%params%newdat=.true.
        shared_data%params%npts8=74736
        shared_data%params%nfa=flow
        shared_data%params%nfsplit=fsplit
        shared_data%params%nfb=fhigh
        shared_data%params%ntol=ntol
        shared_data%params%kin=64800
        shared_data%params%nzhsym=nhsym
        shared_data%params%ndepth=ndepth
        shared_data%params%lft8apon=.true.
        shared_data%params%ljt65apon=.true.
        shared_data%params%napwid=75
        shared_data%params%dttol=3.
        shared_data%params%nmode=mode
        if(mode.eq.0) shared_data%params%nmode=65+9
        shared_data%params%nsubmode=nsubmode
        shared_data%params%mycall=transfer(mycall,shared_data%params%mycall)
        shared_data%params%mygrid=transfer(mygrid,shared_data%params%mygrid)
        shared_data%params%hiscall=transfer(hiscall,shared_data%params%hiscall)
        shared_data%params%hisgrid=transfer(hisgrid,shared_data%params%hisgrid)
        
        call multimode_decoder(shared_data%ss,shared_data%id2, &
             shared_data%params,nfsample)
        
        nutc = nutc + nint(TRperiod/60.0) ! Dummy update
        if (k.lt.npts) exit ! Last block
     enddo
5    close(unit=wav%lun)
  enddo

  call timer('jt9     ',1)
  call timer('jt9     ',101)

999 continue

! Output decoder statistics
  call fini_timer ()

! Save FFTW wisdom and free memory
  if(len(trim(wisfile)).gt.0) iret=fftwf_export_wisdom_to_filename(wisfile)
  call four2a(a,-1,1,1,1)
  call filbig(a,-1,1,0.0,0,0,0,0,0)        !used for all other FFT plans
  call fftwf_cleanup_threads()
  call fftwf_cleanup()

end program jt9
