! FT2Libre: FT8 modulation (8-GFSK, LDPC(174,91), Costas 7x7) at ~33 Bd
! Modulation index h=0.75: tone spacing = 25 Hz, BW ~ 175 Hz
! LDPC (174,91) code
parameter (KK=91)                     !Information bits (77 + CRC14)
parameter (ND=58)                     !Data symbols
parameter (NS=21)                     !Sync symbols (3 @ Costas 7x7)
parameter (NN=NS+ND)                  !Total channel symbols (79)
parameter (NSPS=360)                  !Samples per symbol at 12000 S/s (30ms)
parameter (NZ=NSPS*NN)                !Samples in full waveform (28440)
parameter (NMAX=45000)                !Samples in iwave (3.75s = 60/16)
parameter (NFFT1=960, NH1=NFFT1/2)    !FFTs for symbol spectra (df=12.5 Hz)
parameter (NSTEP=NSPS/4)              !Rough time-sync step size (90)
parameter (NHSYM=NMAX/NSTEP-3)        !Number of symbol spectra (497)
parameter (NDOWN=12)                  !Downsample factor (360/12=30 samp/sym)
parameter (NTBIN=2)                   !Tone spacing in FFT bins (h=0.75, nfos*h)
