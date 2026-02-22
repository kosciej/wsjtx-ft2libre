! FT2Libre: 4-GFSK, LDPC(174,91), 4x4 Costas arrays at ~41.7 Bd
! Modulation index h=1.0: tone spacing = 41.667 Hz, BW ~ 125 Hz
! LDPC (174,91) code
parameter (KK=91)                     !Information bits (77 + CRC14)
parameter (ND=87)                     !Data symbols (174 bits / 2 bits per symbol)
parameter (NS=16)                     !Sync symbols (4 x 4-symbol Costas arrays)
parameter (NN=NS+ND)                  !Total channel symbols (103)
parameter (NN2=NN+2)                  !Total with ramp symbols (105)
parameter (NSPS=288)                  !Samples per symbol at 12000 S/s
parameter (NZ=NSPS*NN)                !Samples in sync+data waveform (29664)
parameter (NZ2=NSPS*NN2)              !Total samples including ramps (30240)
parameter (NMAX=45000)                !Samples in iwave (3.75s = 60/16)
parameter (NFFT1=1152, NH1=NFFT1/2)  !FFTs for symbol spectra (df=10.42 Hz)
parameter (NSTEP=NSPS)                !Rough time-sync step size (288)
parameter (NHSYM=(NMAX-NFFT1)/NSTEP)  !Number of symbol spectra (152)
parameter (NDOWN=9)                   !Downsample factor (288/9=32 samp/sym)
parameter (NTBIN=4)                   !Tone spacing in FFT bins (h=1.0)
