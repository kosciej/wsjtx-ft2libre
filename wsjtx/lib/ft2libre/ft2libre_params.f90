! FT2Libre: FT8 modulation (8-GFSK, LDPC(174,91), Costas 7x7) at 2x FT4 baud rate
! Each message (79 symbols) is transmitted twice per T/R period (~3.8s)
! LDPC (174,91) code
parameter (KK=91)                     !Information bits (77 + CRC14)
parameter (ND=58)                     !Data symbols
parameter (NS=21)                     !Sync symbols (3 @ Costas 7x7)
parameter (NN=NS+ND)                  !Total channel symbols per copy (79)
parameter (NSPS=288)                  !Samples per symbol at 12000 S/s (24ms)
parameter (NZ=NSPS*NN*2)              !Samples in full waveform (45504, 2 copies)
parameter (NMAX=45000)                !Samples in iwave (3.75s = 60/16, divisible by 9)
parameter (NFFT1=2*NSPS, NH1=NFFT1/2) !Length of FFTs for symbol spectra (576, 288)
parameter (NSTEP=NSPS/4)              !Rough time-sync step size (72)
parameter (NHSYM=NMAX/NSTEP-3)        !Number of symbol spectra (622)
parameter (NDOWN=9)                   !Downsample factor (288/9=32 samp/sym)
