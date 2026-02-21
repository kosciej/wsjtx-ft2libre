subroutine genft2libre(msg,i3,n3,msgsent,msgbits,itone)

! Encode an FT2Libre message, producing array itone(NN=103).
! Uses 4-GFSK with LDPC(174,91), scrambling vector, 4x4 Costas arrays.

  use packjt77
  include 'ft2libre_params.f90'
  character msg*37,msgsent*37
  character*77 c77
  integer*1 msgbits(77),codeword(174),scrambled(77)
  integer itone(NN)
  integer icos_a(4),icos_b(4),icos_c(4),icos_d(4)
  integer graymap(0:3)
  logical unpk77_success

  data icos_a/0,1,3,2/
  data icos_b/1,0,2,3/
  data icos_c/2,3,1,0/
  data icos_d/3,2,0,1/
  data graymap/0,1,3,2/

  integer*1 rvec(77)
  data rvec/0,1,0,0,1,0,1,0,0,1,0,1,1,1,1,0,1,0,0,0,1,0,0,1,1,0,1,1,0, &
            1,0,0,1,0,1,1,0,0,0,0,1,0,0,0,1,0,1,0,0,1,1,1,1,0,0,1,0,1, &
            0,1,0,1,0,1,1,0,1,1,1,1,1,0,0,0,1,0,1/

  i3=-1
  n3=-1
  call pack77(msg,i3,n3,c77)
  call unpack77(c77,0,msgsent,unpk77_success)
  read(c77,'(77i1)',err=1) msgbits
  if(unpk77_success) go to 2
1 msgbits=0
  itone=0
  msgsent='*** bad message ***                  '
  go to 900

entry get_ft2libre_tones_from_77bits(msgbits,itone)

! Scramble message bits with rvec before encoding
2 do i=1,77
    scrambled(i)=ieor(msgbits(i),rvec(i))
  enddo
  call encode174_91(scrambled,codeword)

! Insert four 4-symbol Costas arrays
  itone(1:4)=icos_a
  itone(34:37)=icos_b
  itone(67:70)=icos_c
  itone(100:103)=icos_d

! Map 174 coded bits to 87 data symbols (2 bits/symbol, Gray code)
! Data positions: 5-33, 38-66, 71-99
  k=4
  do j=1,ND
     i=2*j-1
     k=k+1
     if(j.eq.30) k=k+4                   !Skip Costas B
     if(j.eq.59) k=k+4                   !Skip Costas C
     indx=codeword(i)*2+codeword(i+1)
     itone(k)=graymap(indx)
  enddo

900 return
end subroutine genft2libre
