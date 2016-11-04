subroutine display2(MX, MY, MZ, NX, NY, NZ, TIME, RE, U, V, W, P, UU, VV, WW, &
                    MY_RANK)
!
!  use params
  implicit none
!
  integer i, j, k
  integer MX, MY, MZ
  integer NX, NY, NZ
  integer MY_RANK
  real(8) TIME, RE
  real(8) U(0:MX+1,0:MY+1,0:MZ+1), V(0:MX+1,0:MY+1,0:MZ+1), W(0:MX+1,0:MY+1,0:MZ+1)
  real(8) P(0:MX+1,0:MY+1,0:MZ+1)
  real(8) UU(0:MX+1,0:MY+1,0:MZ+1), VV(0:MX+1,0:MY+1,0:MZ+1), WW(0:MX+1,0:MY+1,0:MZ+1)
!
  character(len=256) fname
  write(fname,'(a,i6.6,a)') 'uvwp_id',MY_RANK,'.dat'

! /....<<FIELD-DATA OUTPUT>>
  open(1,FILE=fname,FORM='FORMATTED',STATUS='UNKNOWN')
  write(1,'(3I5,2F15.7)')NX,NY,NZ,TIME,RE
  write(1,'(3F15.7)')(((U(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
  write(1,'(3F15.7)')(((V(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
  write(1,'(3F15.7)')(((W(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
  write(1,'(3F15.7)')(((P(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
  write(1,'(3F15.7)')(((UU(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
  write(1,'(3F15.7)')(((VV(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
  write(1,'(3F15.7)')(((WW(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
  close(1)
!
  return
end subroutine display2
