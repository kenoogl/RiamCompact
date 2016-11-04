subroutine copy_s2vex(MX, MY, MZ, NX, NY, NZ, S1, S2, S3, V3DEX);
!
!  use params
  implicit none
!
  integer i, j, k
  integer MX, MY, MZ
  integer NX, NY, NZ
  real(8) S1(0:MX+1,0:MY+1,0:MZ+1), S2(0:MX+1,0:MY+1,0:MZ+1), S3(0:MX+1,0:MY+1,0:MZ+1)
  real(8) V3DEX(3,0:MX+1,0:MY+1,0:MZ+1)
!
  do k=0,NZ+1
  do j=0,NY+1
  do i=0,NX+1
    V3DEX(1,i,j,k) = S1(i,j,k)
    V3DEX(2,i,j,k) = S2(i,j,k)
    V3DEX(3,i,j,k) = S3(i,j,k)
  end do
  end do
  end do
!
  return
end subroutine copy_s2vex
