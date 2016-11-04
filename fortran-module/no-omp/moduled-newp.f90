subroutine calc_new_pressure(MX, MY, MZ, NX, NY, NZ, REI, P, Z, W, C1, C2, C3, C5, C6, C7, C8, C9)
!
!  use params
  implicit none
!
  integer i, j, k
  integer MX, MY, MZ
  integer NX, NY, NZ

  real(8) REI, PP
  real(8) P(0:MX+1,0:MY+1,0:MZ+1),Z(0:MX+1,0:MY+1,0:MZ+1), W(0:MX+1,0:MY+1,0:MZ+1)
  real(8) C1(MX,MY,MZ), C2(MX,MY,MZ), C3(MX,MY,MZ), C5(MX,MY,MZ), C6(MX,MY,MZ)
  real(8) C7(MX,MY,MZ), C8(MX,MY,MZ), C9(MX,MY,MZ)


!....<<AVERAGING & SUBTRACTION>>
  PP = 0.D0
  do k = 1, NZ
    do j = 1, NY
      do i = 1, NX
        PP = PP + P(i,j,k)
      enddo
    enddo
  enddo
!
  do k = 1, NZ
    do j = 1, NY
      do i = 1, NX
        P(i,j,k) = P(i,j,k) - PP / dfloat(NX*NY*NZ)
      enddo
    enddo
  enddo
!
!....<<MODIFY PRESSURE B.C.>>
  k = 1
  do j = 2, NY-1
    do i = 2, NX-1
      P(i,j,1) = P(i,j,2) - .5D0 * (Z(i,j,2) - Z(i,j,0)) &
               * (C1(i,j,k) * (W(i+1,j,k) - 2.D0 * W(i,j,k) + W(i-1,j,k)) &
               +  C2(i,j,k) * (W(i,j+1,k) - 2.D0 * W(i,j,k) + W(i,j-1,k)) &
               +  C3(i,j,k) * (W(i,j,k+1) - 2.D0 * W(i,j,k) + W(i,j,k-1)) &
               +  C5(i,j,k) * (W(i,j+1,k+1) - W(i,j+1,k-1) - W(i,j-1,k+1) + W(i,j-1,k-1)) * .25D0 &
               +  C6(i,j,k) * (W(i+1,j,k+1) - W(i+1,j,k-1) - W(i-1,j,k+1) + W(i-1,j,k-1)) * .25D0 &
               +  C7(i,j,k) * (W(i+1,j,k) - W(i-1,j,k)) * .5D0 &
               +  C8(i,j,k) * (W(i,j+1,k) - W(i,j-1,k)) * .5D0 &
               +  C9(i,j,k) * (W(i,j,k+1) - W(i,j,k-1)) * .5D0) * REI
      P(i,j,NZ) = P(i,j,NZ-1)
    enddo
  enddo
!
  do  k =1, NZ
    do i = 1, NX
      P(i,1,k) = P(i,2,k)
      P(i,NY,k) = P(i,NY-1,k)
    enddo
  enddo
!
  do k = 1, NZ
    do j = 1, NY
      P(1,j,k) = P(2,j,k)
      P(NX,j,k) = P(NX-1,j,k)
    enddo
  enddo
!
  return
end subroutine calc_new_pressure
