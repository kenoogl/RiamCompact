subroutine calc_contravariant_v(MX, MY, MZ, NX, NY, NZ, DT, U, V, W, GX, EY, TX, TY, TZ, &
                                UU, VV, WW, YJA, RHS)
!
!  use params
  implicit none
!
  integer i, j, k
  integer MX, MY, MZ
  integer NX, NY, NZ
!
  real(8) DT
  real(8) U(0:MX+1,0:MY+1,0:MZ+1), V(0:MX+1,0:MY+1,0:MZ+1), W(0:MX+1,0:MY+1,0:MZ+1)
  real(8) UU(0:MX+1,0:MY+1,0:MZ+1), VV(0:MX+1,0:MY+1,0:MZ+1), WW(0:MX+1,0:MY+1,0:MZ+1)
  real(8) CU(0:MX+1,0:MY+1,0:MZ+1), CV(0:MX+1,0:MY+1,0:MZ+1), CW(0:MX+1,0:MY+1,0:MZ+1)
  real(8) GX(MX,MY,MZ), EY(MX,MY,MZ), TX(MX,MY,MZ), TY(MX,MY,MZ), TZ(MX,MY,MZ)
  real(8) YJA(MX,MY,MZ), RHS(MX,MY,MZ)
!
!
!==========================================================
!....<<NEW PRESSURE, PHYSICAL & CONTRAVARIANT VELOCITY>>
!==========================================================
!....<<CONTRAVARIANT VELOCITY OF A CELL CENTER>>
  do k = 1, NZ
    do j = 1, NY
      do i = 1, NX
        CU(i,j,k) = YJA(i,j,k) * U(i,j,k) * GX(i,j,k)
        CV(i,j,k) = YJA(i,j,k) * V (i,j,k) * EY(i,j,k)
        CW(i,j,k) = YJA(i,j,k) * (U(i,j,k) * TX(i,j,k) + V(i,j,k) * TY(i,j,k) + W(i,j,k) * TZ(i,j,k))
      enddo
    enddo
  enddo
!
!....<<CONTRAVARIANT VELOCITY OF A CELL SURFACE>>
  do k = 2, NZ-1
    do j = 2, NY-1
      do i = 1, NX-1
        UU(i,j,k) = .5D0 * (CU(i+1,j,k) + CU(i,j,k))
      enddo
    enddo
  enddo
!
  do k= 2, NZ-1
    do j= 1, NY-1
      do i= 2, NX-1
        VV(i,j,k) = .5D0 * ( CV(i,j+1, k) + CV(i,j,k))
      enddo
    enddo
  enddo
!
  do k = 1, NZ-1
    do j = 2, NY-1
      do i = 2, NX-1
        WW(i,j,k) = .5D0 * (CW(i,j,k+1) + CW(i,j,k))
      enddo
    enddo
  enddo
!
!....<<RIGHT HAND SIDE>>
  do k = 2, NZ-1
    do j = 2, NY-1
      do i = 2, NX-1
        RHS(i,j,k) = -((UU(i,j,k) - UU(i-1,j,k) &
                    + VV(i,j,k)   - VV(i,j-1,k) &
                    + WW(i,j,k)   -  WW(i,j,k-1)) / YJA(i,j,k)) / DT
      enddo
    enddo
  enddo
!
  return
end subroutine calc_contravariant_v
