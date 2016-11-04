subroutine calc_new_velocity(MX, MY, MZ, NX, NY, NZ, KA, KB, JA, JB, IA, IB, DT, X, P, U, V, W, &
                             UD, VD, WD, UU, VV, WW, GX, EY, TX, TY, TZ, YJA, D1, D2, D3, D5, D6 )
!
!  use params
  implicit none
!
  integer i, j, k
  integer MX, MY, MZ
  integer NX, NY, NZ
  integer KA, KB, JA, JB, IA, IB 
!
  real(8) DT
  real(8) PG, PE, PT, PX, PY, PZ, DXI
  real(8) DUDX, DVDX, DWDX
!
  real(8) X(0:MX+1,0:MY+1,0:MZ+1), YJA(MX,MY,MZ)
  real(8) P(0:MX+1,0:MY+1,0:MZ+1), U(0:MX+1,0:MY+1,0:MZ+1), V(0:MX+1,0:MY+1,0:MZ+1), W(0:MX+1,0:MY+1,0:MZ+1)
  real(8) UU(0:MX+1,0:MY+1,0:MZ+1), VV(0:MX+1,0:MY+1,0:MZ+1), WW(0:MX+1,0:MY+1,0:MZ+1)
  real(8) UD(0:MX+1,0:MY+1,0:MZ+1), VD(0:MX+1,0:MY+1,0:MZ+1), WD(0:MX+1,0:MY+1,0:MZ+1)
  real(8) CU(0:MX+1,0:MY+1,0:MZ+1), CV(0:MX+1,0:MY+1,0:MZ+1), CW(0:MX+1,0:MY+1,0:MZ+1)
  real(8) GX(MX,MY,MZ), EY(MX,MY,MZ), TX(MX,MY,MZ), TY(MX,MY,MZ), TZ(MX,MY,MZ)
  real(8) D1(MX,MY,MZ), D2(MX,MY,MZ), D3(MX,MY,MZ), D5(MX,MY,MZ), D6(MX,MY,MZ)
!
!....<<NEW PHYSICAL VELOCITY>>
  do k = 2, NZ-1
    do j = 2, NY-1
      do i = 2, NX-1
        PG = .5D0 * (P(i+1,j,k) - P(i-1,j,k))
        PE = .5D0 * (P(i,j+1,k) - P(i,j-1,k))
        PT = .5D0 * (P(i,j,k+1) - P(i,j,k-1))
        U(i,j,k) = U(i,j,k) + DT * (-PG * GX(i,j,k) - PT * TX(i,j,k))
        V(i,j,k) = V(i,j,k) + DT * (-PE * EY(i,j,k) - PT * TY(i,j,k))
        W(i,j,k) = W(i,j,k) + DT * (-PT * TZ(i,j,k))
      enddo
    enddo
  enddo
!
!....<<NEW CONTRAVARIANT VELOCITY>>
  do k = 2, NZ-1
    do j = 2, NY-1
      do i = 1, NX-1
        PX = .5D0 * (D1(i+1,j,k) + D1(i,j,k)) * (P(i+1,j,k) - P(i,j,k)) &
           + .5D0 * (D6(i+1,j,k) + D6(i,j,k)) * (.25D0 * (P(i+1,j,k+1) + P(i,j,k+1) - P(i+1,j,k-1) - P(i,j,k-1)))
        UU(i,j,k) = UU(i,j,k) - DT * PX
      enddo
    enddo
  enddo
!
  do k = 2 ,NZ-1
    do j = 1, NY-1
      do i = 2, NX-1
        PY= .5D0 * (D2(i,j+1,k) + D2(i,j,k)) * (P(i,j+1,k) - P(i,j,k)) &
           +.5D0 * (D5(i,j+1,k) + D5(i,j,k)) * (.25D0 * (P(i,j+1,k+1) + P(i,j,k+1) - P(i,j,k-1) - P(i,j+1,k-1)))
        VV(i,j,k) = VV(i,j,k) - DT * PY
      enddo
    enddo
  enddo
!
  do k = 1, NZ-1
    do j = 2, NY-1
      do i = 2, NX-1
        PZ = .5D0 * (D6(i,j,k+1) + D6(i,j,k)) * (.25D0 * (P(i+1,j,k) + P(i+1,j,k+1) - P(i-1,j,k+1) - P(i-1,j,k))) &
           + .5D0 * (D5(i,j,k+1) + D5(i,j,k)) * (.25D0 * (P(i,j+1,k) + P(i,j+1,k+1) - P(i,j-1,k+1) - P(i,j-1,k))) &
           + .5D0 * (D3(i,j,k+1) + D3(i,j,k)) * (P(i,j,k+1) - P(i,j,k))
        WW(i,j,k) = WW(i,j,k) - DT * PZ
       enddo
     enddo
   enddo
!
!....<<MODIFY PHYSICAL VELOCITY B.C.>>
  do k = 1, NZ
    do j = 1, NY
      U(0,j,k) = 1.D0
      V(0,j,k) = 0.D0
      W(0,j,k) = 0.D0
      U(1,j,k) = 1.D0
      V(1,j,k) = 0.D0
      W(1,j,k) = 0.D0
    enddo
  enddo
!
  do j = 0, NY+1
    do i = 0, NX+1
      U(i,j,1) = 0.D0
      V(i,j,1) = 0.D0
      W(i,j,1) = 0.D0
      U(i,j,0) = -U(i,j,2)
      V(i,j,0) = -V(i,j,2)
      W(i,j,0) = -W(i,j,2)
      U(i,j,NZ) = U(i,j,NZ-1)
      V(i,j,NZ) = V(i,j,NZ-1)
      W(i,j,NZ) = 0.D0
      U(i,j,NZ+1) = 2.D0 * U(i,j,NZ) - U(i,j,NZ-1)
      V(i,j,NZ+1) = 2.D0 * V(i,j,NZ) - V(i,j,NZ-1)
      W(i,j,NZ+1) = -W(i,j,NZ-1)
    enddo
  enddo
!
  do k = 0, NZ+1
    do i = 0, NX+1
      U(i,1,k) = U(i,2,k)
      V(i,1,k) = 0.D0
      W(i,1,k) = W(i,2,k)
      U(i,0,k) = 2.D0*U(i,1,k) - U(i,2,k)
      V(i,0,k) = -V(i,2,k)
      W(i,0,k) =2.D0 * W(i,1,k) - W(i,2,k)
      U(i,NY,k) = U(i,NY-1,k)
      V(i,NY,k) = 0.D0
      W(i,NY,k) = W(i,NY-1,k)
      U(i,NY+1,k) = 2.D0 * U(i,NY,k) - U(i,NY-1,k)
      V(i,NY+1,k) = -V(i,NY-1,k)
      W(i,NY+1,k) = 2.D0 * W(i,NY,k) - W(i,NY-1,k)
    enddo
  enddo
!
  do k = 2, NZ-1
    do j = 2, NY-1
      DXI = 1.D0 / (X(NX,1,1) - X(NX-1,1,1))
      DUDX = .5D0 * (3.D0 * UD(NX,j,k) - 4.D0 * UD(NX-1,j,k) + UD(NX-2,j,k)) * DXI
      DVDX = .5D0 * (3.D0 * VD(NX,j,k) - 4.D0 * VD(NX-1,j,k) + VD(NX-2,j,k)) * DXI
      DWDX = .5D0 * (3.D0 * WD(NX,j,k) - 4.D0 * WD(NX-1,j,k) + WD(NX-2,j,k)) * DXI
      U(NX,j,k) = UD(NX,j,k) - DT * DUDX
      V(NX,j,k) = VD(NX,j,k) - DT * DVDX
      W(NX,j,k) = WD(NX,j,k) - DT * DWDX
    enddo
  enddo
!
  do k = 0, NZ+1
    do j = 0, NY+1
      U(NX+1,j,k) = 2.D0 * U(NX,j,k) - U(NX-1,j,k)
      V(NX+1,j,k) = 2.D0 * V(NX,j,k) - V(NX-1,j,k)
      W(NX+1,j,k) = 2.D0 * W(NX,j,k) - W(NX-1,j,k)
    enddo
  enddo
!
  do k = KA, KB
    do j = JA, JB
      do i = IA, IB
        U(i,j,k) = 0.D0
        V(i,j,k) = 0.D0
        W(i,j,k) = 0.D0
      enddo
    enddo
  enddo
!
  return
end subroutine calc_new_velocity
