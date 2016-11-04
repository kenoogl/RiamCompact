subroutine eddy_viscosity(MX, MY, MZ, NX, NY, NZ, RE, REI, CS2, Z, U, V, W, GX, EY, TX, TY, TZ, YJA, SGS)
!
!  use params
  implicit none
!
  integer i, j, k
  integer MX, MY, MZ
  integer NX, NY, NZ
  real(8) RE, REI, ZP, CS2
  real(8) FILG, DCON
  real(8) DUDX, DUDY, DUDZ
  real(8) DVDX, DVDY, DVDZ
  real(8) DWDX, DWDY, DWDZ

  real(8) Z(0:MX+1,0:MY+1,0:MZ+1)
  real(8) U(0:MX+1,0:MY+1,0:MZ+1), V(0:MX+1,0:MY+1,0:MZ+1), W(0:MX+1,0:MY+1,0:MZ+1)
  real(8) GX(MX,MY,MZ), EY(MX,MY,MZ), TX(MX,MY,MZ), TY(MX,MY,MZ), TZ(MX,MY,MZ)
  real(8) VAN(MX,MY,MZ), SGS(MX,MY,MZ), YJA(MX,MY,MZ)

!
!=========================================
!....<<SGS EDDY VISCOSITY COEFFICIENT>>
!=========================================
!
!....<<WALL DAMPING FUNCTION>>
    do k = 2, NZ-1
      do j = 2, NY-1
        do i = 2, NX-1
          DUDZ = 1.D0 / (Z(i,j,2) - Z(i,j,1)) * dabs(U(i,j,2))
          ZP = (Z(i,j,k) - Z(i,j,1)) * RE * dsqrt(REI*DUDZ)
          VAN(i,j,k) = 1.D0 - dexp(-ZP/25.D0)
        enddo
      enddo
    enddo
!
!....<<INNER POINTS>>
    do k = 2, NZ-1
      do j = 2, NY-1
        do i = 2, NX-1
!....<<SPACING DERIVATIVES>>
          DUDX = .5D0 * (U(i+1,j,k) - U(i-1,j,k)) * GX(i,j,k) + .5D0 * (U(i,j,k+1) - U(i,j,k-1)) * TX(i,j,k)
          DVDX = .5D0 * (V(i+1,j,k) - V(i-1,j,k)) * GX(i,j,k) + .5D0 * (V(i,j,k+1) - V(i,j,k-1)) * TX(i,j,k)
          DWDX = .5D0 * (W(i+1,j,k) - W(i-1,j,k)) * GX(i,j,k) + .5D0 * (W(i,j,k+1) - W(i,j,k-1)) * TX(i,j,k)
          DUDY = .5D0 * (U(i,j+1,k) - U(i,j-1,k)) * EY(i,j,k) + .5D0 * (U(i,j,k+1) - U(i,j,k-1)) * TY(i,j,k)
          DVDY = .5D0 * (V(i,j+1,k) - V(i,j-1,k)) * EY(i,j,k) + .5D0 * (V(i,j,k+1) - V(i,j,k-1)) * TY(i,j,k)
          DWDY = .5D0 * (W(i,j+1,k) - W(i,j-1,k)) * EY(i,j,k) + .5D0 * (W(i,j,k+1) - W(i,j,k-1)) * TY(i,j,k)
          DUDZ = .5D0 * (U(i,j,k+1) - U(i,j,k-1)) * TZ(i,j,k)
          DVDZ = .5D0 * (V(i,j,k+1) - V(i,j,k-1)) * TZ(i,j,k)
          DWDZ = .5D0 * (W(i,j,k+1) - W(i,j,k-1)) * TZ(i,j,k)
!....<<FILTER LENGTH IN SMAGORINSKY MODEL>>
          FILG = YJA(i,j,k) ** (2.D0/3.D0)
!....<<SGS EDDY VISCOSITY COEFFICIENT>>
          DCON = dsqrt(2.D0 * (DUDX ** 2 + DVDY ** 2 + DWDZ ** 2) &
                 + (DWDY + DVDZ) ** 2 + (DUDZ + DWDX) ** 2 + (DUDY + DVDX) ** 2 )
          SGS(i,j,k) = (CS2 * VAN(i,j,k) ** 2) * FILG * DCON
        enddo
      enddo
    enddo
!
!....<<MODIFY B.C.>>
    do j = 1, NY
      do i = 1, NX
        SGS(i,j,1) = 0.D0
        SGS(i,j,NZ) = SGS(i,j,NZ-1)
      enddo
    enddo
!
    do k = 1, NZ
      do i = 1, NX
        SGS(i,1,k) = SGS(i,2,k)
        SGS(i,NY,k) = SGS(i,NY-1,k)
      enddo
    enddo
!
    do k = 1, NZ
      do j = 1, NY
        SGS(1,j,k) = 0.D0
        SGS(NX,j,k) = SGS(NX-1,j,k)
      enddo
    enddo
end subroutine eddy_viscosity
