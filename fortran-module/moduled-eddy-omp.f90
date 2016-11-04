subroutine eddy_viscosity(MX, MY, MZ, NX, NY, NZ, RE, REI, CS2, Z, U, V, W, &
                          GX, EY, TX, TY, TZ, YJA, SGS, VAN, &
                          nID )
!
!$ use omp_lib
!  use params
  implicit none
!
  include 'cpm_fparam.fi'
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
  real(8) GX(0:MX+1,0:MY+1,0:MZ+1), EY(0:MX+1,0:MY+1,0:MZ+1)
  real(8) TX(0:MX+1,0:MY+1,0:MZ+1), TY(0:MX+1,0:MY+1,0:MZ+1), TZ(0:MX+1,0:MY+1,0:MZ+1)
  real(8) VAN(MX,MY,MZ), SGS(0:MX+1,0:MY+1,0:MZ+1), YJA(0:MX+1,0:MY+1,0:MZ+1)
!
  integer nID(0:5)
  integer ist,jst,kst
  integer iend,jend,kend
!
!  計算開始位置,終了位置の設定
  ist = 2
  jst = 2
  kst = 2
  iend = NX-1
  jend = NY-1
  kend = NZ-1
  if( nID(X_MINUS) >= 0 ) ist = 1;  ! -X方向に隣接ランクがある場合は開始位置を変更
  if( nID(X_PLUS ) >= 0 ) iend= NX; ! +X方向に隣接ランクがある場合は終了位置を変更
  if( nID(Y_MINUS) >= 0 ) jst = 1;  ! -Y方向に隣接ランクがある場合は開始位置を変更
  if( nID(Y_PLUS ) >= 0 ) jend= NY; ! +Y方向に隣接ランクがある場合は終了位置を変更
  if( nID(Z_MINUS) >= 0 ) kst = 1;  ! -Z方向に隣接ランクがある場合は開始位置を変更
  if( nID(Z_PLUS ) >= 0 ) kend= NZ; ! +Z方向に隣接ランクがある場合は終了位置を変更

!
!=========================================
!....<<SGS EDDY VISCOSITY COEFFICIENT>>
!=========================================
!
!....<<WALL DAMPING FUNCTION>>

!$OMP PARALLEL
!$OMP DO &
!$OMP PRIVATE(DUDZ, ZP) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = kst, kend
  do j = jst, jend
  do i = ist, iend
    DUDZ = 1.D0 / (Z(i,j,2) - Z(i,j,1)) * dabs(U(i,j,2))
    ZP = (Z(i,j,k) - Z(i,j,1)) * RE * dsqrt(REI*DUDZ)
    VAN(i,j,k) = 1.D0 - dexp(-ZP/25.D0)
  enddo
  enddo
  enddo
!$OMP END DO

!
!....<<INNER POINTS>>

!$OMP DO &
!$OMP PRIVATE(DUDX, DVDX, DWDX, DUDY, DVDY, DWDY) &
!$OMP PRIVATE(DUDZ, DVDZ, DWDZ, FILG, DCON) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = kst, kend
  do j = jst, jend
  do i = ist, iend
!....<<SPACING DERIVATIVES>>
    DUDX = .5D0 * (U(i+1,j,k) - U(i-1,j,k)) * GX(i,j,k) &
         + .5D0 * (U(i,j,k+1) - U(i,j,k-1)) * TX(i,j,k)
    DVDX = .5D0 * (V(i+1,j,k) - V(i-1,j,k)) * GX(i,j,k) &
         + .5D0 * (V(i,j,k+1) - V(i,j,k-1)) * TX(i,j,k)
    DWDX = .5D0 * (W(i+1,j,k) - W(i-1,j,k)) * GX(i,j,k) &
         + .5D0 * (W(i,j,k+1) - W(i,j,k-1)) * TX(i,j,k)
    DUDY = .5D0 * (U(i,j+1,k) - U(i,j-1,k)) * EY(i,j,k) &
         + .5D0 * (U(i,j,k+1) - U(i,j,k-1)) * TY(i,j,k)
    DVDY = .5D0 * (V(i,j+1,k) - V(i,j-1,k)) * EY(i,j,k) &
         + .5D0 * (V(i,j,k+1) - V(i,j,k-1)) * TY(i,j,k)
    DWDY = .5D0 * (W(i,j+1,k) - W(i,j-1,k)) * EY(i,j,k) &
         + .5D0 * (W(i,j,k+1) - W(i,j,k-1)) * TY(i,j,k)
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
!$OMP END DO
!
!....<<MODIFY B.C.>>
! -Z方向に隣接ランクがない場合境界条件を与える
  if( nID(Z_MINUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ)
    do j = 1, NY
    do i = 1, NX
      SGS(i,j,1) = 0.D0
    enddo
    enddo
!$OMP END DO
  endif
!
! +Z方向に隣接ランクがない場合境界条件を与える
  if( nID(Z_PLUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ)
    do j = 1, NY
    do i = 1, NX
      SGS(i,j,NZ) = SGS(i,j,NZ-1)
    enddo
    enddo
!$OMP END DO
  endif
!
! -Y方向に隣接ランクがない場合境界条件を与える
  if( nID(Y_MINUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ)
    do k = 1, NZ
    do i = 1, NX
      SGS(i,1,k) = SGS(i,2,k)
    enddo
    enddo
!$OMP END DO
  endif
! +Y方向に隣接ランクがない場合境界条件を与える
  if( nID(Y_PLUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ)
    do k = 1, NZ
    do i = 1, NX
      SGS(i,NY,k) = SGS(i,NY-1,k)
    enddo
    enddo
!$OMP END DO
  endif

! -X方向に隣接ランクがない場合境界条件を与える
  if( nID(X_MINUS) < 0 ) then 
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ)
    do k = 1, NZ
    do j = 1, NY
      SGS(1,j,k) = 0.D0
    enddo
    enddo
!$OMP END DO
  endif
! +X方向に隣接ランクがない場合境界条件を与える
  if( nID(X_PLUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ)
    do k = 1, NZ
    do j = 1, NY
      SGS(NX,j,k) = SGS(NX-1,j,k)
    enddo
    enddo
!$OMP END DO
  endif
!$OMP END PARALLEL

end subroutine eddy_viscosity
