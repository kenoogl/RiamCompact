subroutine calc_new_velocity(MX, MY, MZ, NX, NY, NZ, KA, KB, JA, JB, IA, IB, DT, X, P, U, V, W, &
                             UD, VD, WD, UU, VV, WW, GX, EY, TX, TY, TZ, YJA, D1, D2, D3, D5, D6, &
                             nID , head, tail)
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
  integer KA, KB, JA, JB, IA, IB 
!
  real(8) DT
  real(8) PG, PE, PT, PX, PY, PZ, DXI
  real(8) DUDX, DVDX, DWDX
!
  real(8) X(0:MX+1,0:MY+1,0:MZ+1), YJA(0:MX+1,0:MY+1,0:MZ+1)
  real(8) P(0:MX+1,0:MY+1,0:MZ+1)
  real(8) U(0:MX+1,0:MY+1,0:MZ+1), V(0:MX+1,0:MY+1,0:MZ+1), W(0:MX+1,0:MY+1,0:MZ+1)
  real(8) UU(0:MX+1,0:MY+1,0:MZ+1), VV(0:MX+1,0:MY+1,0:MZ+1), WW(0:MX+1,0:MY+1,0:MZ+1)
  real(8) UD(-1:MX+2,-1:MY+2,-1:MZ+2), VD(-1:MX+2,-1:MY+2,-1:MZ+2), WD(-1:MX+2,-1:MY+2,-1:MZ+2)
  real(8) GX(0:MX+1,0:MY+1,0:MZ+1), EY(0:MX+1,0:MY+1,0:MZ+1)
  real(8) TX(0:MX+1,0:MY+1,0:MZ+1), TY(0:MX+1,0:MY+1,0:MZ+1), TZ(0:MX+1,0:MY+1,0:MZ+1)
  real(8) D1(0:MX+1,0:MY+1,0:MZ+1), D2(0:MX+1,0:MY+1,0:MZ+1), D3(0:MX+1,0:MY+1,0:MZ+1)
  real(8) D5(0:MX+1,0:MY+1,0:MZ+1), D6(0:MX+1,0:MY+1,0:MZ+1)
!
  integer nID(0:5), head(0:2), tail(0:2)

  integer ist, jst, kst
  integer iend, jend, kend

  integer l_IA, l_IB, l_JA, l_JB, l_KA, l_KB

  ist = 2
  jst = 2
  kst = 2
  iend = NX-1
  jend = NY-1
  kend = NZ-1
  if( nID(X_MINUS) >= 0 ) ist  = 1  !-X方向に隣接ランクがある場合は開始位置を変更
  if( nID(X_PLUS ) >= 0 ) iend = NX !+X方向に隣接ランクがある場合は終了位置を変更
  if( nID(Y_MINUS) >= 0 ) jst  = 1  !-Y方向に隣接ランクがある場合は開始位置を変更
  if( nID(Y_PLUS ) >= 0 ) jend = NY !+Y方向に隣接ランクがある場合は終了位置を変更
  if( nID(Z_MINUS) >= 0 ) kst  = 1  !-Z方向に隣接ランクがある場合は開始位置を変更
  if( nID(Z_PLUS ) >= 0 ) kend = NZ !+Z方向に隣接ランクがある場合は終了位置を変更  
!
!....<<NEW PHYSICAL VELOCITY>>
!
!$OMP PARALLEL
!$OMP DO &
!$OMP PRIVATE(PG, PE, PT) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = kst, kend 
  do j = jst, jend
  do i = ist, iend
    PG = .5D0 * (P(i+1,j,k) - P(i-1,j,k))
    PE = .5D0 * (P(i,j+1,k) - P(i,j-1,k))
    PT = .5D0 * (P(i,j,k+1) - P(i,j,k-1))
    U(i,j,k) = U(i,j,k) + DT * (-PG * GX(i,j,k) - PT * TX(i,j,k))
    V(i,j,k) = V(i,j,k) + DT * (-PE * EY(i,j,k) - PT * TY(i,j,k))
    W(i,j,k) = W(i,j,k) + DT * (-PT * TZ(i,j,k))
  enddo
  enddo
  enddo
!$OMP END DO
!
!....<<NEW CONTRAVARIANT VELOCITY>>
!
!$OMP DO &
!$OMP PRIVATE(PX) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = kst  , kend
  do j = jst  , jend
  do i = ist-1, iend
    PX = .5D0 * (D1(i+1,j,k) + D1(i,j,k)) * (P(i+1,j,k) - P(i,j,k)) &
       + .5D0 * (D6(i+1,j,k) + D6(i,j,k)) &
       * (.25D0 * (P(i+1,j,k+1) + P(i,j,k+1) - P(i+1,j,k-1) - P(i,j,k-1)))
    UU(i,j,k) = UU(i,j,k) - DT * PX
  enddo
  enddo
  enddo
!$OMP END DO
!
!$OMP DO &
!$OMP PRIVATE(PY) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = kst  , kend
  do j = jst-1, jend
  do i = ist  , iend
    PY= .5D0 * (D2(i,j+1,k) + D2(i,j,k)) * (P(i,j+1,k) - P(i,j,k)) &
      +.5D0 * (D5(i,j+1,k) + D5(i,j,k))  &
      * (.25D0 * (P(i,j+1,k+1) + P(i,j,k+1) - P(i,j,k-1) - P(i,j+1,k-1)))
    VV(i,j,k) = VV(i,j,k) - DT * PY
  enddo
  enddo
  enddo
!$OMP END DO
!
!$OMP DO &
!$OMP PRIVATE(PZ) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = kst-1, kend
  do j = jst  , jend
  do i = ist  , iend
    PZ = .5D0 * (D6(i,j,k+1) + D6(i,j,k)) &
       * (.25D0 * (P(i+1,j,k) + P(i+1,j,k+1) - P(i-1,j,k+1) - P(i-1,j,k))) &
       + .5D0 * (D5(i,j,k+1) + D5(i,j,k)) &
       * (.25D0 * (P(i,j+1,k) + P(i,j+1,k+1) - P(i,j-1,k+1) - P(i,j-1,k))) &
       + .5D0 * (D3(i,j,k+1) + D3(i,j,k)) * (P(i,j,k+1) - P(i,j,k))
    WW(i,j,k) = WW(i,j,k) - DT * PZ
   enddo
   enddo
   enddo
!$OMP END DO
!
!....<<MODIFY PHYSICAL VELOCITY B.C.>>
!
! -X方向に隣接ランクがない場合境界条件を与える
  if( nID(X_MINUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NY, NZ)
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
!$OMP END DO
  endif
!
! -Z方向に隣接ランクがない場合境界条件を与える
  if( nID(Z_MINUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ)
    do j = 0, NY+1
    do i = 0, NX+1
      U(i,j,1) = 0.D0
      V(i,j,1) = 0.D0
      W(i,j,1) = 0.D0
      U(i,j,0) = -U(i,j,2)
      V(i,j,0) = -V(i,j,2)
      W(i,j,0) = -W(i,j,2)
    enddo
    enddo
!$OMP END DO
  endif
!
! +Z方向に隣接ランクがない場合境界条件を与える
  if( nID(Z_PLUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ)
    do j = 0, NY+1
    do i = 0, NX+1
      U(i,j,NZ) = U(i,j,NZ-1)
      V(i,j,NZ) = V(i,j,NZ-1)
      W(i,j,NZ) = 0.D0
      U(i,j,NZ+1) = 2.D0 * U(i,j,NZ) - U(i,j,NZ-1)
      V(i,j,NZ+1) = 2.D0 * V(i,j,NZ) - V(i,j,NZ-1)
      W(i,j,NZ+1) = -W(i,j,NZ-1)
    enddo
    enddo
!$OMP END DO
  endif

! -Y方向に隣接ランクがない場合境界条件を与える
  if( nID(Y_MINUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ)
    do k = 0, NZ+1
    do i = 0, NX+1
      U(i,1,k) = U(i,2,k)
      V(i,1,k) = 0.D0
      W(i,1,k) = W(i,2,k)
      U(i,0,k) = 2.D0*U(i,1,k) - U(i,2,k)
      V(i,0,k) = -V(i,2,k)
      W(i,0,k) =2.D0 * W(i,1,k) - W(i,2,k)
    enddo
    enddo
!$OMP END DO
  endif
!
! +Y方向に隣接ランクがない場合境界条件を与える
  if( nID(Y_PLUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ)
    do k = 0, NZ+1
    do i = 0, NX+1
      U(i,NY,k) = U(i,NY-1,k)
      V(i,NY,k) = 0.D0
      W(i,NY,k) = W(i,NY-1,k)
      U(i,NY+1,k) = 2.D0 * U(i,NY,k) - U(i,NY-1,k)
      V(i,NY+1,k) = -V(i,NY-1,k)
      W(i,NY+1,k) = 2.D0 * W(i,NY,k) - W(i,NY-1,k)
    enddo
    enddo
!$OMP END DO
  endif

!
! +X方向に隣接ランクがない場合境界条件を与える
  if( nID(X_PLUS) < 0 ) then
!$OMP DO &
!$OMP PRIVATE(DXI, DUDX, DVDX, DWDX)
   !do k = 2, NZ-1
   !do j = 2, NY-1
    do k = kst, kend
    do j = jst, jend
      DXI = 1.D0 / (X(NX,1,1) - X(NX-1,1,1))
      DUDX = .5D0 * (3.D0 * UD(NX,j,k) - 4.D0 * UD(NX-1,j,k) + UD(NX-2,j,k)) * DXI
      DVDX = .5D0 * (3.D0 * VD(NX,j,k) - 4.D0 * VD(NX-1,j,k) + VD(NX-2,j,k)) * DXI
      DWDX = .5D0 * (3.D0 * WD(NX,j,k) - 4.D0 * WD(NX-1,j,k) + WD(NX-2,j,k)) * DXI
      U(NX,j,k) = UD(NX,j,k) - DT * DUDX
      V(NX,j,k) = VD(NX,j,k) - DT * DVDX
      W(NX,j,k) = WD(NX,j,k) - DT * DWDX
    enddo
    enddo
!$OMP END DO
  endif

!
! +X方向に隣接ランクがない場合境界条件を与える
  if( nID(X_PLUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ)
    do k = 0, NZ+1
    do j = 0, NY+1
      U(NX+1,j,k) = 2.D0 * U(NX,j,k) - U(NX-1,j,k)
      V(NX+1,j,k) = 2.D0 * V(NX,j,k) - V(NX-1,j,k)
      W(NX+1,j,k) = 2.D0 * W(NX,j,k) - W(NX-1,j,k)
    enddo
    enddo
!$OMP END DO
  endif

! 立方体が含まれないランクのときは ゼロクリアしないで return する

  l_IA = max(1, IA-head(X_DIR))
  l_IB = min(NX,IB-head(X_DIR))
  l_JA = max(1, JA-head(Y_DIR))
  l_JB = min(NY,JB-head(Y_DIR))
  l_KA = max(1, KA-head(Z_DIR))
  l_KB = min(NY,KB-head(Z_DIR))
!
!$OMP DO 
  do k = l_KA, l_KB
  do j = l_JA, l_JB
  do i = l_IA, l_IB
    U(i,j,k) = 0.D0
    V(i,j,k) = 0.D0
    W(i,j,k) = 0.D0
  enddo
  enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL
!
  return
end subroutine calc_new_velocity
