subroutine calc_contravariant_v(MX, MY, MZ, NX, NY, NZ, DT, U, V, W, GX, EY, TX, TY, TZ, &
                                UU, VV, WW, YJA, RHS, &
                                CU, CV, CW, nID)
!
!  use params
  implicit none
!
  include 'cpm_fparam.fi'
!
  integer i, j, k
  integer MX, MY, MZ
  integer NX, NY, NZ
!
  real(8) DT
  real(8) U(0:MX+1,0:MY+1,0:MZ+1), V(0:MX+1,0:MY+1,0:MZ+1), W(0:MX+1,0:MY+1,0:MZ+1)
  real(8) UU(0:MX+1,0:MY+1,0:MZ+1), VV(0:MX+1,0:MY+1,0:MZ+1), WW(0:MX+1,0:MY+1,0:MZ+1)
  real(8) CU(0:MX+1,0:MY+1,0:MZ+1), CV(0:MX+1,0:MY+1,0:MZ+1), CW(0:MX+1,0:MY+1,0:MZ+1)
  real(8) GX(0:MX+1,0:MY+1,0:MZ+1), EY(0:MX+1,0:MY+1,0:MZ+1)
  real(8) TX(0:MX+1,0:MY+1,0:MZ+1), TY(0:MX+1,0:MY+1,0:MZ+1), TZ(0:MX+1,0:MY+1,0:MZ+1)
  real(8) YJA(0:MX+1,0:MY+1,0:MZ+1), RHS(MX,MY,MZ)
!
  integer nID(0:5)
  integer ist, jst, kst
  integer iend, jend, kend

!
!
!==========================================================
!....<<NEW PRESSURE, PHYSICAL & CONTRAVARIANT VELOCITY>>
!==========================================================
!....<<CONTRAVARIANT VELOCITY OF A CELL CENTER>>
!

! 計算開始、終了位置の設定
  ist = 2
  jst = 2
  kst = 2
  iend = NX - 1
  jend = NY - 1
  kend = NZ - 1
  if( nID(X_MINUS ) >= 0 ) ist  = 1;  !-X方向に隣接ランクがある場合は開始位置を変更
  if( nID(X_PLUS  ) >= 0 ) iend = NX; !+X方向に隣接ランクがある場合は終了位置を変更
  if( nID(Y_MINUS ) >= 0 ) jst  = 1;  !-Y方向に隣接ランクがある場合は開始位置を変更
  if( nID(Y_PLUS  ) >= 0 ) jend = NY; !+Y方向に隣接ランクがある場合は終了位置を変更
  if( nID(Z_MINUS ) >= 0 ) kst  = 1;  !-Z方向に隣接ランクがある場合は開始位置を変更
  if( nID(Z_PLUS  ) >= 0 ) kend = NZ; !+Z方向に隣接ランクがある場合は終了位置を変更

!$OMP PARALLEL
!$OMP DO &
!$OMP FIRSTPRIVATE(ist, jst, kst, iend, jend, kend) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = kst-1, kend+1
  do j = jst-1, jend+1
  do i = ist-1, iend+1
    CU(i,j,k) = YJA(i,j,k) * U(i,j,k) * GX(i,j,k)
    CV(i,j,k) = YJA(i,j,k) * V (i,j,k) * EY(i,j,k)
    CW(i,j,k) = YJA(i,j,k) * (U(i,j,k) * TX(i,j,k) + V(i,j,k) * TY(i,j,k) + W(i,j,k) * TZ(i,j,k))
  enddo
  enddo
  enddo
!$OMP END DO
!
!....<<CONTRAVARIANT VELOCITY OF A CELL SURFACE>>
!$OMP DO &
!$OMP FIRSTPRIVATE(ist, jst, kst, iend, jend, kend) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = kst  , kend
  do j = jst  , jend
  do i = ist-1, iend
    UU(i,j,k) = .5D0 * (CU(i+1,j,k) + CU(i,j,k))
  enddo
  enddo
  enddo
!$OMP END DO
!
!$OMP DO &
!$OMP FIRSTPRIVATE(ist, jst, kst, iend, jend, kend) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k= kst  , kend
  do j= jst-1, jend
  do i= ist  , iend
    VV(i,j,k) = .5D0 * ( CV(i,j+1, k) + CV(i,j,k))
  enddo
  enddo
  enddo
!$OMP END DO
!
!$OMP DO &
!$OMP FIRSTPRIVATE(ist, jst, kst, iend, jend, kend) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = kst-1, kend
  do j = jst  , jend
  do i = ist  , iend
    WW(i,j,k) = .5D0 * (CW(i,j,k+1) + CW(i,j,k))
  enddo
  enddo
  enddo
!$OMP END DO
!
!....<<RIGHT HAND SIDE>>
!$OMP DO &
!$OMP FIRSTPRIVATE(ist, jst, kst, iend, jend, kend, DT) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = kst, kend
  do j = jst, jend
  do i = ist, iend
     RHS(i,j,k) = -((UU(i,j,k) - UU(i-1,j,k) &
                + VV(i,j,k)   - VV(i,j-1,k) &
                + WW(i,j,k)   -  WW(i,j,k-1)) / YJA(i,j,k)) / DT
  enddo
  enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL
!
  return
end subroutine calc_contravariant_v
