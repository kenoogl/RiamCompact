subroutine calc_new_pressure(MX, MY, MZ, NX, NY, NZ, REI, P, Z, W, C1, C2, C3, &
                             C5, C6, C7, C8, C9, nID, GNX, GNY, GNZ, PP)
!
!  use params
  implicit none
!
  include 'cpm_fparam.fi'
!
  integer i, j, k
  integer MX, MY, MZ
  integer NX, NY, NZ
  integer GNX, GNY, GNZ

  real(8) REI, PP
  real(8) P(0:MX+1,0:MY+1,0:MZ+1),Z(0:MX+1,0:MY+1,0:MZ+1), W(0:MX+1,0:MY+1,0:MZ+1)
  real(8) C1(MX,MY,MZ), C2(MX,MY,MZ), C3(MX,MY,MZ), C5(MX,MY,MZ), C6(MX,MY,MZ)
  real(8) C7(MX,MY,MZ), C8(MX,MY,MZ), C9(MX,MY,MZ)

  integer nID(0:5)
  integer ist,jst,kst
  integer iend, jend, kend

  integer ierr, pg
  real(8) PP_BUF

  pg=0

  iend = NX-1
  jend = NY-1
  kend = NZ-1
  if( nID(X_PLUS) < 0 ) iend = NX
  if( nID(Y_PLUS) < 0 ) jend = NY
  if( nID(Z_PLUS) < 0 ) kend = NZ

!....<<AVERAGING & SUBTRACTION>>
  PP = 0.D0
!
!$OMP PARALLEL
!$OMP DO &
!$OMP REDUCTION(+:PP) 
  do k = 1, kend 
  do j = 1, jend
  do i = 1, iend
    PP = PP + P(i,j,k)
  enddo
  enddo
  enddo
!$OMP END DO


!$OMP MASTER
  pg = 0
  PP_BUF = PP
  call cpm_Allreduce(PP_BUF,PP,1,CPM_DOUBLE,CPM_SUM,pg,ierr)
!$OMP END MASTER
!$OMP BARRIER

!
!$OMP DO &
!$OMP PRIVATE(PP) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = 1, NZ
    do j = 1, NY
      do i = 1, NX
        P(i,j,k) = P(i,j,k) - PP / dfloat(GNX*GNY*GNZ)
      enddo
    enddo
  enddo
!$OMP END DO
!
!....<<MODIFY PRESSURE B.C.>>

! 計算位置の設定
  ist  = 2
  jst  = 2
  kst  = 2
  iend = NX - 1
  jend = NY - 1
  kend = NZ - 1
  if( nID(X_MINUS) >= 0 ) ist = 1  ! -X方向に隣接ランクがある場合は開始位置を変更
  if( nID(X_PLUS ) >= 0 ) iend= NX ! +X方向に隣接ランクがある場合は終了位置を変更
  if( nID(Y_MINUS) >= 0 ) jst = 1  ! -Y方向に隣接ランクがある場合は開始位置を変更
  if( nID(Y_PLUS ) >= 0 ) jend= NY ! +Y方向に隣接ランクがある場合は終了位置を変更
  if( nID(Z_MINUS) >= 0 ) kst = 1  ! -Z方向に隣接ランクがある場合は開始位置を変更
  if( nID(Z_PLUS ) >= 0 ) kend= NZ ! +Z方向に隣接ランクがある場合は終了位置を変更

  k = 1
! -Z方向に隣接ランクがない場合境界条件を与える
  if( nID(Z_MINUS) < 0 ) then
!$OMP DO
    do j = jst, jend
    do i = ist, iend
      P(i,j,1) = P(i,j,2) - .5D0 * (Z(i,j,2) - Z(i,j,0)) &
               * (C1(i,j,k) * (W(i+1,j,k) - 2.D0 * W(i,j,k) + W(i-1,j,k)) &
               +  C2(i,j,k) * (W(i,j+1,k) - 2.D0 * W(i,j,k) + W(i,j-1,k)) &
               +  C3(i,j,k) * (W(i,j,k+1) - 2.D0 * W(i,j,k) + W(i,j,k-1)) &
               +  C5(i,j,k) * (W(i,j+1,k+1) - W(i,j+1,k-1) - W(i,j-1,k+1) &
               +  W(i,j-1,k-1)) * .25D0 &
               +  C6(i,j,k) * (W(i+1,j,k+1) - W(i+1,j,k-1) - W(i-1,j,k+1) &
               +  W(i-1,j,k-1)) * .25D0 &
               +  C7(i,j,k) * (W(i+1,j,k) - W(i-1,j,k)) * .5D0 &
               +  C8(i,j,k) * (W(i,j+1,k) - W(i,j-1,k)) * .5D0 &
               +  C9(i,j,k) * (W(i,j,k+1) - W(i,j,k-1)) * .5D0) * REI
    enddo
    enddo
!$OMP END DO
  endif
!
! +Z方向に隣接ランクがない場合境界条件を与える
  if( nID(Z_PLUS) < 0 ) then
!$OMP DO
    do j = jst, jend
    do i = ist, iend
      P(i,j,NZ) = P(i,j,NZ-1)
    enddo
    enddo
!$OMP END DO
  endif
!

! -Y方向に隣接ランクがない場合境界条件を与える
  if( nID(Y_MINUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NZ)
    do  k =1, NZ
    do i = 1, NX
      P(i,1,k) = P(i,2,k)
    enddo
    enddo
!$OMP END DO
  endif
!
! +Y方向に隣接ランクがない場合境界条件を与える
  if( nID(Y_PLUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ)
    do  k =1, NZ
    do i = 1, NX
      P(i,NY,k) = P(i,NY-1,k)
    enddo
    enddo
!$OMP END DO
  endif

! -X方向に隣接ランクがない場合境界条件を与える
  if( nID(X_MINUS) < 0 ) then
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY)
    do k = 1, NZ
    do j = 1, NY
      P(1,j,k) = P(2,j,k)
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
      P(NX,j,k) = P(NX-1,j,k)
    enddo
    enddo
!$OMP END DO
  endif
!$OMP END PARALLEL

  return
end subroutine calc_new_pressure
