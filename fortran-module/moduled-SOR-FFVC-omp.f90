subroutine sor(MX, MY, MZ, NX, NY, NZ, ILAP, ITER, NOXYZ, OMEGA, EPS, RMSP, P, &
               C1, C2, C3, C5, C6, C7, C8, C9, RHS, ERRP, nID, head)
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
  integer ILAP, ITER, NOXYZ

  real(8) ERR, ERR1, ERR2, ERR3, ERR4, RP, OMEGA, EPS, RMSP
  real(8) P(0:MX+1,0:MY+1,0:MZ+1), RHS(MX,MY,MZ)
  real(8) C1(MX,MY,MZ), C2(MX,MY,MZ), C3(MX,MY,MZ), C5(MX,MY,MZ), C6(MX,MY,MZ)
  real(8) C7(MX,MY,MZ), C8(MX,MY,MZ), C9(MX,MY,MZ)
  real(8) ERRP(0:MX+1,0:MY+1,0:MZ+1)
  real(8) ERR_BUF

  integer nID(0:5), head(0:2)

  integer ist, jst, kst
  integer iend, jend, kend
  integer ihead, jhead, khead

  integer istp, jstp, kstp

  integer pg, ierr
  pg=0

! 計算開始位置の設定
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

! エラー値計算開始位置の設定
  istp = 1
  jstp = 1
  kstp = 1
  if( nID(X_MINUS) < 0 ) istp = 2
  if( nID(Y_MINUS) < 0 ) jstp = 2
  if( nID(Z_MINUS) < 0 ) kstp = 2

! グローバルの開始位置が奇数か偶数かを設定
  ihead = mod(ist+head(X_DIR), 2) ! 0:偶数、1:奇数
  jhead = mod(jst+head(Y_DIR), 2) ! 0:偶数、1:奇数
  khead = mod(kst+head(Z_DIR), 2) ! 0:偶数、1:奇数

!....<<SOR METHOD>>
  do ILAP = 1, ITER
   ERR1 = 0.D0
   ERR2 = 0.D0
   ERR3 = 0.D0
   ERR4 = 0.D0

   !do k = 2, NZ-1, 2
    do k = kst+khead, kend, 2
!$OMP PARALLEL
!$OMP DO &
!$OMP REDUCTION(+:ERR1) &
!$OMP PRIVATE(RP) 
   !do j = 2, NY-1
   !do i = 2 + mod(j , 2) , NX-1, 2
    do j = jst, jend
    do i = ist + mod(j+head(Y_DIR)+ihead , 2) , iend, 2
       RP = .5D0 * ((C1(i,j,k) * (P(i+1,j,k)   + P(i-1,j,k)) &
          +   C2(i,j,k) * (P(i,j+1,k)   + P(i,j-1,k)) &
          +   C3(i,j,k) * (P(i,j,k+1)   + P(i,j,k-1)) &
          +   C5(i,j,k) * (P(i,j+1,k+1) - P(i,j+1,k-1) - P(i,j-1,k+1) + P(i,j-1,k-1)) * .25D0 &
          +   C6(i,j,k) * (P(i+1,j,k+1) - P(i+1,j,k-1) - P(i-1,j,k+1) + P(i-1,j,k-1)) * .25D0 &
          +   C7(i,j,k) * (P(i+1,j,k)   - P(i-1,j,k))*.5D0 &
          +   C8(i,j,k) * (P(i,j+1,k)   - P(i,j-1,k))*.5D0 &
          +   C9(i,j,k) * (P(i,j,k+1)   - P(i,j,k-1))*.5D0) &
          + RHS(i,j,k)) / ((C1(i,j,k) + C2(i,j,k) + C3(i,j,k))) - P(i,j,k)
       ERR1 = ERR1 + RP * RP
       ERRP(i,j,k) = RP * RP
       P(i,j,k) = P(i,j,k) + OMEGA * RP
    enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL
    enddo

    call cpm_BndCommS3D(P,MX,MY,MZ,1,1,CPM_DOUBLE,pg,ierr)

    do k = kst+khead, kend, 2
!
!$OMP PARALLEL
!$OMP DO &
!$OMP REDUCTION(+:ERR2) &
!$OMP PRIVATE(RP) 
   !do j = 2, NY-1
   !do i = 2 + mod(j+1 , 2) , NX-1, 2
    do j = jst, jend
    do i = ist + mod(j+head(Y_DIR)+ihead+1 , 2) , iend, 2
       RP = .5D0 * ((C1(i,j,k) * (P(i+1,j,k)   + P(i-1,j,k)) &
          +   C2(i,j,k) * (P(i,j+1,k)   + P(i,j-1,k)) &
          +   C3(i,j,k) * (P(i,j,k+1)   + P(i,j,k-1)) &
          +   C5(i,j,k) * (P(i,j+1,k+1) - P(i,j+1,k-1) - P(i,j-1,k+1) + P(i,j-1,k-1)) * .25D0 &
          +   C6(i,j,k) * (P(i+1,j,k+1) - P(i+1,j,k-1) - P(i-1,j,k+1) + P(i-1,j,k-1)) * .25D0 &
          +   C7(i,j,k) * (P(i+1,j,k)   - P(i-1,j,k))*.5D0 &
          +   C8(i,j,k) * (P(i,j+1,k)   - P(i,j-1,k))*.5D0 &
          +   C9(i,j,k) * (P(i,j,k+1)   - P(i,j,k-1))*.5D0) &
          + RHS(i,j,k)) / ((C1(i,j,k) + C2(i,j,k) + C3(i,j,k))) - P(i,j,k)
       ERR2 = ERR2 + RP * RP
       ERRP(i,j,k) = RP * RP
       P(i,j,k) = P(i,j,k) + OMEGA * RP
    enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL
!
    enddo

    call cpm_BndCommS3D(P,MX,MY,MZ,1,1,CPM_DOUBLE,pg,ierr)

   !do k = 3, NZ-1, 2
    do k = kst+1-khead, kend, 2
!$OMP PARALLEL
!$OMP DO &
!$OMP REDUCTION(+:ERR3) &
!$OMP PRIVATE(RP) 
   !do j = 2, NY-1
   !do i = 2 + mod(j , 2) , NX-1, 2
    do j = jst, jend
    do i = ist + mod(j+head(Y_DIR)+ihead , 2) , iend, 2
       RP = .5D0 * ((C1(i,j,k) * (P(i+1,j,k)   + P(i-1,j,k)) &
          +   C2(i,j,k) * (P(i,j+1,k)   + P(i,j-1,k)) &
          +   C3(i,j,k) * (P(i,j,k+1)   + P(i,j,k-1)) &
          +   C5(i,j,k) * (P(i,j+1,k+1) - P(i,j+1,k-1) - P(i,j-1,k+1) + P(i,j-1,k-1)) * .25D0 &
          +   C6(i,j,k) * (P(i+1,j,k+1) - P(i+1,j,k-1) - P(i-1,j,k+1) + P(i-1,j,k-1)) * .25D0 &
          +   C7(i,j,k) * (P(i+1,j,k)   - P(i-1,j,k))*.5D0 &
          +   C8(i,j,k) * (P(i,j+1,k)   - P(i,j-1,k))*.5D0 &
          +   C9(i,j,k) * (P(i,j,k+1)   - P(i,j,k-1))*.5D0) &
          + RHS(i,j,k)) / ((C1(i,j,k) + C2(i,j,k) + C3(i,j,k))) - P(i,j,k)
       ERR3 = ERR3 + RP * RP
       ERRP(i,j,k) = RP * RP
       P(i,j,k) = P(i,j,k) + OMEGA * RP
    enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL
    enddo

    call cpm_BndCommS3D(P,MX,MY,MZ,1,1,CPM_DOUBLE,pg,ierr)

    do k = kst+1-khead, kend, 2
!
!$OMP PARALLEL
!$OMP DO &
!$OMP REDUCTION(+:ERR4) &
!$OMP PRIVATE(RP) 
   !do j = 2, NY-1
   !do i = 2 + mod(j+1 , 2) , NX-1, 2
    do j = jst, jend
    do i = ist + mod(j+head(Y_DIR)+ihead+1 , 2) , iend, 2
       RP = .5D0 * ((C1(i,j,k) * (P(i+1,j,k)   + P(i-1,j,k)) &
          +   C2(i,j,k) * (P(i,j+1,k)   + P(i,j-1,k)) &
          +   C3(i,j,k) * (P(i,j,k+1)   + P(i,j,k-1)) &
          +   C5(i,j,k) * (P(i,j+1,k+1) - P(i,j+1,k-1) - P(i,j-1,k+1) + P(i,j-1,k-1)) * .25D0 &
          +   C6(i,j,k) * (P(i+1,j,k+1) - P(i+1,j,k-1) - P(i-1,j,k+1) + P(i-1,j,k-1)) * .25D0 &
          +   C7(i,j,k) * (P(i+1,j,k)   - P(i-1,j,k))*.5D0 &
          +   C8(i,j,k) * (P(i,j+1,k)   - P(i,j-1,k))*.5D0 &
          +   C9(i,j,k) * (P(i,j,k+1)   - P(i,j,k-1))*.5D0) &
          + RHS(i,j,k)) / ((C1(i,j,k) + C2(i,j,k) + C3(i,j,k))) - P(i,j,k)
       ERR4 = ERR4 + RP * RP
       ERRP(i,j,k) = RP * RP
       P(i,j,k) = P(i,j,k) + OMEGA * RP
    enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL
!
    enddo

    call cpm_BndCommS3D(P,MX,MY,MZ,1,1,CPM_DOUBLE,pg,ierr)

!....<<CONVERGE OR NOT ?>>
   !ERR = ERR1 + ERR2 + ERR3 + ERR4
    ERR = 0.0D0
    do k=kstp,kend
    do j=jstp,jend
    do i=istp,iend
      ERR = ERR + ERRP(i,j,k)
    enddo
    enddo
    enddo

    pg = 0
    ERR_BUF = ERR
    call cpm_Allreduce(ERR_BUF,ERR,1,CPM_DOUBLE,CPM_SUM,pg,ierr)

    RMSP = dsqrt( ERR / dfloat(NOXYZ))
    if(RMSP <= EPS) exit
  enddo
  return
!
end subroutine SOR
