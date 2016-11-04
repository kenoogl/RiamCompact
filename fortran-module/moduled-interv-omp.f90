subroutine calc_intermediate_v(MX, MY, MZ, NX, NY, NZ, ALPHA, DT, REI, U, V, W, &
                               UD, VD, WD, UU, VV, WW, YJA, SGS, & 
                               GX, EY, TX, TY, TZ, C1, C2, C3, C5, C6, C7, C8, C9, &
                               nID, ierr)

!$ use omp_lib
!  use params
  implicit none

  include 'cpm_fparam.fi'

  integer i, j, k
  integer MX, MY, MZ
  integer NX, NY, NZ

  real(8) ALPHA, DT, REI
  real(8) DABSCU, DABSCV, DABSCW
  real(8) DSGSDX, DSGSDY, DSGSDZ
  real(8) ADVX, ADVY, ADVZ, VISX, VISY, VISZ
  real(8) DUDY, DUDZ
  real(8) DVDX, DVDZ
  real(8) DWDX, DWDY

  real(8) U(0:MX+1,0:MY+1,0:MZ+1), V(0:MX+1,0:MY+1,0:MZ+1), W(0:MX+1,0:MY+1,0:MZ+1)
  real(8) UD(-1:MX+2,-1:MY+2,-1:MZ+2), VD(-1:MX+2,-1:MY+2,-1:MZ+2), WD(-1:MX+2,-1:MY+2,-1:MZ+2)
  real(8) UU(0:MX+1,0:MY+1,0:MZ+1), VV(0:MX+1,0:MY+1,0:MZ+1), WW(0:MX+1,0:MY+1,0:MZ+1)
  real(8) GX(0:MX+1,0:MY+1,0:MZ+1), EY(0:MX+1,0:MY+1,0:MZ+1)
  real(8) TX(0:MX+1,0:MY+1,0:MZ+1), TY(0:MX+1,0:MY+1,0:MZ+1), TZ(0:MX+1,0:MY+1,0:MZ+1)
  real(8) YJA(0:MX+1,0:MY+1,0:MZ+1), SGS(0:MX+1,0:MY+1,0:MZ+1)
  real(8) C1(MX,MY,MZ), C2(MX,MY,MZ), C3(MX,MY,MZ), C5(MX,MY,MZ), C6(MX,MY,MZ)
  real(8) C7(MX,MY,MZ), C8(MX,MY,MZ), C9(MX,MY,MZ)
!
  integer nID(0:5)
  integer ierr, pg
  integer ist, jst, kst
  integer iend, jend, kend

  ierr=0
!
!=========================================
!....<<INTERMEDIATE PHYSICAL VELOCITY>>
!=========================================
!....<<STORE PREVIOUS VELOCITY>>

!$OMP PARALLEL
!$OMP DO &
!$OMP FIRSTPRIVATE(NX, NY, NZ) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = 0, NZ+1
  do j = 0, NY+1
  do i = 0, NX+1
     UD(i,j,k) = U(i,j,k)
     VD(i,j,k) = V(i,j,k)
     WD(i,j,k) = W(i,j,k)
  enddo
  enddo
  enddo
!$OMP END DO

! UV,VD,WDの袖通信をする
!$OMP MASTER
  pg = 0
  call cpm_BndCommS3D(UD,MX,MY,MZ,2,2,CPM_DOUBLE,pg,ierr)
  if( ierr.ne.0 ) then
    write(*,*) "Error at UD Comm Band Cell"
  endif 
  call cpm_BndCommS3D(VD,MX,MY,MZ,2,2,CPM_DOUBLE,pg,ierr)
  if( ierr.ne.0 ) then
    write(*,*) "Error at VD Comm Band Cell"
  endif 
  call cpm_BndCommS3D(WD,MX,MY,MZ,2,2,CPM_DOUBLE,pg,ierr)
  if( ierr.ne.0 ) then
    write(*,*) "Error at WD Comm Band Cell"
  endif 
!$OMP END MASTER
!$OMP BARRIER

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

!
!....<<INNER POINTS>>
!$OMP DO &
!$OMP PRIVATE(DABSCU, DABSCV, DABSCW, ADVX, ADVY, ADVZ) &
!$OMP PRIVATE(VISX, VISY, VISZ, DVDX, DWDX, DSGSDX) &
!$OMP PRIVATE(DUDY, DWDY, DSGSDY, DUDZ, DVDZ, DSGSDZ) &
!$OMP SCHEDULE(static) COLLAPSE(2)
  do k = kst, kend
  do j = jst, jend
  do i = ist, iend
!....<<CONVECTIVE TERMS(3RD-ORDER UPWIND SCHEME)>>
    DABSCU = dabs(UD(i,j,k) * GX(i,j,k))
    DABSCV = dabs(VD(i,j,k) * EY(i,j,k))
    DABSCW = dabs(UD(i,j,k) * TX(i,j,k) + VD(i,j,k) * TY(i,j,k) + WD(i,j,k) * TZ(i,j,k))
    ADVX = .5D0 / 24.D0 *  &
       ((UU(i-1,j,k)  * (-UD(i+1,j,k) + 27.D0 * (UD(i,j,k)   - UD(i-1,j,k)) + UD(i-2,j,k)) &
       +  UU(i,j,k)   * (-UD(i+2,j,k) + 27.D0 * (UD(i+1,j,k) - UD(i,j,k))   + UD(i-1,j,k))) &
       + (VV(i,j-1,k) * (-UD(i,j+1,k) + 27.D0 * (UD(i,j,k)   - UD(i,j-1,k)) + UD(i,j-2,k)) &
       +  VV(i,j,k)   * (-UD(i,j+2,k) + 27.D0 * (UD(i,j+1,k) - UD(i,j,k))   + UD(i,j-1,k))) &
       + (WW(i,j,k-1) * (-UD(i,j,k+1) + 27.D0 * (UD(i,j,k)   - UD(i,j,k-1)) + UD(i,j,k-2)) &
       +  WW(i,j,k)   * (-UD(i,j,k+2) + 27.D0 * (UD(i,j,k+1) - UD(i,j,k))   + UD(i,j,k-1)))) / YJA(i,j,k) &
       + ALPHA * (DABSCU * (UD(i+2,j,k) - 4.D0 *  (UD(i+1,j,k) + UD(i-1,j,k)) &
       + UD(i-2,j,k) + 6.D0 * UD(i,j,k)) + DABSCV*(UD(i,j+2,k) - 4.D0 * (UD(i,j+1,k) + UD(i,j-1,k)) &
       + UD(i,j-2,k) + 6.D0 * UD(i,j,k)) + DABSCW*(UD(i,j,k+2) - 4.D0 * (UD(i,j,k+1) + UD(i,j,k-1)) &
       + UD(i,j,k-2) + 6.D0 * UD(i,j,k))) / 12.D0
    ADVY = .5D0 / 24.D0 *  &
       ((UU(i-1,j,k) * (-VD(i+1,j,k) + 27.D0 * (VD(i,j,k)   - VD(i-1,j,k)) + VD(i-2,j,k)) &
       +  UU(i,j,k)   * (-VD(i+2,j,k) + 27.D0 * (VD(i+1,j,k) - VD(i,j,k))   + VD(i-1,j,k))) &
       + (VV(i,j-1,k) * (-VD(i,j+1,k) + 27.D0 * (VD(i,j,k)   - VD(i,j-1,k)) + VD(i,j-2,k)) &
       +  VV(i,j,k)   * (-VD(i,j+2,k) + 27.D0 * (VD(i,j+1,k) - VD(i,j,k))   + VD(i,j-1,k))) &
       + (WW(i,j,k-1) * (-VD(i,j,k+1) + 27.D0 * (VD(i,j,k)   - VD(i,j,k-1)) + VD(i,j,k-2)) &
       +  WW(i,j,k)   * (-VD(i,j,k+2) + 27.D0 * (VD(i,j,k+1) - VD(i,j,k))   + VD(i,j,k-1)))) / YJA(i,j,k) &
       + ALPHA * (DABSCU * (VD(i+2,j,k) - 4.D0 * (VD(i+1,j,k) + VD(i-1,j,k)) &
       + VD(i-2,j,k)  + 6.D0 * VD(i,j,k)) + DABSCV * (VD(i,j+2,k) - 4.D0 * (VD(i,j+1,k) + VD(i,j-1,k)) &
       + VD(i,j-2,k)  + 6.D0 * VD(i,j,k)) + DABSCW * (VD(i,j,k+2) - 4.D0 * (VD(i,j,k+1) + VD(i,j,k-1)) &
       + VD(i,j,k-2)  + 6.D0 * VD(i,j,k))) / 12.D0
    ADVZ = .5D0 / 24.D0 *  &
       ((UU(i-1,j,k)  * (-WD(i+1,j,k) + 27.D0 * (WD(i,j,k)   - WD(i-1,j,k)) + WD(i-2,j,k)) &
       +  UU(i,j,k)   * (-WD(i+2,j,k) + 27.D0 * (WD(i+1,j,k) - WD(i,j,k))   + WD(i-1,j,k))) &
       + (VV(i,j-1,k) * (-WD(i,j+1,k) + 27.D0 * (WD(i,j,k)   - WD(i,j-1,k)) + WD(i,j-2,k)) &
       +  VV(i,j,k)   * (-WD(i,j+2,k) + 27.D0 * (WD(i,j+1,k) - WD(i,j,k))   + WD(i,j-1,k))) &
       + (WW(i,j,k-1) * (-WD(i,j,k+1) + 27.D0 * (WD(i,j,k)   - WD(i,j,k-1)) + WD(i,j,k-2)) &
       +  WW(i,j,k)   * (-WD(i,j,k+2) + 27.D0 * (WD(i,j,k+1) - WD(i,j,k))   + WD(i,j,k-1)))) /YJA(i,j,k) &
       + ALPHA * (DABSCU * (WD(i+2,j,k) - 4.D0 * (WD(i+1,j,k) + WD(i-1,j,k)) &
       + WD(i-2,j,k ) + 6.D0 * WD(i,j,k)) + DABSCV * (WD(i,j+2,k) - 4.D0 * (WD(i,j+1,k) + WD(i,j-1,k)) &
       + WD(i,j-2,k)  + 6.D0 * WD(i,j,k)) + DABSCW * (WD(i,j,k+2) - 4.D0 * (WD(i,j,k+1) + WD(i,j,k-1)) &
       + WD(i,j,k-2)  + 6.D0 * WD(i,j,k))) / 12.D0
!....<<VISCOUS TERMS OF LAPLACIAN PART>>
    VISX = (C1(i,j,k) * (UD(i+1,j,k) - 2.D0 * UD(i,j,k) + UD(i-1,j,k)) &
       + C2(i,j,k) * (UD(i,j+1,k) - 2.D0 * UD(i,j,k) + UD(i,j-1,k)) &
       + C3(i,j,k) * (UD(i,j,k+1) - 2.D0 * UD(i,j,k) + UD(i,j,k-1)) &
       + C5(i,j,k) * (UD(i,j+1,k+1) - UD(i,j+1,k-1) - UD(i,j-1,k+1) + UD(i,j-1,k-1)) * .25D0 &
       + C6(i,j,k) * (UD(i+1,j,k+1) - UD(i+1,j,k-1) - UD(i-1,j,k+1) + UD(i-1,j,k-1)) * .25D0 &
       + C7(i,j,k) * (UD(i+1,j,k)   - UD(i-1,j,k))*.5D0 &
       + C8(i,j,k) * (UD(i,j+1,k) - UD(i,j-1,k)) * .5D0 &
       + C9(i,j,k) * (UD(i,j,k+1) - UD(i,j,k-1)) * .5D0) * (REI + SGS(i,j,k))
    VISY = (C1(i,j,k) * (VD(i+1,j,k) - 2.D0 * VD(i,j,k) + VD(i-1,j,k)) &
       + C2(i,j,k) * (VD(i,j+1,k) - 2.D0 * VD(i,j,k) + VD(i,j-1,k)) &
       + C3(i,j,k) * (VD(i,j,k+1) - 2.D0 * VD(i,j,k) + VD(i,j,k-1)) &
       + C5(i,j,k) * (VD(i,j+1,k+1) - VD(i,j+1,k-1) - VD(i,j-1,k+1) + VD(i,j-1,k-1)) * .25D0 &
       + C6(i,j,k) * (VD(i+1,j,k+1) - VD(i+1,j,k-1) - VD(i-1,j,k+1) + VD(i-1,j,k-1)) * .25D0 &
       + C7(i,j,k) * (VD(i+1,j,k) - VD(i-1,j,k)) * .5D0 &
       + C8(i,j,k) * (VD(i,j+1,k) - VD(i,j-1,k)) * .5D0 &
       + C9(i,j,k) * (VD(i,j,k+1) - VD(i,j,k-1)) * .5D0) * (REI+SGS(i,j,k))
    VISZ = (C1(i,j,k) * (WD(i+1,j,k) - 2.D0 * WD(i,j,k) + WD(i-1,j,k)) &
       + C2(i,j,k) * (WD(i,j+1,k) - 2.D0 * WD(i,j,k) + WD(i,j-1,k)) &
       + C3(i,j,k) * (WD(i,j,k+1) - 2.D0 * WD(i,j,k) + WD(i,j,k-1)) &
       + C5(i,j,k) * (WD(i,j+1,k+1) - WD(i,j+1,k-1) - WD(i,j-1,k+1) + WD(i,j-1,k-1)) * .25D0 &
       + C6(i,j,k) * (WD(i+1,j,k+1) - WD(i+1,j,k-1) - WD(i-1,j,k+1) + WD(i-1,j,k-1)) * .25D0 &
       + C7(i,j,k) * (WD(i+1,j,k) - WD(i-1,j,k)) * .5D0 &
       + C8(i,j,k) * (WD(i,j+1,k) - WD(i,j-1,k)) * .5D0 &
       + C9(i,j,k) * (WD(i,j,k+1) - WD(i,j,k-1)) * .5D0) * (REI + SGS(i,j,k))
!....<<RESIDUAL SGS VISCOSITY TERMS>>
    DVDX = .5D0 * (VD(i+1,j,k) - VD(i-1,j,k)) * GX(i,j,k) + .5D0 * (VD(i,j,k+1) - VD(i,j,k-1)) * TX(i,j,k)
    DWDX = .5D0 * (WD(i+1,j,k) - WD(i-1,j,k)) * GX(i,j,k) + .5D0 * (WD(i,j,k+1) - WD(i,j,k-1)) * TX(i,j,k)
    DSGSDX = .5D0 * (SGS(i+1,j,k) - SGS(i-1,j,k)) * GX(i,j,k) + .5D0 * (SGS(i,j,k+1) - SGS(i,j,k-1)) * TX(i,j,k)
    DUDY = .5D0 * (UD(i,j+1,k) - UD(i,j-1,k)) * EY(i,j,k) + .5D0 * (UD(i,j,k+1)-UD(i,j,k-1)) * TY(i,j,k)
    DWDY = .5D0 * (WD(i,j+1,k) - WD(i,j-1,k)) * EY(i,j,k) + .5D0 * (WD(i,j,k+1) - WD(i,j,k-1)) * TY(i,j,k)
    DSGSDY =.5D0 * (SGS(i,j+1,k) - SGS(i,j-1,k)) * EY(i,j,k) + .5D0 * (SGS(i,j,k+1) - SGS(i,j,k-1)) * TY(i,j,k)
    DUDZ = .5D0 * (UD(i,j,k+1) - UD(i,j,k-1)) * TZ(i,j,k)
    DVDZ = .5D0 * (VD(i,j,k+1) - VD(i,j,k-1)) * TZ(i,j,k)
    DSGSDZ= .5D0 * (SGS(i,j,k+1) - SGS(i,j,k-1)) * TZ(i,j,k)
!....<<INTERMEDIATE PHYSICAL VELOCITY>>
    U(i,j,k) = UD(i,j,k) + DT * (-ADVX + VISX &
       + DSGSDX * (2.D0 * ((UD(i+1,j,k) - UD(i-1,j,k)) * .5D0 * GX(i,j,k) &
       +  (UD(i,j,k+1) - UD(i,j,k-1)) * .5D0 * TX(i,j,k))) &
       + DSGSDY * (DUDY + DVDX) &
       + DSGSDZ * (DUDZ + DWDX))
    V(i,j,k) = VD(i,j,k) + DT * (-ADVY + VISY &
       + DSGSDX * ( DVDX + DUDY ) &
       + DSGSDY * (2.D0 * ((VD(i,j+1,k) - VD(i,j-1,k)) * .5D0 * EY(i,j,k) &
       +  (VD(i,j,k+1) - VD(i,j,k-1)) * .5D0 * TY(i,j,k))) &
       + DSGSDZ * (DVDZ + DWDY))
    W(i,j,k) = WD(i,j,k) + DT * (-ADVZ + VISZ &
       + DSGSDX * (DWDX + DUDZ) &
       + DSGSDY * (DWDY + DVDZ) &
       + DSGSDZ * (2.D0*(WD(i,j,k+1) - WD(i,j,k-1)) * .5D0 * TZ(i,j,k)))
  enddo
  enddo
  enddo
!$OMP END DO
!$OMP END PARALLEL

!
end subroutine calc_intermediate_v
