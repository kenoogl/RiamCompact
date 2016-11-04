subroutine init(MX, MY, MZ, NX, NY, NZ, ICON, XLEN, YLEN, ZLEN, TIME, RE, X, Y, Z, &
                U, V, W, P, UU, VV, WW, SGS, YJA, GX, EY, TX, TY, TZ, C1, C2, C3, &
                C5, C6, C7, C8, C9, D1, D2, D3, D5, D6, ORG, GSIZE, head, nID, ierr )
!
!  use params
  implicit none
!
  include 'cpm_fparam.fi'
!
  integer i, j, k
  integer MX, MY, MZ
  integer NX, NY, NZ, ICON
  real(8) XLEN, YLEN, ZLEN
  real(8) TIME, RE
  real(8) XG, YE, ZE, ZG, ZT
  real(8) XGG, YEE, ZEE, ZET, ZGG, ZGT, ZTT
  real(8) GXX, EYY, TXX, TYY, TZZ
  real(8) YJGA, YJEA, YJTA, YJAI
  real(8) X(0:MX+1,0:MY+1,0:MZ+1), Y(0:MX+1,0:MY+1,0:MZ+1), Z(0:MX+1,0:MY+1,0:MZ+1)
  real(8) U(0:MX+1,0:MY+1,0:MZ+1), V(0:MX+1,0:MY+1,0:MZ+1), W(0:MX+1,0:MY+1,0:MZ+1)
  real(8) P(0:MX+1,0:MY+1,0:MZ+1), SGS(0:MX+1,0:MY+1,0:MZ+1), YJA(0:MX+1,0:MY+1,0:MZ+1)
  real(8) UU(0:MX+1,0:MY+1,0:MZ+1), VV(0:MX+1,0:MY+1,0:MZ+1), WW(0:MX+1,0:MY+1,0:MZ+1)
  real(8) GX(0:MX+1,0:MY+1,0:MZ+1), EY(0:MX+1,0:MY+1,0:MZ+1)
  real(8) TX(0:MX+1,0:MY+1,0:MZ+1), TY(0:MX+1,0:MY+1,0:MZ+1), TZ(0:MX+1,0:MY+1,0:MZ+1)
  real(8) C1(MX,MY,MZ), C2(MX,MY,MZ), C3(MX,MY,MZ), C5(MX,MY,MZ), C6(MX,MY,MZ)
  real(8) C7(MX,MY,MZ), C8(MX,MY,MZ), C9(MX,MY,MZ)
  real(8) D1(0:MX+1,0:MY+1,0:MZ+1), D2(0:MX+1,0:MY+1,0:MZ+1), D3(0:MX+1,0:MY+1,0:MZ+1)
  real(8) D5(0:MX+1,0:MY+1,0:MZ+1), D6(0:MX+1,0:MY+1,0:MZ+1)

  real(8) ORG(3)
  integer GSIZE(3), head(0:2), nID(0:5)

  integer myRank,pg,ierr
  integer NX2, NY2, NZ2
  character(len=256) fname

  ierr=0
  pg=0

!....<<3D-GRID DATA MAKING>>
  do k = 0, NZ+1
  do j = 0, NY+1
  do i = 0, NX+1
    X(i,j,k) = ORG(1) + (XLEN / dfloat(GSIZE(1)-1) * dfloat(i-1+head(X_DIR)))
    Y(i,j,k) = ORG(2) + (YLEN / dfloat(GSIZE(2)-1) * dfloat(j-1+head(Y_DIR)))
    Z(i,j,k) = ORG(3) + (ZLEN / dfloat(GSIZE(3)-1) * dfloat(k-1+head(Z_DIR)))
  enddo
  enddo
  enddo

!....<<ZERO CLEAR>>
  U(:,:,:) = 1.D0
  V(:,:,:) = 0.D0
  W(:,:,:) = 0.D0
  P(:,:,:) = 0.D0
  UU(:,:,:) = 0.D0
  VV(:,:,:) = 0.D0
  WW(:,:,:) = 0.D0
  SGS(:,:,:) = 0.D0

!....<<CONTINUATION CALCULATION>>
  if (ICON == 1) then

! get my rank No.
    call cpm_GetMyRankID(myrank, pg, ierr)

    write(fname,'(a,i6.6,a)') 'uvwp_id',myrank,'.dat'
    open(2,FILE=fname,FORM='FORMATTED',STATUS='OLD',ERR=99)
    read(2,'(3I5,2F15.7)')NX2,NY2,NZ2,TIME,RE

    if( NX==NX2 .and. NY==NY2 .and. NZ==NZ2 ) then
      read(2,'(3F15.7)')(((U(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
      read(2,'(3F15.7)')(((V(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
      read(2,'(3F15.7)')(((W(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
      read(2,'(3F15.7)')(((P(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
      read(2,'(3F15.7)')(((UU(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
      read(2,'(3F15.7)')(((VV(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
      read(2,'(3F15.7)')(((WW(i,j,k),i=0,NX+1),j=0,NY+1),k=0,NZ+1)
    else
      write(*,*) "Error Ther number of grids discords."
      if( ierr == 0 ) ierr=1
    endif
    close(2)
  endif

!....<<EXTRAPOLATION>>
! -Y方向にランクがないとき境界の座標値計算
  if( nID(Y_MINUS) < 0 ) then
    do k = 1, NZ
    do i = 1, NX
      X(i,0,k) = 2.D0 * X(i,1,k) - X(i,2,k)
      Y(i,0,k) = 2.D0 * Y(i,1,k) - Y(i,2,k)
      Z(i,0,k) = 2.D0 * Z(i,1,k) - Z(i,2,k)
    enddo
    enddo
  endif
! +Y方向にランクがないとき境界の座標値計算
  if( nID(Y_PLUS) < 0 ) then
    do k = 1, NZ
    do i = 1, NX
      X(i,NY+1,k) = 2.D0 * X(i,NY,k) - X(i,NY-1,k)
      Y(i,NY+1,k) = 2.D0 * Y(i,NY,k) - Y(i,NY-1,k)
      Z(i,NY+1,k) = 2.D0 * Z(i,NY,k) - Z(i,NY-1,k)
    enddo
    enddo
  endif
!
! -X方向にランクがないとき境界の座標値計算
  if( nID(X_MINUS) < 0 ) then
    do k = 1, NZ
    do j = 0, NY+1
      X(0,j,k) = 2.D0 * X(1,j,k) - X(2,j,k)
      Y(0,j,k) = 2.D0 * Y(1,j,k) - Y(2,j,k)
      Z(0,j,k) = 2.D0 * Z(1,j,k) - Z(2,j,k)
    enddo
    enddo
  endif
! +X方向にランクがないとき境界の座標値計算
  if( nID(X_PLUS) < 0 ) then
    do k = 1, NZ
    do j = 0, NY+1
      X(NX+1,j,k) = 2.D0 * X(NX,j,k) - X(NX-1,j,k)
      Y(NX+1,j,k) = 2.D0 * Y(NX,j,k) - Y(NX-1,j,k)
      Z(NX+1,j,k) = 2.D0 * Z(NX,j,k) - Z(NX-1,j,k)
    enddo
    enddo
  endif
!
! -Z方向にランクがないとき境界の座標値計算
  if( nID(Z_MINUS) < 0 ) then
    do j = 0, NY+1
    do i = 0, NX+1
      X(i,j,0) = 2.D0 * X(i,j,1) - X(i,j,2)
      Y(i,j,0) = 2.D0 * Y(i,j,1) - Y(i,j,2)
      Z(i,j,0) = 2.D0 * Z(i,j,1) - Z(i,j,2)
    enddo
    enddo
  endif
!
! +Z方向にランクがないとき境界の座標値計算
  if( nID(Z_PLUS) < 0 ) then
    do j = 0, NY+1
    do i = 0, NX+1
      X(i,j,NZ+1) = 2.D0 * X(i,j,NZ) - X(i,j,NZ-1)
      Y(i,j,NZ+1) = 2.D0 * Y(i,j,NZ) - Y(i,j,NZ-1)
      Z(i,j,NZ+1) = 2.D0 * Z(i,j,NZ) - Z(i,j,NZ-1)
    enddo
    enddo
  endif
!
!....<<METRICS ETC>>
  do k = 1, NZ
    do j = 1, NY
      do i = 1, NX
        XG = .5D0 * (X(i+1,j,k) - X(i-1,j,k))
        ZG = .5D0 * (Z(i+1,j,k) - Z(i-1,j,k))
        YE = .5D0 * (Y(i,j+1,k) - Y(i,j-1,k))
        ZE = .5D0 * (Z(i,j+1,k) - Z(i,j-1,k))
        ZT = .5D0 * (Z(i,j,k+1) - Z(i,j,k-1))
        XGG = X(i+1,j,k) - 2.D0 * X(i,j,k) + X(i-1,j,k)
        ZGG = Z(i+1,j,k) - 2.D0 * Z(i,j,k) + Z(i-1,j,k)
        YEE = Y(i,j+1,k) - 2.D0 * Y(i,j,k) + Y(i,j-1,k)
        ZEE = Z(i,j+1,k) - 2.D0 * Z(i,j,k) + Z(i,j-1,k)
        ZTT = Z(i,j,k+1) - 2.D0 * Z(i,j,k) + Z(i,j,k-1)
        ZET = .25D0 * (Z(i,j+1,k+1) - Z(i,j+1,k-1) - Z(i,j-1,k+1) + Z(i,j-1,k-1))
        ZGT = .25D0 * (Z(i+1,j,k+1) - Z(i+1,j,k-1) - Z(i-1,j,k+1) + Z(i-1,j,k-1))
!....<<JACOBIAN>>
        YJA(i,j,k) = XG * YE * ZT
        YJGA = XGG * YE * ZT + ZGT * XG * YE
        YJEA = YEE * ZT * XG + ZET * XG * YE
        YJTA = ZTT * XG * YE
        YJAI = 1.D0 / YJA(i,j,k)
!....<<1ST-ORDER METRICS>>
        GX(i,j,k) =  YE * ZT * YJAI
        EY(i,j,k) =  XG * ZT * YJAI
        TX(i,j,k) = -YE * ZG * YJAI
        TY(i,j,k) = -XG * ZE * YJAI
        TZ(i,j,k) =  XG * YE * YJAI
!....<<2ND-ORDER METRICS>>
        GXX =  YE * ZGT * GX(i,j,k) * YJAI - YJGA * GX(i,j,k) * GX(i,j,k) * YJAI &
             + YE * ZTT * TX(i,j,k) * YJAI - YJTA * GX(i,j,k) * TX(i,j,k) * YJAI
        EYY =  XG * ZET * EY(i,j,k) * YJAI - YJEA * EY(i,j,k) * EY(i,j,k) * YJAI &
             + XG * ZTT * TY(i,j,k) * YJAI - YJTA * EY(i,j,k) * TY(i,j,k) * YJAI
        TXX = -YE * ZGG * GX(i,j,k) * YJAI - YJGA * GX(i,j,k) * TX(i,j,k) * YJAI &
              -YE * ZGT * TX(i,j,k) * YJAI - YJTA * TX(i,j,k) * TX(i,j,k) * YJAI
        TYY = -XG * ZEE * EY(i,j,k) * YJAI - YJEA * EY(i,j,k) * TY(i,j,k) * YJAI &
              -XG * ZET * TY(i,j,k) * YJAI - YJTA * TY(i,j,k) * TY(i,j,k) * YJAI
        TZZ=                               - YJTA * TZ(i,j,k) * TZ(i,j,k) * YJAI
!....<<LAPLACIAN COEFFICIENT>>
        C1(i,j,k) =  GX(i,j,k) ** 2
        C2(i,j,k) =  EY(i,j,k) ** 2
        C3(i,j,k) =  TX(i,j,k) ** 2 + TY(i,j,k) ** 2 + TZ(i,j,k) ** 2
        C5(i,j,k) = (EY(i,j,k) * TY(i,j,k)) * 2.D0
        C6(i,j,k) = (TX(i,j,k) * GX(i,j,k)) * 2.D0
        C7(i,j,k) =  GXX
        C8(i,j,k) =  EYY
        C9(i,j,k) =  TXX + TYY + TZZ
!....<<PRESSURE GRADIENT COEFFICIENT>>
        D1(i,j,k) = C1(i,j,k) * YJA(i,j,k)
        D2(i,j,k) = C2(i,j,k) * YJA(i,j,k)
        D3(i,j,k) = C3(i,j,k) * YJA(i,j,k)
        D5(i,j,k) = C5(i,j,k) * YJA(i,j,k) * .5D0
        D6(i,j,k) = C6(i,j,k) * YJA(i,j,k) * .5D0
      enddo
    enddo
  enddo
!
  return
!
99 continue
  write(*,'(2a)') " Error restart file open error : ",fname
  if( ierr == 0 ) ierr=1
  return

end subroutine init
