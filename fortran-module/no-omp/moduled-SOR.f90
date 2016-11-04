subroutine sor(MX, MY, MZ, NX, NY, NZ, ILAP, ITER, NOXYZ, OMEGA, EPS, RMSP, P, C1, C2, C3, C5, C6, C7, C8, C9, RHS)
!
!  use params
  implicit none
!
  integer i, j, k
  integer MX, MY, MZ
  integer NX, NY, NZ
  integer ILAP, ITER, NOXYZ

  real(8) ERR, RP, OMEGA, EPS, RMSP
  real(8) P(0:MX+1,0:MY+1,0:MZ+1), RHS(MX,MY,MZ)
  real(8) C1(MX,MY,MZ), C2(MX,MY,MZ), C3(MX,MY,MZ), C5(MX,MY,MZ), C6(MX,MY,MZ)
  real(8) C7(MX,MY,MZ), C8(MX,MY,MZ), C9(MX,MY,MZ)

!....<<SOR METHOD>>
  do ILAP = 1, ITER
    ERR = 0.D0
!....<<INNER POINTS>>
    do k = 2, NZ-1
      do j = 2, NY-1
        do i = 2, NX-1
!....<<NEW PRESSURE>>
          RP = .5D0 * ((C1(i,j,k) * (P(i+1,j,k)   + P(i-1,j,k)) &
                    +   C2(i,j,k) * (P(i,j+1,k)   + P(i,j-1,k)) &
                    +   C3(i,j,k) * (P(i,j,k+1)   + P(i,j,k-1)) &
                    +   C5(i,j,k) * (P(i,j+1,k+1) - P(i,j+1,k-1) - P(i,j-1,k+1) + P(i,j-1,k-1)) * .25D0 &
                    +   C6(i,j,k) * (P(i+1,j,k+1) - P(i+1,j,k-1) - P(i-1,j,k+1) + P(i-1,j,k-1)) * .25D0 &
                    +   C7(i,j,k) * (P(i+1,j,k)   - P(i-1,j,k))*.5D0 &
                    +   C8(i,j,k) * (P(i,j+1,k)   - P(i,j-1,k))*.5D0 &
                    +   C9(i,j,k) * (P(i,j,k+1)   - P(i,j,k-1))*.5D0) &
                    + RHS(i,j,k)) / ((C1(i,j,k) + C2(i,j,k) + C3(i,j,k))) - P(i,j,k)
          ERR = ERR + RP * RP
          P(i,j,k) = P(i,j,k) + OMEGA * RP
        enddo
      enddo
    enddo
!....<<CONVERGE OR NOT ?>>

    RMSP = dsqrt( ERR / dfloat(NOXYZ))
    if(RMSP <= EPS) exit
  enddo
  return
!
end subroutine SOR
