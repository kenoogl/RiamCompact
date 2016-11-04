subroutine display(MX, MY, MZ, NX, NY, NZ, ISTEP, NSTEP, ILAP, ICOUNT, NSOR, &
                   U, V, W, &
                   IPROG1, IPROG2, TIME, RE, RMSP, MY_RANK)
!
!  use params
  implicit none
!
  integer i, j, k, ii
  integer MX, MY, MZ
  integer NX, NY, NZ
  integer ISTEP, NSTEP, ILAP, ICOUNT
  integer NSOR, IPROG1, IPROG2
  integer MY_RANK

  real(8) TIME, RE, RMSP
  real(8) A(2)
  real(8) U(0:MX+1,0:MY+1,0:MZ+1), V(0:MX+1,0:MY+1,0:MZ+1), W(0:MX+1,0:MY+1,0:MZ+1)

!
  character(LEN=80) FMT1, FMT2, FMT3, FMT4

  character(len=256) fname

!....<<DISPLAY>> 
  if (mod(ISTEP,NSOR) == 0 .and. MY_RANK == 0 ) then
    FMT1="(' Time=',F8.3,'  Re=',F9.1,' Istep=',I7)"
    FMT2="(' Iteration number of Poisson equation=',I7)"
    write(*,*)'**************************************************'
    write(*,FMT1) TIME,RE,ISTEP
    write(*,FMT2) ILAP
    write(*,*)'Root mean square error=',RMSP
    write(*,*)'**************************************************'
  endif
!
  if (mod(ISTEP,IPROG1) == 0 .and. MY_RANK == 0 ) then
    if ( ISTEP == IPROG1 ) then
        open(10,FILE='cal-log.dat',FORM='FORMATTED')
    else
        open(10,FILE='cal-log.dat',FORM='FORMATTED', position = 'append')
    endif

    FMT3="(' Time=',F8.3,'  Re=',F9.1,' Istep=',I7)"
    FMT4="(' Poisson equation iteration=',I7)"
    write(10,*)'*******************************************'
    write(10,FMT3)TIME,RE,ISTEP
    write(10,FMT4)ILAP
    write(10,*)'RMS error=',RMSP
    write(10,*)'*******************************************'

    close(10)

  endif
!
!....<<SUMMATION OF FIELD-DATA>>
!    do k=1,NZ
!      do j=1,NY
!        do i=1,NX
!          UAVE(i,j,k)=UAVE(i,j,k)+U(i,j,k)
!          VAVE(i,j,k)=VAVE(i,j,k)+V(i,j,k)
!          WAVE(i,j,k)=WAVE(i,j,k)+W(i,j,k)
!        enddo
!      enddo
!    enddo
!
!....<<FLOW VISUALIZATION FOR REAL SCALE INSTANTANEOUS DATA>>
  if (mod(ISTEP,IPROG2) == 0) then

    write(fname,'(a,i6.6,a)') '3d-inst-vis_id',MY_RANK,'.dat'

    if ( ISTEP == IPROG2 ) then
        open(20,FILE=fname,FORM='UNFORMATTED')
    else
        open(20,FILE=fname,FORM='UNFORMATTED', position = 'append')
    endif
!
    ICOUNT=ICOUNT+1
    write(20)ISTEP,sngl(TIME)
    write(20)(sngl(A(ii)),ii=1,2)
    write(20)(((sngl(U(i,j,k)),i=1,NX),j=1,NY),k=1,NZ)
    write(20)(((sngl(V(i,j,k)),i=1,NX),j=1,NY),k=1,NZ)
    write(20)(((sngl(W(i,j,k)),i=1,NX),j=1,NY),k=1,NZ)

    if( MY_RANK == 0 ) then
    write(*,*)'**************************************************'
    write(*,*)'Istep=',ISTEP,'  ,RC-Scope for (U) ==> 3d-inst-vis.dat'
    write(*,*)'Istep=',ISTEP,'  ,RC-Scope for (V) ==> 3d-inst-vis.dat'
    write(*,*)'Istep=',ISTEP,'  ,RC-Scope for (W) ==> 3d-inst-vis.dat'
    write(*,*)'Number =',ICOUNT,'/',NSTEP/IPROG2,'(Total)'
    write(*,*)'**************************************************'
    endif

    close(20)
    close(21)

  endif
!
  return
end subroutine display
