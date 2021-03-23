MODULE TRLTOM_MOD
  CONTAINS
  SUBROUTINE TRLTOM(PFBUF_IN,PFBUF,KFIELD)
  
  !**** *TRLTOM * - transposition in Fourierspace
  
  !     Purpose.
  !     --------
  !              Transpose Fourier coefficients from partitioning
  !              over latitudes to partitioning over wave numbers
  !              This is done between inverse Legendre Transform
  !              and inverse FFT.
  !              This is the inverse routine of TRMTOL.
  
  !**   Interface.
  !     ----------
  !        *CALL* *TRLTOM(...)*
  
  !        Explicit arguments : PFBUF  - Fourier coefficient buffer. It is
  !        --------------------          used for both input and output.
  
  !                             KFIELD - Number of fields communicated
  
  !        Implicit arguments :
  !        --------------------
  
  !     Method.
  !     -------
  !        See documentation
  
  !     Externals.
  !     ----------
  
  !     Reference.
  !     ----------
  !        ECMWF Research Department documentation of the IFS
  
  !     Author.
  !     -------
  !        MPP Group *ECMWF*
  
  !     Modifications.
  !     --------------
  !        Original : 95-10-01
  !        Modified : 97-06-18 G. Mozdzynski - control MPI mailbox use
  !                                            (NCOMBFLEN) for nphase.eq.1
  !        Modified : 99-05-28  D.Salmond - Optimise copies.
  !        Modified : 00-02-02  M.Hamrud  - Remove NPHASE
  !        D.Salmond : 01-11-23 LIMP_NOOLAP Option for non-overlapping message
  !                             passing and buffer packing
  !        G.Mozdzynski : 08-01-01 Cleanup
  !        Y.Seity   : 07-08-30 Add barrier synchonisation under LSYNC_TRANS
  !     ------------------------------------------------------------------
  
  USE PARKIND1  ,ONLY : JPIM     ,JPRBT, JPRB
  USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK
  
  USE MPL_MODULE  ,ONLY : MPL_ALLTOALLV, MPL_BARRIER, MPL_ALL_MS_COMM, MPL_MYRANK, MPL_WAIT, JP_NON_BLOCKING_STANDARD
  
  USE TPM_DISTR       ,ONLY : D, MTAGLM, MYSETW, NPRTRW, NPROC, MYPROC
  USE TPM_GEN         ,ONLY : LSYNC_TRANS
  
  USE MPI
  
  !USE SET2PE_MOD
  !USE MYSENDSET_MOD
  !USE MYRECVSET_MOD
  !USE ABORT_TRANS_MOD
  !
  
  IMPLICIT NONE
  
  
  INTERFACE
  
    FUNCTION ALLTOALLV_CUDAIPC(input,len,soff,output,roff,mtol_or_ltom) BIND(C,name='Alltoallv_CUDAIPC')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      real(c_double), dimension(*) :: input,output
      integer(c_int), dimension(*) :: len,soff,roff
      integer(c_int),value :: mtol_or_ltom
      integer(c_int) :: ALLTOALLV_CUDAIPC
    END FUNCTION ALLTOALLV_CUDAIPC
  
  END INTERFACE
  
  
  INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
  REAL(KIND=JPRBT)   ,INTENT(INOUT)  :: PFBUF(:)
  REAL(KIND=JPRBT)   ,INTENT(INOUT)  :: PFBUF_IN(:)
  
  INTEGER(KIND=JPIM) :: ILENS(NPRTRW),IOFFS(NPRTRW),ILENR(NPRTRW),IOFFR(NPRTRW)
  
  INTEGER(KIND=JPIM) :: ITAG, J, ILEN, ISTA
  
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  REAL(KIND=JPRB) :: ZHOOK_HANDLE_BAR
  REAL(KIND=JPRB) :: ZHOOK_HANDLE_BAR2
  
  REAL(KIND=JPRBT)    :: ZDUM(1)
  INTEGER(KIND=JPIM) :: IREQ
  INTEGER(KIND=JPIM) :: IERROR
  !     ------------------------------------------------------------------
  
  REAL(KIND=JPRBT) :: T1, T2, TIMEF, tc
  INTEGER(KIND=JPIM) :: MTOL_OR_LTOM, NOFULLPEERACCESS
  INTEGER(KIND=JPIM) :: IRANK,iunit

  IF (LHOOK) CALL DR_HOOK('TRLTOM',0,ZHOOK_HANDLE)
  
  ITAG = MTAGLM
  
  DO J=1,NPRTRW
    ILENS(J) = D%NLTSGTB(J)*KFIELD
    IOFFS(J) = D%NSTAGT1B(D%MSTABF(J))*KFIELD
    ILENR(J) = D%NLTSFTB(J)*KFIELD
    IOFFR(J) = D%NSTAGT1B(J)*KFIELD
  ENDDO
  
  IF(NPROC > 1) THEN
    CALL GSTATS(806,0)
  
 !!! #ifdef USE_CUDA_AWARE_MPI
    !IERROR=0
    !NOFULLPEERACCESS=0
    !MTOL_OR_LTOM=1 ! 0 if called from TRMTOL, 1 from TRLTOM
    !NOFULLPEERACCESS=ALLTOALLV_CUDAIPC(PFBUF_IN,ILENS,IOFFS,PFBUF,IOFFR,MTOL_OR_LTOM)
    !if (NOFULLPEERACCESS) then
    !   CALL MPI_ALLTOALLV(PFBUF_IN,ILENS,IOFFS,MPI_DOUBLE_PRECISION, &
    !        & PFBUF,ILENR,IOFFR,MPI_DOUBLE_PRECISION,MPL_ALL_MS_COMM,IERROR)
    !end if
    !!!!$ACC end host_data
  
  !!#else
  
    !$ACC update host(PFBUF_IN)
    CALL MPL_ALLTOALLV(PSENDBUF=PFBUF_IN,KSENDCOUNTS=ILENS,&
     & PRECVBUF=PFBUF,KRECVCOUNTS=ILENR,KSENDDISPL=IOFFS,KRECVDISPL=IOFFR,&
     & KCOMM=MPL_ALL_MS_COMM,CDSTRING='TRLTOM:')
    !$ACC update device(PFBUF)

    call MPI_BARRIER(MPI_COMM_WORLD,IERROR)
   ! debug
   !IF( MYPROC==1 ) THEN
   !  iunit=300+myproc
     ! what is sent by proc==2 
    !     ISTA=IOFFR(2)
    !     ILEN=ILENR(2)
    !     DO J=ISTA,ISTA+ILEN-1
    !       write(iunit,*) 'mpi from proc 2 ',J,ISTA,ILEN,PFBUF(J)
    !     ENDDO
    !     ILEN = D%NLTSGTB(MYSETW)*KFIELD
    !     ISTA = D%NSTAGT1B(MYSETW)*KFIELD+1
    !     DO J=ISTA,ISTA+ILEN-1
    !       write(iunit,*) 'mpi from proc 1 ',J,ISTA,ILEN,PFBUF(J)
    !     ENDDO
   !ENDIF

  !!#endif
  
  CALL GSTATS(806,1)
  ELSE

    ILEN = D%NLTSGTB(MYSETW)*KFIELD
    ISTA = D%NSTAGT1B(MYSETW)*KFIELD+1

#ifdef UNDEF
WRITE (*, *) __FILE__, ':', __LINE__ 
PRINT *, " ISTA = ", ISTA, " ILEN = ", ILEN
!$acc parallel copyin (ISTA, ILEN) num_gangs(1) num_workers(1) vector_length(1) present(PFBUF_IN)
PRINT *, " ISTA = ", ISTA, " ILEN = ", ILEN
DO J = ISTA, ISTA+ILEN-1
IF (ABS (PFBUF_IN (J)) < 1E-15) THEN
PRINT *, J, 0._8
ELSE
PRINT *, J, PFBUF_IN (J)
ENDIF
ENDDO
!$acc end parallel
#endif


    CALL GSTATS(1607,0)
    !$ACC data present(PFBUF_IN,PFBUF)
    !$ACC parallel loop
    DO J=ISTA,ISTA+ILEN-1
      PFBUF(J) = PFBUF_IN(J)
    ENDDO
    !$ACC end data

#ifdef UNDEF
WRITE (*, *) __FILE__, ':', __LINE__ 
!$acc parallel copyin (ISTA, ILEN) num_gangs(1) num_workers(1) vector_length(1) present(PFBUF)
DO J = ISTA, ISTA+ILEN-1
IF (ABS (PFBUF (J)) < 1E-15) THEN
PRINT *, J, 0._8
ELSE
PRINT *, J, PFBUF (J)
ENDIF
ENDDO
!$acc end parallel
#endif

    CALL GSTATS(1607,1)
  ENDIF
  
  IF (LHOOK) CALL DR_HOOK('TRLTOM',1,ZHOOK_HANDLE)
  !     ------------------------------------------------------------------
  END SUBROUTINE TRLTOM
  END MODULE TRLTOM_MOD
