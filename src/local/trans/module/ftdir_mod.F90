MODULE FTDIR_MOD
CONTAINS
SUBROUTINE FTDIR(KFIELDS)


!**** *FTDIR - Direct Fourier transform

!     Purpose. Routine for Grid-point to Fourier transform
!     --------

!**   Interface.
!     ----------
!        CALL FTDIR(..)

!        Explicit arguments :  PREEL   - Fourier/grid-point array
!        --------------------  KFIELDS - number of fields

!     Method.
!     -------

!     Externals.  FFT992 - FFT routine
!     ----------
!

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-03-03
!        G. Radnoti 01-04-24 2D model (NLOEN=1)
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        G. Mozdzynski (Oct 2014): support for FFTW transforms
!        G. Mozdzynski (Jun 2015): Support alternative FFTs to FFTW

!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM, JPIB, JPRBT

USE TPM_DISTR       ,ONLY : D, MYSETW, MYPROC, NPROC,D_NSTAGTF,D_NPTRLS
USE TPM_TRANS       ,ONLY : ZGTF
USE TPM_GEOMETRY    ,ONLY : G,G_NMEN,G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX
USE TPM_FFT         ,ONLY : T, TB
#ifdef WITH_FFTW
USE TPM_FFTW        ,ONLY : TW, EXEC_FFTW
#endif
USE TPM_FFTC        ,ONLY : CREATE_PLAN_FFT
USE TPM_DIM         ,ONLY : R,R_NNOEXTZL
USE CUDA_DEVICE_MOD
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELDS
INTEGER(KIND=JPIM)  :: KGL
!REAL(KIND=JPRBT), INTENT(INOUT) :: PREEL(:,:)

INTEGER(KIND=JPIM) :: IGLG,IST,ILEN,IJUMP,JJ,JF,IST1
INTEGER(KIND=JPIM) :: IOFF,IRLEN,ICLEN, ITYPE
INTEGER(KIND=JPIM) :: IPLAN_R2C
INTEGER(KIND=JPIM) :: JMAX
REAL(KIND=JPRBT)    :: SCAL
LOGICAL :: LL_ALL=.FALSE. ! T=do kfields ffts in one batch, F=do kfields ffts one at a time

INTEGER(KIND=JPIM) :: IBEG,IEND,IINC,ISCAL
INTEGER(KIND=JPIM) :: OFFSET_VAR, IUNIT, ISIZE
integer :: istat

!     ------------------------------------------------------------------

#ifdef UNDEF
BLOCK
  INTEGER :: ISTR, JFLD, JGL
  IF (MODULO (D%NLENGTF, D%NDGL_FS)) & 
    & CALL ABOR1 ('D%NLENGTF IS NOT A MULTIPLE OF D%NDGL_FS')
  ISTR = D%NLENGTF / D%NDGL_FS
  DO JFLD = 1, 1
  WRITE (88, *) " JFLD = ", JFLD
  DO JGL = 1, D%NDGL_FS
    WRITE (88, '(100E12.4)') ZGTF (JFLD,1+(JGL-1)*ISTR:JGL*ISTR)
  ENDDO
  ENDDO
ENDBLOCK
#endif

IF(MYPROC > NPROC/2)THEN
  IBEG=1
  IEND=D%NDGL_FS
  IINC=1
ELSE
  IBEG=D%NDGL_FS
  IEND=1
  IINC=-1
ENDIF

OFFSET_VAR=D_NPTRLS(MYSETW)

! debug
!iunit=myproc+300
!!$ACC update host(ZGTF)
!!$ACC wait
!DO IGLG=IBEG+OFFSET_VAR-1,IEND+OFFSET_VAR-1,IINC
!  DO JJ=1, G_NLOEN_MAX+2
!    JMAX = G_NLOEN(IGLG)+2
!    if (JJ .le. JMAX) then
!      IF (G_NLOEN(IGLG)>1) THEN
!      KGL=IGLG-OFFSET_VAR+1
!      IOFF=D_NSTAGTF(KGL)+1
!      !JF=41+137-1
!      !if( IGLG .EQ. 41 .and. JJ.LE.10 ) write(iunit,*) 'fft_in ',JJ,G_NLOEN(IGLG),IGLG,ZGTF(JF,IOFF+JJ-1)
!!      !if( IGLG .EQ. 41 ) write(iunit,*) 'fft_in ',JJ,G_NLOEN(IGLG),IGLG,ZGTF(1,IOFF+JJ-1)
!!      !if( IGLG .EQ. 43 ) write(iunit,*) 'fft_in ',JJ,G_NLOEN(IGLG),IGLG,ZGTF(1,IOFF+JJ-1)
!      endif
!    endif
!  ENDDO
!ENDDO

!istat = cuda_Synchronize()
DO KGL=IBEG,IEND,IINC

   !IF( T%LUSEFFT992(KGL) )THEN
   !   stop 'Error: code path not (yet) supported in GPU version'
   !END IF
   
  ITYPE=-1
  IJUMP= 1
  IGLG = D_NPTRLS(MYSETW)+KGL-1
  IST  = 2*(G_NMEN(IGLG)+1)+1
  !ILEN = G_NLOEN(IGLG)+R_NNOEXTZL+3-IST

  !IF (G_NLOEN(IGLG)>1) THEN
    IOFF=D_NSTAGTF(KGL)+1
    !IRLEN=G_NLOEN(IGLG)+R_NNOEXTZL
    !ICLEN=(IRLEN/2+1)*2

    CALL CREATE_PLAN_FFT(IPLAN_R2C,-1,KN=G_NLOEN(IGLG),KLOT=KFIELDS)
    !$ACC host_data use_device(ZGTF)
    CALL EXECUTE_PLAN_FFTC(IPLAN_R2C,-1,ZGTF(1, IOFF))
    !$ACC end host_data

   !ENDIF
END DO

#ifdef UNDEF
!$acc update host (ZGTF)
BLOCK
  INTEGER :: ISTR, JFLD, JGL
  REAL (KIND=JPRBT) :: ZZGTF (UBOUND (ZGTF, 1), UBOUND (ZGTF, 2))
! REAL (KIND=JPRBT) :: ZZGTF (1,D%NLENGTF)
  WHERE (ABS (ZGTF) < 1E-15) 
    ZZGTF = 0._JPRB
  ELSEWHERE
    ZZGTF = ZGTF
  ENDWHERE
  WRITE (88, *) __FILE__, ':', __LINE__ 
  WRITE (88, *) "-------------"
  WRITE (88, *) " KFIELDS = ", KFIELDS
  ISTR = D%NLENGTF / D%NDGL_FS
  DO JFLD = 1, 1
  WRITE (88, *) " JFLD = ", JFLD
  DO JGL = 1, D%NDGL_FS
    WRITE (88, '(100E12.4)') ZZGTF (JFLD,1+(JGL-1)*ISTR:JGL*ISTR)
  ENDDO
  ENDDO
  CALL FLUSH (88)
ENDBLOCK
#endif

istat = cuda_Synchronize()

!$ACC data &
!$ACC& COPY(D,D_NSTAGTF,D_NPTRLS,G_NMEN,G_NMEN_MAX,G_NLOEN,G_NLOEN_MAX,R_NNOEXTZL) &
!$ACC& present(ZGTF)

!$ACC parallel loop collapse(3) private(JMAX,KGL,IOFF,SCAL)
DO IGLG=IBEG+OFFSET_VAR-1,IEND+OFFSET_VAR-1,IINC
   DO JJ=1, G_NLOEN_MAX+2
      DO JF=1,KFIELDS
         JMAX = G_NLOEN(IGLG)
         if (JJ .le. JMAX) then
           KGL=IGLG-OFFSET_VAR+1
           IOFF=D_NSTAGTF(KGL)+1
           SCAL = 1._JPRBT/REAL(G_NLOEN(IGLG),JPRBT)
           ZGTF(JF,IOFF+JJ-1)= SCAL * ZGTF(JF, IOFF+JJ-1)
         end if
      ENDDO
   ENDDO
END DO

! debug
!iunit=myproc+300
!!$ACC update host(ZGTF)
!!$ACC wait
!DO IGLG=IBEG+OFFSET_VAR-1,IEND+OFFSET_VAR-1,IINC
!  DO JJ=1, G_NLOEN_MAX+2
!      DO JF=1,KFIELDS
!!
!    JMAX = G_NLOEN(IGLG)+2
!    if (JJ .le. JMAX) then
!      IF (G_NLOEN(IGLG)>1) THEN
!      KGL=IGLG-OFFSET_VAR+1
!      IOFF=D_NSTAGTF(KGL)+1
!      !if( JF.eq.(41+137-1) .and. JJ .LT. 10 ) write(iunit,*) 'fft_01 ',JF,IGLG,KGL,IOFF+JJ-1,G_NLOEN(IGLG),G_NMEN(IGLG),ZGTF(JF,IOFF+JJ-1)
!      endif
!    endif
!  ENDDO
!  ENDDO
!ENDDO

!$ACC parallel loop !collapse(3) private(IGLG,IST,ILEN,IST1)
DO KGL=IBEG,IEND,IINC

   IGLG = D_NPTRLS(MYSETW)+KGL-1
   IST  = 2*(G_NMEN(IGLG)+1)+1
   ILEN = G_NLOEN(IGLG)+R_NNOEXTZL+3-IST
   
   IST1=1
   IF (G_NLOEN(IGLG)==1) IST1=0

   !$ACC loop collapse(2)
   DO JJ=IST1, ILEN
      DO JF=1,KFIELDS
               ZGTF(JF,IST+D_NSTAGTF(KGL)+JJ-1) = 0.0_JPRBT
      ENDDO
   ENDDO
END DO

!$ACC end data

#ifdef UNDEF
!$acc update host (ZGTF)
BLOCK
  INTEGER :: ISTR, JFLD, JGL
  REAL (KIND=JPRBT) :: ZZGTF (UBOUND (ZGTF, 1), UBOUND (ZGTF, 2))
! REAL (KIND=JPRBT) :: ZZGTF (1,D%NLENGTF)
  WHERE (ABS (ZGTF) < 1E-15) 
    ZZGTF = 0._JPRB
  ELSEWHERE
    ZZGTF = ZGTF
  ENDWHERE
  WRITE (88, *) __FILE__, ':', __LINE__ 
  WRITE (88, *) "-------------"
  WRITE (88, *) " KFIELDS = ", KFIELDS
  ISTR = D%NLENGTF / D%NDGL_FS
  DO JFLD = 1, 1
  WRITE (88, *) " JFLD = ", JFLD
  DO JGL = 1, D%NDGL_FS
    WRITE (88, '(100E12.4)') ZZGTF (JFLD,1+(JGL-1)*ISTR:JGL*ISTR)
  ENDDO
  ENDDO
  CALL FLUSH (88)
ENDBLOCK
#endif

!     ------------------------------------------------------------------

END SUBROUTINE FTDIR
END MODULE FTDIR_MOD
