MODULE EPRFI2B_MOD
CONTAINS
SUBROUTINE EPRFI2B(KFIELD,PFFA)

!**** *EPRFI2B* - Prepare input work arrays for direct transform

!     Purpose.
!     --------
!        To extract the Fourier fields for a specific zonal wavenumber
!        and put them in an order suitable for the direct Legendre
!        tranforms, i.e. split into symmetric and anti-symmetric part.

!**   Interface.
!     ----------
!     *CALL* *EPRFI2B(..)

!        Explicit arguments :
!        -------------------   KFIELD - number of fields
!                              KM - zonal wavenumber
!                              KMLOC - local zonal wavenumber
!                              PAOA - antisymmetric part of Fourier
!                              fields for zonal wavenumber KM
!                              PSOA - symmetric part of Fourier
!                              fields for zonal wavenumber KM

!        Implicit arguments :  FOUBUF in TPM_TRANS
!        --------------------

!     Method.
!     -------

!     Externals.   None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 90-07-01
!        MPP Group: 95-10-01 Support for Distributed Memory version
!        Modified : 04/06/99 D.Salmond : change order of AIA and SIA
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB       ,JPRBT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE TPM_DIM         ,ONLY : R, R_NDGNH, R_NDGL
USE TPM_TRANS       ,ONLY : FOUBUF
USE TPM_GEOMETRY    ,ONLY : G, G_NDGLU
USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS,D_NSTAGT1B,D_NPROCL,D_NPNTGTB1,MYPROC
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN)  :: KFIELD
INTEGER(KIND=JPIM) :: KM, KMLOC
REAL(KIND=JPRBT)  , INTENT(OUT) :: PFFA(:,:,:)

INTEGER(KIND=JPIM) ::  ISTAN, JF, JGL
INTEGER(KIND=JPIM) :: IJR,IJI
REAL(KIND=JPRB) :: ZHOOK_HANDLE

INTEGER :: J, N

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EPRFI2B_MOD:EPRFI2B',0,ZHOOK_HANDLE)

!*       1.    EXTRACT SYM./ANTISYM. FIELDS FROM FOURIER ARRAY.
!              ------------------------------------------------

#ifdef UNDEF
WRITE (0, *) __FILE__, ':', __LINE__

N = SIZE (FOUBUF)

!$acc parallel num_gangs(1) num_workers(1) vector_length(1) present (FOUBUF) copyin (N)
DO J = 1, N
  PRINT *, J, FOUBUF (J)
ENDDO
!$acc end parallel

WRITE (0, *) __FILE__, ':', __LINE__
!$acc parallel num_gangs(1) num_workers(1) vector_length(1) present (FOUBUF) COPY(R_NDGL,D_NSTAGT1B,D_NPNTGTB1,D_NPROCL,D_NUMP,D_MYMS,G_NDGLU)
PRINT *, "KMLOC", "JGL", "JF", "FOUBUF(ISTAN+IJR)", "FOUBUF(ISTAN+IJI)"
DO KMLOC = 1, D_NUMP
  DO JGL=1,R_NDGL
    DO JF =1,KFIELD
      KM = D_MYMS(KMLOC)
      IJR = 2*(JF-1)+1
      IJI = IJR+1
      ISTAN = (D_NSTAGT1B(D_NPROCL(JGL) )+D_NPNTGTB1(KMLOC,JGL ))*2*KFIELD
      PRINT *, KMLOC, JGL, JF, FOUBUF(ISTAN+IJR), FOUBUF(ISTAN+IJI)
!     PFFA(JGL,IJR,KMLOC) = FOUBUF(ISTAN+IJR)
!     PFFA(JGL,IJI,KMLOC) = FOUBUF(ISTAN+IJI)
    ENDDO
  ENDDO
ENDDO
!$acc end parallel

#endif

!$ACC data &
!$ACC& present(PFFA) &
!$ACC& present(FOUBUF) &
!$ACC& COPY(R_NDGL,D_NSTAGT1B,D_NPNTGTB1,D_NPROCL,D_NUMP,D_MYMS,G_NDGLU)
  
!loop over wavenumber
!$ACC parallel loop collapse(3) private(ISTAN,KM,IJR,IJI)
DO KMLOC = 1, D_NUMP
  DO JGL=1,R_NDGL
    DO JF =1,KFIELD
      KM = D_MYMS(KMLOC)
      IJR = 2*(JF-1)+1
      IJI = IJR+1
      ISTAN = (D_NSTAGT1B(D_NPROCL(JGL) )+D_NPNTGTB1(KMLOC,JGL ))*2*KFIELD
      PFFA(JGL,IJR,KMLOC) = FOUBUF(ISTAN+IJR)
      PFFA(JGL,IJI,KMLOC) = FOUBUF(ISTAN+IJI)
    ENDDO
  ENDDO
ENDDO

!$ACC end data


#ifdef UNDEF
WRITE (0, *) __FILE__, ':', __LINE__
!$acc parallel num_gangs(1) num_workers(1) vector_length(1) present (PFFA) COPY(R_NDGL,D_NSTAGT1B,D_NPNTGTB1,D_NPROCL,D_NUMP,D_MYMS,G_NDGLU)
DO KMLOC = 1, D_NUMP
  DO JGL=1,R_NDGL
    DO JF =1,KFIELD
      IJR = 2*(JF-1)+1
      IJI = IJR+1
      PRINT *, KMLOC, JGL, JF, PFFA(JGL,IJR,KMLOC), PFFA(JGL,IJI,KMLOC)
    ENDDO
  ENDDO
ENDDO
!$acc end parallel
#endif


IF (LHOOK) CALL DR_HOOK('EPRFI2B_MOD:EPRFI2B',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE EPRFI2B
END MODULE EPRFI2B_MOD
