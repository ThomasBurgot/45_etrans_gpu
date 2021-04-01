MODULE EUVTVD_MOD
CONTAINS
SUBROUTINE EUVTVD(KFIELD)
!SUBROUTINE EUVTVD(KM,KMLOC,KFIELD,KFLDPTR,PU,PV,PVOR,PDIV,PSPMEANU,PSPMEANV)

!**** *EUVTVD* - Compute vor/div from u and v in spectral space

!     Purpose.
!     --------
!        To compute vorticity and divergence from u and v in spectral
!       space. Input u and v from KM to NTMAX+1, output vorticity and
!       divergence from KM to NTMAX - calculation part.

!**   Interface.
!     ----------
!        CALL EUVTVD(KM,KFIELD,PEPSNM,PU,PV,PVOR,PDIV)

!        Explicit arguments :  KM - zonal wave-number
!        --------------------  KFIELD - number of fields (levels)
!                              KFLDPTR - fields pointers
!                              PEPSNM - REPSNM for wavenumber KM
!                              PU - u wind component for zonal
!                                   wavenumber KM
!                              PV - v wind component for zonal
!                                   wavenumber KM
!                              PVOR - vorticity for zonal
!                                     wavenumber KM
!                              PDIV - divergence for zonal
!                                     wavenumber KM

!     Method.  See ref.
!     -------

!     Externals.  None.
!     ----------

!     Reference.
!     ----------
!        ECMWF Research Department documentation of the IFS

!     Author.
!     -------
!        Mats Hamrud and Philippe Courtier  *ECMWF*

!     Modifications.
!     --------------
!        Original : 91-07-01
!        D. Giard : NTMAX instead of NSMAX
!        01-08-27 : R. El Khatib Fix for NPROMATR /= 0
!        03-03-03 : G. Radnoti: b-level conform mean-wind distribution
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        F. Vana + NEC 28-Apr-2009 MPI-OpenMP fix
!        D. Degrauwe  (Feb 2012): Alternative extension zone (E')
!        N. Lopes & R. El Khatib 15-Jun-2012 Scalability enhancement
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRBT
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE TPM_DIM         ,ONLY : R
!USE TPM_FIELDS
USE TPM_DISTR       ,ONLY : D, NPRCIDS, NPRTRW, MYSETV, MYSETW, MYPROC, NPROC
USE TPMALD_GEO      ,ONLY : GALD
USE TPMALD_DISTR    ,ONLY : DALD
USE MPL_MODULE  ,ONLY : MPL_RECV, MPL_SEND
USE SET2PE_MOD      ,ONLY : SET2PE
USE TPM_FIELDS      ,ONLY : ZOA1,ZOA2,ZEPSNM
USE TPM_DIM         ,ONLY : R, R_NTMAX

USE TPM_DISTR       ,ONLY : D,D_NUMP,D_MYMS



USE ABORT_TRANS_MOD ,ONLY : ABORT_TRANS
!

IMPLICIT NONE

INTEGER(KIND=JPIM), INTENT(IN) :: KFIELD 
INTEGER(KIND=JPIM) :: KM, KMLOC

!REAL(KIND=JPRB), INTENT(INOUT) :: PU  (:,:),PV  (:,:)
!REAL(KIND=JPRB), INTENT(OUT)   :: PVOR(:,:),PDIV(:,:)

!INTEGER(KIND=JPIM), OPTIONAL, INTENT(IN)  :: KFLDPTR(:)
!REAL(KIND=JPRB),    OPTIONAL, INTENT(IN) :: PSPMEANU(:),PSPMEANV(:)

INTEGER(KIND=JPIM) :: II, IN, IR, J, JN, JNMAX
INTEGER(KIND=JPIM) :: IUS, IUE, IVS, IVE, IVORS, IVORE, IDIVS, IDIVE


REAL(KIND=JPRBT), POINTER :: PU(:,:,:),PV(:,:,:),PVOR(:,:,:),PDIV(:,:,:)


REAL(KIND=JPRBT) :: ZKM
REAL(KIND=JPRBT) :: ZSPU(2*KFIELD)
REAL(KIND=JPRBT) :: ZIN
INTEGER(KIND=JPIM) :: JA,ITAG,ILEN,IFLD,ISND
REAL(KIND=JPRBT) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('EUVTVD_MOD:EUVTVD',0,ZHOOK_HANDLE)

IUS = 1
IUE = 2*KFIELD
IVS = 2*KFIELD+1
IVE = 4*KFIELD
IVORS = 1
IVORE = 2*KFIELD
IDIVS = 2*KFIELD+1
IDIVE = 4*KFIELD


!$acc data&
!$acc& copy(D_MYMS,D_NUMP,R_NTMAX) &
!$acc& present(ZOA1,ZOA2)




!*       1.    COMPUTE U V FROM VORTICITY AND DIVERGENCE.
!              ------------------------------------------

PU => ZOA1(IUS:IUE,:,:)
PV => ZOA1(IVS:IVE,:,:)
PVOR => ZOA2(IVORS:IVORE,:,:)
PDIV => ZOA2(IDIVS:IDIVE,:,:)

!$ACC parallel loop collapse(3) private(IR,II,KM,ZKM)
DO KMLOC=1,D_NUMP
  DO JN=1,R%NDGL+R%NNOEXTZG
     DO J=1,KFIELD
       IR=2*J-1
       II=IR+1
       KM = D_MYMS(KMLOC)
       ZKM=REAL(KM,JPRBT)*GALD%EXWN
       PDIV(IR,JN,KMLOC)=-ZKM*PU(II,JN,KMLOC)
       PDIV(II,JN,KMLOC)=ZKM*PU(IR,JN,KMLOC)
       PVOR(IR,JN,KMLOC)=-ZKM*PV(II,JN,KMLOC)
       PVOR(II,JN,KMLOC)=ZKM*PV(IR,JN,KMLOC)
     ENDDO
  ENDDO
ENDDO




!ZKM=REAL(KM,JPRB)*GALD%EXWN
!DO J=1,KFIELD
!  IR=2*J-1
!  II=IR+1
!  DO JN=1,R%NDGL+R%NNOEXTZG
!    PDIV(JN,IR)=-ZKM*PU(JN,II)
!    PDIV(JN,II)= ZKM*PU(JN,IR)
!    PVOR(JN,IR)=-ZKM*PV(JN,II)
!    PVOR(JN,II)= ZKM*PV(JN,IR)
!  ENDDO
!ENDDO



JNMAX = MAXVAL(DALD%NCPL2M)


!$ACC parallel loop collapse(3) private(IN,KM,ZIN)
DO KMLOC=1,D_NUMP
  DO JN=1,JNMAX,2
    DO J=1,2*KFIELD 
      KM = D_MYMS(KMLOC)
      IF (JN <= DALD%NCPL2M(KM)) THEN
        IN=(JN-1)/2
        ZIN=REAL(IN,JPRBT)*GALD%EYWN
        PVOR(J,JN,KMLOC )=PVOR(J,JN,KMLOC)+ZIN*PU(J,JN+1,KMLOC)
        PVOR(J,JN+1,KMLOC)=PVOR(J,JN+1,KMLOC)-ZIN*PU(J,JN,KMLOC)
        PDIV(J,JN,KMLOC )=PDIV(J,JN,KMLOC)-ZIN*PV(J,JN+1,KMLOC)
        PDIV(J,JN+1,KMLOC)=PDIV(J,JN+1,KMLOC)+ZIN*PV(J,JN,KMLOC)
      ENDIF
    ENDDO
  ENDDO
ENDDO






!!$ACC parallel loop collapse(3) private(IN,KM,ZIN)
!DO KMLOC=1,D_NUMP
!  KM = D_MYMS(KMLOC)
!  DO JN=1,DALD%NCPL2M(KM),2
!    DO J=1,2*KFIELD 
!      IN=(JN-1)/2
!      ZIN=REAL(IN,JPRBT)*GALD%EYWN
!      PVOR(J,JN,KMLOC )=PVOR(J,JN,KMLOC)+ZIN*PU(J,JN+1,KMLOC)
!      PVOR(J,JN+1,KMLOC)=PVOR(J,JN+1,KMLOC)-ZIN*PU(J,JN,KMLOC)
!      PDIV(J,JN,KMLOC )=PDIV(J,JN,KMLOC)-ZIN*PV(J,JN+1,KMLOC)
!      PDIV(J,JN+1,KMLOC)=PDIV(J,JN+1,KMLOC)+ZIN*PV(J,JN,KMLOC)
!    ENDDO
!  ENDDO
!ENDDO




!DO J=1,2*KFIELD
!  DO JN=1,DALD%NCPL2M(KM),2
!    IN=(JN-1)/2
!    ZIN=REAL(IN,JPRB)*GALD%EYWN
!    PVOR(JN,J  )=PVOR(JN  ,J)+ZIN*PU(JN+1,J)
!    PVOR(JN+1,J)=PVOR(JN+1,J)-ZIN*PU(JN  ,J)
!    PDIV(JN,J  )=PDIV(JN  ,J)-ZIN*PV(JN+1,J)
!    PDIV(JN+1,J)=PDIV(JN+1,J)+ZIN*PV(JN  ,J)
!  ENDDO
!ENDDO

!$acc end data


WRITE(*,*) 'PDIV_THOMAS', PDIV


END SUBROUTINE EUVTVD
END MODULE EUVTVD_MOD
