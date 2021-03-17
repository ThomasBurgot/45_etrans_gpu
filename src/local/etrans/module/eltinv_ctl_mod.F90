MODULE ELTINV_CTL_MOD
CONTAINS
SUBROUTINE ELTINV_CTL(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,&
 & PSPVOR,PSPDIV,PSPSCALAR,&
 & PSPSC3A,PSPSC3B,PSPSC2,&
 & KFLDPTRUV,KFLDPTRSC,PSPMEANU,PSPMEANV,FSPGL_PROC)

!**** *ELTINV_CTL* - Control routine for inverse Legandre transform.

!     Purpose.
!     --------
!        Control routine for the inverse LEGENDRE transform

!**   Interface.
!     ----------
!     CALL EINV_TRANS_CTL(...)
!     KF_OUT_LT    - number of fields coming out from inverse LT
!     KF_UV        - local number of spectral u-v fields
!     KF_SCALARS   - local number of scalar spectral fields
!     KF_SCDERS    - local number of derivatives of scalar spectral fields
!     PSPVOR(:,:)  - spectral vorticity (input)
!     PSPDIV(:,:)  - spectral divergence (input)
!     PSPSCALAR(:,:) - spectral scalarvalued fields (input)
!     KFLDPTRUV(:) - field pointer array for vor./div.
!     KFLDPTRSC(:) - field pointer array for PSPSCALAR
!     FSPGL_PROC  - external procedure to be executed in fourier space
!                   before transposition

!     Method.
!     -------

!     Externals.
!     ----------

!     Author.
!     -------
!        Mats Hamrud *ECMWF*

!     Modifications.
!     --------------
!        Original : 00-06-03
!        M.Hamrud      01-Oct-2003 CY28 Cleaning
!        O.Spaniel     Oct-2004 phasing for AL29
!     ------------------------------------------------------------------

USE PARKIND1  ,ONLY : JPIM     ,JPRB
USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

USE TPM_GEN         ,ONLY : LALLOPERM
!USE TPM_DIM
USE TPM_TRANS       ,ONLY : FOUBUF, FOUBUF_IN
USE TPM_DISTR       ,ONLY : D

USE ELTINV_MOD      ,ONLY : ELTINV
USE TRMTOL_MOD      ,ONLY : TRMTOL
!

IMPLICIT NONE

INTEGER(KIND=JPIM),INTENT(IN) :: KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPVOR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPDIV(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSCALAR(:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSC3A(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSC3B(:,:,:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN)  :: PSPSC2(:,:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRUV(:)
INTEGER(KIND=JPIM),OPTIONAL,INTENT(IN) :: KFLDPTRSC(:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN) :: PSPMEANU(:)
REAL(KIND=JPRB) ,OPTIONAL, INTENT(IN) :: PSPMEANV(:)

EXTERNAL  FSPGL_PROC
OPTIONAL  FSPGL_PROC

INTEGER(KIND=JPIM) :: JM,IM,IBLEN,ILEI2,IDIM1
REAL(KIND=JPRB) :: ZHOOK_HANDLE

!     ------------------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('ELTINV_CTL_MOD:ELTINV_CTL',0,ZHOOK_HANDLE)
ILEI2 = 8*KF_UV + 2*KF_SCALARS + 2*KF_SCDERS
IDIM1 = 2*KF_OUT_LT
IBLEN = D%NLENGT0B*2*KF_OUT_LT

IF(KF_OUT_LT > 0) THEN
CALL GSTATS(1647,0)
    !!$ACC DATA &
    !!$ACC      COPY   (ZIA)                    &
    !!$ACC      COPY   (ZEPSNM)                 &
    !!$ACC      COPYOUT(ZAOA1,ZSOA1) &
    !!$ACC COPYIN(D,D%NASM0,D%NUMP,D%MYMS,D%NPROCL,D%NPMT,D%NPNTGTB1,D%NSTAGT0B)  &
    !!$ACC COPYIN(R,F,G,F%REPSNM,D%NPNTGTB1,D%NSTAGT0B,R%NDGNH,G%NDGLU)           &
    !!$ACC  COPYIN(KFLDPTRSC)

    !!$ACC DATA IF( PRESENT(PSPVOR)) COPYIN(PSPVOR)
    !!$ACC DATA IF( PRESENT(PSPDIV)) COPYIN(PSPDIV)
    !!$ACC DATA IF( PRESENT(PSPSCALAR))COPYIN(PSPSCALAR)
    !!$ACC DATA IF( PRESENT(PSPSC3A)) COPYIN(PSPSC3A)
    !!$ACC DATA IF( PRESENT(PSPSC3B)) COPYIN(PSPSC3B)
    !!$ACC DATA IF( PRESENT(PSPSC2)) COPYIN(PSPSC2)

      CALL ABOR1 ('BROKEN ELTINV')
!     CALL ELTINV(KF_OUT_LT,KF_UV,KF_SCALARS,KF_SCDERS,ILEI2,IDIM1,&
!         & PSPVOR,PSPDIV,PSPSCALAR ,&
!         & PSPSC3A,PSPSC3B,PSPSC2 , &
!         & KFLDPTRUV,KFLDPTRSC,FSPGL_PROC)
    
    !!$ACC end data
    !!$ACC end data
    !!$ACC end data
    !!$ACC end data
    !!$ACC end data
    !!$ACC end data
  
CALL GSTATS(1647,1)
ENDIF

CALL GSTATS(152,0)
CALL TRMTOL(FOUBUF_IN,FOUBUF,2*KF_OUT_LT)
CALL GSTATS(152,1)
IF (LHOOK) CALL DR_HOOK('ELTINV_CTL_MOD:ELTINV_CTL',1,ZHOOK_HANDLE)

!     ------------------------------------------------------------------

END SUBROUTINE ELTINV_CTL
END MODULE ELTINV_CTL_MOD
