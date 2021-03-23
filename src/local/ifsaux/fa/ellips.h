! Jan-2011 P. Marguinaud Interface to thread-safe FA
SUBROUTINE _ELLIPS_ (KSMAX,KMSMAX,KNTMP,KMTMP)
USE PARKIND1, ONLY : JPRB
USE YOMHOOK , ONLY : LHOOK, DR_HOOK
USE LFI_PRECISION
IMPLICIT NONE
!
! ***ELLIPS*** - General routine for computing elliptic truncation
!
!    Purpose.
!    --------
!       Computation of zonal and meridional limit wavenumbers within the ellipse
!    Interface:
!    ----------
!                   *CALL* *ELLIPS *
!
!        Explicit arguments :
!        --------------------
!
!        Implicit arguments :
!        --------------------
!
!
!     Method.
!     -------
!        See documentation
!
!     Externals.   NONE.
!     ----------
!
!     Reference.
!     ----------
!        ARPEGE/ALADIN documentation
!
!     Author.
!     -------
!        G. Radnoti LACE 97/04/04
!
!     Modifications.
!-------------------------------------------------------------
!        J.Vivoda, 99/05/19  treating NSMAX=0 and NMSMAX=0
!
!
INTEGER (KIND=JLIK) KSMAX, KMSMAX
INTEGER (KIND=JLIK) KNTMP(0:KMSMAX),KMTMP(0:KSMAX)
!
INTEGER (KIND=JLIK) JM, JN
!
REAL (KIND=JPDBLR) ZEPS, ZKN, ZKM, ZAUXIL
!
REAL(KIND=JPRB) :: ZHOOK_HANDLE
IF (LHOOK) CALL DR_HOOK('ELLIPS',0,ZHOOK_HANDLE)

KNTMP = KSMAX
KMTMP = KMSMAX

#ifdef UNDEF
ZEPS=1.E-10
ZAUXIL=0.
!
! 1. Computing meridional limit wavenumbers along zonal wavenumbers
!
DO JM=1,KMSMAX-1
ZKN = REAL(KSMAX,JPDBLR)/REAL(KMSMAX,JPDBLR)* &
& SQRT(MAX(ZAUXIL,REAL(KMSMAX**2-JM**2,JPDBLR)))
  KNTMP(JM)=INT(ZKN+ZEPS, JLIK)
ENDDO

IF( KMSMAX.EQ.0 )THEN
   KNTMP(0)=KSMAX
ELSE
   KNTMP(0)=KSMAX
   KNTMP(KMSMAX)=0
ENDIF
!
! 2. Computing zonal limit wavenumbers along meridional wavenumbers
!             
DO JN=1,KSMAX-1
 ZKM = REAL(KMSMAX,JPDBLR)/REAL(KSMAX,JPDBLR)* &
     & SQRT(MAX(ZAUXIL,REAL(KSMAX**2-JN**2,JPDBLR)))
  KMTMP(JN)=INT(ZKM+ZEPS, JLIK)
ENDDO   

IF( KSMAX.EQ.0 )THEN
   KMTMP(0)=KMSMAX
ELSE
   KMTMP(0)=KMSMAX
   KMTMP(KSMAX)=0
ENDIF
#endif
!
IF (LHOOK) CALL DR_HOOK('ELLIPS',1,ZHOOK_HANDLE)
END      
