!----------------------------------------------------------------------
!----------------------------------------------------------------------
!   NextGen : Next-Generation Neural Mass Model for a cortical neuron 
!	with dynamics synapses.
!       Here we analyze the homogeneous component of properly normalized
!       network, which boils down to study a self-coupled NextGen model.   
!----------------------------------------------------------------------
!----------------------------------------------------------------------

      SUBROUTINE FUNC(NDIM,U,ICP,PAR,IJAC,F,DFDU,DFDP)
!     ---------- ----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM, ICP(*), IJAC
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM), PAR(*)
      DOUBLE PRECISION, INTENT(OUT) :: F(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: DFDU(NDIM,NDIM), DFDP(NDIM,*)


      DOUBLE PRECISION r_e,v_e,s_e,r_i,v_i,s_i,tau_e,tau_i,tau_se
      DOUBLE PRECISION tau_si,nu_e,nu_i,Delta_e,Delta_i,Jee,Jii,Jei,Jie,Iext_e,Iext_i,pi,eps
      INTEGER i,j

       pi = 4.D0 * DATAN(1.D0)
       ! print*, "The value of pi is", pi
       r_e=U(1)
       v_e=U(2)
       s_e=U(3)
       r_i=U(4)
       v_i=U(5)
       s_i=U(6)
        
       tau_e = 8
       tau_i = 8
       tau_se=1
       tau_si=5
       nu_e = -5
       nu_i = -5
       Delta_e = 1
       Delta_i = 1
       Jee = 5
       Jei = 13
       Jii = 5
       Jie = 13
       Iext_i = 0

       Iext_e=PAR(1) 
       eps=PAR(2)


       F(1)= (Delta_e/(tau_e*pi)+2*r_e*v_e)/tau_e
       F(2)= (v_e**2+nu_e-(pi*r_e*tau_e)**2+Iext_e+tau_e*(Jee+eps)*s_e-tau_e*Jei*s_i)/tau_e
       F(3)= (-s_e+r_e)/tau_se
       F(4)= (Delta_i/(tau_i*pi)+2*r_i*v_i)/tau_i
       F(5)= (v_i**2+nu_i-(pi*r_i*tau_i)**2+Iext_i+tau_i*(Jie+eps)*s_e-tau_i*Jii*s_i)/tau_i
       F(6)= (-s_i+r_i)/tau_si


      IF(IJAC.EQ.0)RETURN

        do i = 1, 6
                do j = 1, 6
                        DFDU(i,j)=0;
                end do
        end do

        DFDU(1,1)=(2*v_e)/tau_e;
        DFDU(2,1)=-2*r_e*tau_e*(pi)**2;
        DFDU(3,1)=1/tau_se;

        DFDU(1,2)=(2*r_e)/tau_e;
        DFDU(2,2)=(2*v_e)/tau_e;

        DFDU(2,3)=Jee+eps;
        DFDU(3,3)=-1/tau_se;
        DFDU(5,3)=Jie+eps;

        DFDU(4,4)=2*v_i/tau_i;
        DFDU(5,4)=-2*r_i*tau_i*(pi)**2;
        DFDU(6,4)= 1/tau_si;

        DFDU(4,5)=(2*r_i)/tau_i;
        DFDU(5,5)=(2*v_i)/tau_i;

        DFDU(2,6)=-Jei;
        DFDU(5,6)=-Jii;
        DFDU(6,6)=-1/tau_si;


      IF(IJAC.EQ.1)RETURN

        do i = 1, 6
                do j = 1, 2
                        DFDP(i,j)=0;
                end do
        end do

       DFDP(2,1)=1/tau_e
       DFDP(2,2)=s_e
       DFDP(5,2)=s_e

      END SUBROUTINE FUNC

      SUBROUTINE STPNT(NDIM,U,PAR,T)
!     ---------- -----

      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(INOUT) :: U(NDIM),PAR(*)
      DOUBLE PRECISION, INTENT(IN) :: T

      PAR(1)=0.0 ! Iext_e
      PAR(2)=0.0 ! eps

      !U(1)=0.05708775 !0.05145963
      !U(2)=-0.34848748 !-0.38660146
      !U(3)=0.05708775 !0.05145963
      !U(4)=-0.05708775 !0.03343058
      !U(5)=0.34848748 !-0.59509495
      !U(6)=-0.05708775 !0.03343058

      U(1)=0.0080890471116
      U(2)=-2.4594204499
      U(3)=0.0080890471116
      U(4)=0.0096867123919
      U(5)=-2.0537791442
      U(6)=0.0096867123919


      END SUBROUTINE STPNT

      SUBROUTINE BCND
      END SUBROUTINE BCND

      SUBROUTINE ICND
      END SUBROUTINE ICND

      SUBROUTINE FOPT
      END SUBROUTINE FOPT


!-----------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETUY_MAX(U,NDX,NTST,NCOL)
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION  :: AUX_Y(0:NTST*NCOL)

      AUX_Y=U(2,:)

      GETUY_MAX=MAXVAL(AUX_Y)

      END FUNCTION GETUY_MAX
!-----------------------------------------------------------
!-----------------------------------------------------------
      DOUBLE PRECISION FUNCTION GETUY_MIN(U,NDX,NTST,NCOL)
      INTEGER, INTENT(IN) :: NDX,NCOL,NTST
      DOUBLE PRECISION, INTENT(IN) :: U(NDX,0:NCOL*NTST)
      DOUBLE PRECISION  :: AUX_Y(0:NTST*NCOL)

      AUX_Y=U(2,:)

      GETUY_MIN=MINVAL(AUX_Y)

      END FUNCTION GETUY_MIN
!-----------------------------------------------------------

      SUBROUTINE PVLS(NDIM,U,PAR)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NDIM
      DOUBLE PRECISION, INTENT(IN) :: U(NDIM)
      DOUBLE PRECISION, INTENT(INOUT) :: PAR(*)
      DOUBLE PRECISION, EXTERNAL :: GETP,GETUY_MAX,GETUY_MIN
      INTEGER NDX,NCOL,NTST

      NDX=NINT(GETP('NDX',0,U))
      NTST=NINT(GETP('NTST',0,U))
      NCOL=NINT(GETP('NCOL',0,U))

      PAR(3)=GETUY_MIN(U,NDX,NTST,NCOL)
      PAR(4)=GETUY_MAX(U,NDX,NTST,NCOL)

      END SUBROUTINE PVLS
