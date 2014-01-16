%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE WHEEL SET PROBLEM: INDEX = 2.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function WHEEL
global MYtagSTIFF LOUT YPRIME

%FOR USE IN RESMBS%


MYtagSTIFF = 1;

ND=23;
NEQN = 23;
N = 23;

LWORK=(41+3*ND)*ND+3;
LIWORK=ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);
YPRIME = zeros(1,ND);


%MASBND = zeros(1,4);
%MBND = zeros(1,4);
            
         MASBND(1) = 1;
         MASBND(2) = 0;
         MASBND(3) = 0;
         MASBND(4) = 1;
         MBND(1) = NEQN;
         MBND(2) = NEQN;
         MBND(3) = NEQN;
         MBND(4) = NEQN;
         

MYVAR = 0;
RPAR = 0;
IPAR = 0;
IERR = 0;

    
% C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% C     WORK() HOUSES THE FOLLOWING ARRAYS
% C
% C     Y(N,12)   , YHOLD(N,12) , YNHOLD(N,2) , YMAX(N)
% C     ERRORS(N) , SAVE1(N) , SAVE2(N) , SCALE(N) , ARH(N) , PW(MBND(4)*N)
% C     PWCOPY(MBND(4)*N) ,AM(MASBND(4)*N)
% C     IF THE SPARSE OPTION IS NOT BEING USED THEN MBND(4)=N=MASBND(4).
% C     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
IWORK_15 = zeros(1,ND);
WORK_1   = 0; %WORK(1)
WORK_2   = 0; %WORK(2)
WORK_3   = zeros(ND,12); %WORK(3).....ETC
WORK_I1  = zeros(ND,12);
WORK_I2  = zeros(ND,2);
WORK_I3  = zeros(1,ND);
WORK_I4  = zeros(1,ND);
WORK_I5  = zeros(1,ND);
WORK_I6  = zeros(1,ND);
WORK_I7  = zeros(1,ND);
WORK_I8  = zeros(1,ND);
WORK_I9  = zeros(1,MBND(4)*ND);
WORK_I10 = zeros(1,MBND(4)*ND);
WORK_I11 = zeros(1,MASBND(4)*ND);
      
LOUT = fopen('wheel.txt','w');
fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE WHEEL PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r Number of Equations : %g \r',ND);


%... ENDPOINT OF INTEGRATION
         X=0;
         XEND=10;
         XOUT=XEND;
         
%...  REQUIRED TOLERANCE

        %RTOL=.316227766016837972E-06;
        RTOL(1:16)=0.100000000000000009E-06;
        ATOL=RTOL;
        RTOL(17:23) = 1E10;
        ATOL(17:23) = 1E10;
        ITOL=5;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E and %E',RTOL(1),RTOL(17));
        fprintf(LOUT,'\r ATOL = %E and %E',ATOL(1),ATOL(17));
       
%...  INITIAL VALUES
        [N,X,Y,YPRIME] = INIT(N,X,Y,YPRIME);
%...  SET DEFAULT VALUES
        MF=22;
        INDEX=1;
        MAXDER=7;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=100000;
        H=1E-8;
        fprintf(LOUT,'Hstart = %E',H);
         IWORK(1)=15;
         IWORK(2)=8;
         IWORK(3)=0  ;    
      

%...
%... CALL OF THE SUBROUTINE
%...



MYtag = 0;
MYcase = 220;
    while (MYtag == 0)
        
        
        switch MYcase
            
            case 220
          
          [N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,WORK_3,...
              WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,WORK_I7,...
              WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,IWORK_15,...
              MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR] =MEBDF...
              (N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,...
              WORK_3,WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,...
              WORK_I7,WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,...
              IWORK_15,MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR);  

           
           if (INDEX == 1)
               INDEX = 0;
               MYcase = 220;
               continue
           end
           
          MYtag = 34;
          continue
                
        end
    end
    
     TRUE = SOLN(1);


 fprintf(LOUT,'\r\rX = %g',X);
fprintf(LOUT,'\r\rXOUT = %g',XOUT);


for t = 1:N
    
         fprintf(LOUT,'\r\rAPROX Y%g = %E',t,Y(t));
         fprintf(LOUT,'\rTRUE  Y%g = %E\r\r',t,TRUE(t));

        
end      

      fclose(LOUT)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function  y = SOLN(dummy)
        
% c     DASSL applied to Wheelset problem, tend = 10
% c
% c     ATOL = RTOL = 1d-9 for p, v, q and
% c     ATOL = RTOL = 1d10 for lambda to exclude the Lagrange multipliers
% c                            from error control
% c
% c     # steps =             58510
% c     # steps accepted =    54229
% c     # f-eval =            95533
% c     # Jac-eval =          10218
% c
      y( 1) =  0.86355386965811D-02;
      y( 2) =  0.13038281022727D-04;
      y( 3) = -0.93635784016818D-04;
      y( 4) = -0.13642299804033D-01;
      y( 5) =  0.15292895005422D-02;
      y( 6) = -0.76985374142666D-01;
      y( 7) = -0.25151106429207D-03;
      y( 8) =  0.20541188079539D-02;
      y( 9) = -0.23904837703692D+00;
      y(10) = -0.13633468454173D-01;
      y(11) = -0.24421377661131D+00;
      y(12) = -0.33666751972196D-03;
      y(13) = -0.15949425684022D+00;
      y(14) =  0.37839614386969D-03;
      y(15) =  0.14173214964613D+00;
% c....
      y(16) =  -7.81798689824689  ;
      y(17) =  -0.2168949688189130D-01;
      y(18) =  0.277742700815839   ;
      y(19) =  12.6070421367451    ;
      y(20) = -1.38452105711844     ; 
      y(21) =  -5.62411409474332 ;
% c....      
      y(22) = -0.10124044903201D-01;
      y(23) = -0.56285630573753D-02;            
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [T,Y,DY,DELTA,IRES,RPAR,IPAR] =...
    RESWHS(T,Y,DY,DELTA,RPAR,IPAR)
% C====================================================================
% C
% C     RESWHS
% C     ======
% C
% C     PART OF THE TEST SET FOR INITIAL VALUE PROBLEMS, SEE
% C
% C        HTTP://WWW.CWI.NL/CWI/PROJECTS/IVPTESTSET.HTML
% C
% C     VERSION :        NOV. 1995
% C     AUTHORS :        SIMEON, FUEHRER, RENTROP
% C     PURPOSE :        RESIDUAL OF EQS. OF MOTION FOR
% C                      SIMULATION OF WHEELSET WITH DASSL
% C     DESCRIPTION:     SEE THE ABOVE MENTIONED TESTSET
% C     REFERENCE:       B. SIMEON, C. FUEHRER, P. RENTROP:
% C                      INTRODUCTION TO DIFFERENTIAL-ALGEBRAIC
% C                      EQUATIONS IN VEHICLE SYSTEM DYNAMICS,
% C                      SURV. MATH. IND. I: 1-37 (1991)
% C     SUBROUTINES:     WHEELP   WHEEL PROFILE (CONE)
% C                      RAILP     RAIL PROFILE  (CIRCLE)
% C                      CREEP   CREEP FORCES (ANALYTICAL)
% C                      CONSTM  CONSTRAINT JACOBIAN ("G"-MATRIX)
% C
% C====================================================================
% C
% C     PARAMETERS  (I INPUT , O OUTPUT)
% C     ==========
% C
% C       T            SIMULATION TIME                    (I)
% C       Y(1:17)      DEPENDENT VARIABLES IN THE ORDER   (I)
% C                    ( P, V, BETA, LAMBDA, Q), CF. DESCRIPTION
% C       DY(1:17)     DERIVATIVES OF Y                   (I)
% C       RPAR         SYSTEM PARAMETERS (DUMMY)
% C       IPAR(1)      SWITCH FOR THE CONSTRAINT          (I)
% C                    (=0: INDEX-2 CONSTR., =1: INDEX-3 CONSTR.)
% C       DELTA(1:17)  RESIDUALS OF EQS. OF MOTION        (O)
% C       IRES         ERROR FLAG: 0 = THINGS WENT FINE   (O)
% C                               -1 = THE INTEGRATION MUST BE
% C                                    TERMINATED DUE TO AN ERROR
% C                                    (TYPICALLY A DERAILMENT OF
% C                                     THE WHEELSET)
% C
% C       THE NOTATION OF INTERNAL VARIABLES IS CONFORMAL TO THE
% C       DESCRIPTION OF THE PROBLEM IN THE TESTSET, EXCEPT FOR
% C       THE DYNAMIC VARIABLE BETA WHICH IS HERE CALLED BETAP.
% C
% C==============================================================
% C
%       DOUBLE PRECISION
%      *           T, Y(1:23), DY(1:23), DELTA(1:23), RPAR
%       INTEGER    IRES, IPAR(*)
%       DOUBLE PRECISION
%      *           XX, YY, ZZ, XXP, YYP, ZZP, XXPP, YYPP, ZZPP,
%      *           TETA, PHI, TETAP, PHIP, TETAPP, PHIPP, BETAP, BETAPP,
%      *           E1, E2, PSIL, PSIR, XRL, XRR, XGL, XGR, FNL, FNR,
%      *           RXL, RXR, DRXL, DRXR, D2RXL, D2RXR, D3RXL, D3RXR,
%      *           GXL, GXR, DGXL, DGXR, D2GXL, D2GXR, D3GXL, D3GXR,
%      *           TL(1:3), TR(1:3), QL(1:5), QR(1:5),
%      *           DELTAL, DELTAR, DETER,
%      *           SIT, COT, SIP, COP, SIA, COA, SISL, COSL, SISR, COSR,
%      *           SIDL,CODL,SIDR,CODR,W1,W2,SE,S0,
%      *           ALPHA, S, MR, G, V, FA, OMEGA, RN0, LA, LI1, LI2,
%      *           MA, HA, FQ1, FQ2, FQ3, LM1, LM2, LM3, TOL, MU ,TG,
%      *           CX,CZ, XL
%       INTEGER    IERR, I
TL = zeros(1,3); TR = zeros(1,3); QL=zeros(1,5);   QR=zeros(1,5);
XX=0; YY=0; ZZ=0; XXP=0; YYP=0; ZZP=0; XXPP=0; YYPP=0; ZZPP=0;
TETA=0; PHI=0; TETAP=0; PHIP=0; TETAPP=0; PHIPP=0; BETAP=0; BETAPP=0;
E1=0; E2=0; PSIL=0; PSIR=0; XRL=0; XRR=0; XGL=0; XGR=0; FNL=0; FNR=0;
RXL=0; RXR=0; DRXL=0; DRXR=0; D2RXL=0; D2RXR=0; D3RXL=0; D3RXR=0;
GXL=0; GXR=0; DGXL=0; DGXR=0; D2GXL=0; D2GXR=0; D3GXL=0; D3GXR=0;
RXL = 0;
TOL   = 0.00000001d0;
MR    = 16.08d0;
G     = 9.81000d0;
V     = 30.000d0;
RN0   = 0.10000d0;
LI1   = 0.0605d0;
LI2   = 0.366d0;
MA    = 0.d0;
HA    = 0.2d0;
%  FRICTION COEFFICIENT
MU    = 0.120000d0;
XL    = 0.19d0;
CX   = 6400.d0;
CZ   = 6400.d0;

% C
% C =====================================================================
% C
% C     COORDINATES
% C........ P
      XX      = Y(1);
      YY      = Y(2);
      ZZ      = Y(3);
      TETA    = Y(4);
      PHI     = Y(5);
% c....... p' = v
      XXP     = Y(6);
      YYP     = Y(7);
      ZZP     = Y(8);
      TETAP   = Y(9);
      PHIP    = Y(10);
% c..... B
      BETAP   = Y(11);
% c..... lambda
% C     SCALE LAGRANGE MULTIPLIERS
      E1      = Y(12)*1.d+4;
      E2      = Y(13)*1.d+4;
%c....... q
      PSIL    = Y(14);
      XRL     = Y(15);
      PSIR    = Y(22);
      XRR     = Y(23);
% C.......  w
      XXPP    = Y(16);
      YYPP    = Y(17);
      ZZPP    = Y(18);
      TETAPP  = Y(19);
      PHIPP   = Y(20);
% c....... m
      BETAPP  = Y(21);
% C
% C    STRAIGHT TRACK, NO ADDITIONAL PROPULSION FORCES
% C
      TG      = 0.D0;
      S0      = 0.D0;
      SE      = 0.D0;
      FA      = 0.D0;
      LA      = 0.D0;
      S       =  0.0D0;
      ALPHA   =  0.0D0;
%C
      OMEGA = V/RN0;
      SIT   = sin(TETA);
      COT   = cos(TETA);
      SIP   = sin(PHI);
      COP   = cos(PHI);
      SIA   = sin(ALPHA);
      COA   = cos(ALPHA);
      SISL  = sin(PSIL);
      COSL  = cos(PSIL);
      SISR  = sin(PSIR);
      COSR  = cos(PSIR);
%C
      IERR  = 0;
      IRES  = 0;
% C
% C    EVALUATION OF PROFILE FUNCTIONS
% C    WHEEL, LEFT SIDE

      [XRL,RXL,DRXL,D2RXL, D3RXL,IERR] =...
          WHEELP(XRL,RXL,DRXL,D2RXL, D3RXL,IERR);
      
      if (IERR~=0)
          %GO TO 999
%^^^^999^^^^*
      if (IERR<0)
         IRES = -1;
      else
         IRES = 0;
      end
      return
 %^^^^999^^^^*
      end
      
      XGL   = XX + XRL*COT*COP + RXL*(COT*SIP*COSL - SIT*SISL);
      
%C    RAIL,  LEFT SIDE

      [XGL,GXL,DGXL,D2GXL,D3GXL,IERR] =...
          RAILP(XGL,GXL,DGXL,D2GXL,D3GXL,IERR);
      
      if (IERR~=0)
          %GO TO 999
 %^^^^999^^^^*
      if (IERR<0)
         IRES = -1;
      else
         IRES = 0;
      end
      return
 %^^^^999^^^^*          
      end
      
%C    WHEEL, RIGHT SIDE
 
      [XRR,RXR,DRXR,D2RXR,D3RXR,IERR] =...
          WHEELP(XRR,RXR,DRXR,D2RXR,D3RXR,IERR);
      
      if (IERR~= 0)
          %GO TO 999
 %^^^^999^^^^*
      if (IERR<0)
         IRES = -1;
      else
         IRES = 0;
      end
      return
 %^^^^999^^^^*          
      end
      
      XGR   = XX + XRR*COT*COP + RXR*(COT*SIP*COSR - SIT*SISR);
      
%C    RAIL,  RIGHT SIDE

      [XGR,GXR,DGXR,D2GXR,D3GXR,IERR] =...
          RAILP(XGR,GXR,DGXR,D2GXR,D3GXR,IERR);
      
      if (IERR~=0)
          %GO TO 999
 %^^^^999^^^^*
      if (IERR<0)
         IRES = -1;
      else
         IRES = 0;
      end
      return
 %^^^^999^^^^*          
      end
      
% C
% C     BUILD UP CONSTRAINT JACOBIAN MATRIX "G"
% C     LEFT CONSTRAINT
     
      [XRL,RXL,DGXL,SIT,COT,SIP,COP,SISL,COSL,QL] =...
          CONSTM (XRL,RXL,DGXL,SIT,COT,SIP,COP,SISL,COSL,QL);

%C     RIGHT CONSTRAINT
      [XRR,RXR,DGXR,SIT,COT,SIP,COP,SISR,COSR,QR] =...
          CONSTM (XRR,RXR,DGXR,SIT,COT,SIP,COP,SISR,COSR,QR);
      
% C
% C     ANGLE OF CONTACT PLANE
% C


      W1     = (DRXL*COP - SIP*COSL)*COT + SISL*SIT;
      W2     = -DRXL*SIP - COSL*COP;
      DELTAL = atan( W1/W2 );
      W1     = (DRXR*COP - SIP*COSR)*COT + SISR*SIT;
      W2     =  DRXR*SIP + COSR*COP;
      DELTAR = atan( W1/W2 );
      SIDL   = sin(DELTAL);
      CODL   = cos(DELTAL);
      SIDR   = sin(DELTAR);
      CODR   = cos(DELTAR);
      
% C
% C     NORMAL FORCES N(P,Q,LAMBDA)
% C

      DETER    = -SIDL*CODR - SIDR*CODL;
      if(abs(DETER)<TOL) 
      
          IERR = - 20;
          
          %GO TO 999
 %^^^^999^^^^*
      if (IERR<0)
         IRES = -1;
      else
         IRES = 0;
      end
      return
 %^^^^999^^^^*          
      end
      
      W1     = QL(1)*E1 + QR(1)*E2;
      W2     = QL(2)*E1 + QR(2)*E2;
      FNL    = ( CODR*W1 - SIDR*W2) / DETER;
      FNR    = (-CODL*W1 - SIDL*W2) / DETER;

% C
% C     CREEPAGE FORCES
% C

      [Y, FNL, FNR, V, S, ALPHA, OMEGA,RXL, RXR, DRXL, DRXR,...
          D2RXL, D2RXR,DGXL, DGXR, D2GXL, D2GXR,DELTAL, DELTAR, MU, TL,...
          TR, IERR ] =...
          CREEP ( Y, FNL, FNR, V, S, ALPHA, OMEGA,RXL, RXR, DRXL, DRXR,...
          D2RXL, D2RXR,DGXL, DGXR, D2GXL, D2GXR,DELTAL, DELTAR, MU, TL,...
          TR, IERR );
      
% C
% C     FORCES OF CHASSIS
% C

      FQ1    = MA*G*( V*V*S/G - tan(ALPHA) ) / COA;
      FQ2    = -MA*G*COA*(V*V*S*tan(ALPHA)/G + 1.0d0);
      FQ3    = -2.0d0*CZ*ZZ;
      LM1    = 0.0d0;
      LM2    = -2.0D0*XL^2*CZ*TETA;
      LM3    = -HA*FQ1;
% C
% C --------------------------------------------------------------------
% C
% C     KINEMATICS
% C
% c..... P' = v

      DELTA(1)  = XXP ;
      DELTA(2)  = YYP ;
      DELTA(3)  = ZZP ;
      DELTA(4)  = TETAP; 
      DELTA(5)  = PHIP;
% c..... v' = w
      DELTA(6)  = XXPP;
      DELTA(7)  = YYPP;
      DELTA(8)  = ZZPP;
      DELTA(9)  = TETAPP;
      DELTA(10) = PHIPP;
%c.... B' = m
      DELTA(11) = BETAPP;
% C
% C     DYNAMICS: NEWTON'S LAW
% C
      DELTA(12) = MR*(-XXPP+ V*V*S*COA*(1.0d0 + (XX*COA-YY*SIA)*S)+...
          2.0D0*V*S*COA*ZZP)+ TL(1) + TR(1) + FQ1 - MR*G*SIA+ QL(1)*E1...
          + QR(1)*E2 - 2.0d0*CX*XX;
      
      DELTA(13) = MR*( -YYPP-V*V*S*SIA*(1.0d0 + (XX*COA-YY*SIA)*S)-...
          2.0D0*V*S*SIA*ZZP)+ TL(2) + TR(2) + FQ2 - MR*G*COA+ QL(2)*E1...
          + QR(2)*E2;
      
      DELTA(14)  = MR*( -ZZPP- 2.0d0*V*S*(XXP*COA-YYP*SIA)+ V*V*S*S*ZZ)+...
          TL(3) + TR(3) + FA + FQ3+ QL(3)*E1 + QR(3)*E2;
      
% C
% C     DYNAMICS: EULER'S LAW
% C
      W1     = -(XRL*SIT+RXL*SISL*COT*COP)*TL(1) - RXL*SISL*SIP*TL(2)+...
          (-XRL*COT+RXL*SISL*SIT*COP)*TL(3);
      
      W2     = -(XRR*SIT+RXR*SISR*COT*COP)*TR(1) - RXR*SISR*SIP*TR(2)+...
          (-XRR*COT+RXR*SISR*SIT*COP)*TR(3);
      
      DELTA(15) =  - LI2*( TETAPP*COP - TETAP*PHIP*SIP +...
V*S*(PHIP*(SIA*COT*COP + COA*SIP)- TETAP*SIA*SIT*SIP))-...
LI1*(OMEGA+BETAP)*(PHIP - V*S*SIT*SIA)- (LI1 - LI2)*(TETAP*SIP...
- V*S*(COT*COP*SIA + SIP*COA))*(PHIP - V*S*SIT*SIA)+ W1 + W2 + COP*LM2...
- COT*SIP*LM1 + SIT*SIP*LM3+ QL(4)*E1 + QR(4)*E2;

      W1        = -(XRL*COT*SIP-RXL*COSL*COT*COP)*TL(1)+...
          (XRL*COP+RXL*COSL*SIP)*TL(2)+ (XRL*SIT*SIP-RXL*...
          COSL*SIT*COP)*TL(3);
      
      W2        = -(XRR*COT*SIP-RXR*COSR*COT*COP)*TR(1)+...
          (XRR*COP+RXR*COSR*SIP)*TR(2)+ (XRR*SIT*SIP-RXR...
          *COSR*SIT*COP)*TR(3);
      
      DELTA(16) = - LI2*(PHIPP - TETAP*V*S*SIA*COT)+ LI1*(OMEGA+BETAP)*...
          (TETAP*COP + V*S*(COT*SIP*SIA-COP*COA))+ (LI1-LI2)*(TETAP*SIP...
          - V*S*(COT*COP*SIA + SIP*COA))*(TETAP*COP + V*S*(COT*SIP*SIA-...
          COP*COA))+ W1 + W2 + LM3 + QL(5)*E1 + QR(5)*E2;
      
      W1     = -RXL*(COSL*SIT+SISL*COT*SIP)*TL(1) + RXL*SISL*COP*TL(2)-...
          RXL*(COSL*COT-SISL*SIT*SIP)*TL(3);
      
      W2     = -RXR*(COSR*SIT+SISR*COT*SIP)*TR(1) + RXR*SISR*COP*TR(2)-...
          RXR*(COSR*COT-SISR*SIT*SIP)*TR(3);
      
      DELTA(17) =  - LI1*(BETAPP + TETAPP*SIP + TETAP*PHIP*COP-...
          V*S*(PHIP*(COA*COP-SIA*COT*SIP)- TETAP*SIA*SIT*COP))+...
          W1 + W2 + SIP*LM2 + LA;
      
% C
% C     CONSTRAINT EQUATIONS
% C
% C     CONTACT CONDITION "G_1(P,Q)"
% C

      for I = 18:19
         DELTA(I) = 0.0D0;
      end
      
      if (IPAR(1)==1)
% C
% C     INDEX-3 FORMULATION
         DELTA(18) = GXL - YY - XRL*SIP + RXL*COP*COSL;
         DELTA(19) = GXR - YY - XRR*SIP + RXR*COP*COSR;
      else
% C
% C     INDEX-2 FORMULATION
         for I = 1:5
            DELTA(18) =  DELTA(18)+QL(I)*Y(5+I);
            DELTA(19) =  DELTA(19)+QR(I)*Y(5+I);
         end
      end
% C
%ADDITIONAL INDEX- 1 EQUATIONS "G_2(P,Q)"
%(NORMAL VECTORS OF CONTACT PLANE ARE PARALLEL, NONINTERSECTION CONDITION)
% C
      DELTA(20) = DGXL*(DRXL*SIP + COP*COSL) + DRXL*COT*COP-...
          COT*SIP*COSL + SIT*SISL;
      
      DELTA(21) = DRXL*SIT*COP - SIT*SIP*COSL - COT*SISL;
      
      DELTA(22) = DGXR*(DRXR*SIP + COP*COSR) + DRXR*COT*COP-...
          COT*SIP*COSR + SIT*SISR;
      
      DELTA(23) = DRXR*SIT*COP - SIT*SIP*COSR - COT*SISR;
      
% C
% C     ERROR HANDLING
% C

 %^^^^999^^^^*
      if (IERR<0)
         IRES = -1;
      else
         IRES = 0;
      end
      return
 %^^^^999^^^^*
       
%==========================================================================
%==========================================================================
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%                               END OF RESWHS                             %
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%==========================================================================
%==========================================================================
end     


 
function [XR,RX,DGX,SIT,COT,SIP,COP,SIPS,COPS,Q] =...
    CONSTM(XR,RX,DGX,SIT,COT,SIP,COP,SIPS,COPS,Q)

%       DOUBLE PRECISION XR, RX, DGX, SIT, COT, SIP, COP, SIPS, COPS, Q(5)
% C=======================================================================
% C
% C     constm
% C     ======
% C
% C     Part of the test set for Initial Value Problems, see
% C
% C        http://www.cwi.nl/cwi/projects/IVPtestset.html
% C
% C     Version :        Nov. 1995
% C     Authors :        SIMEON, FUEHRER, RENTROP
% C     Purpose :        Computation of constraint Jacobian for
% C                      simulation of wheelset with DASSL
% C     Description:     See the above mentioned testset
% C     Reference:       B. SIMEON, C. FUEHRER, P. RENTROP:
% C                      Introduction to differential-algebraic
% C                      equations in vehicle system dynamics,
% C                      Surv. Math. Ind. I: 1-37 (1991)
% C
% C     PARAMETER  (I=input, O=output)
% C
% C     XR                      displacement xi                       (I)
% c     rx                      wheel profile                         (I)
% C     DGX                     derivative of rail profile            (I)
% C     SIT                     SIN(TETA),TETA=Y(4)                   (I)
% C     COT                     COS(TETA)                             (I)
% C     SIP                     SIN(PHI),PHI=Y(5)                     (I)
% C     COP                     COS(PHI)                              (I)
% C     SIPS                    SIN(PSIL/R),PSIL/R=Y(14/16)           (I)
% C     COPS                    COS(PSIL/R)                           (I)
% C     Q(1:5)                  constraint matrix (left/right)        (O)
% C
% C
% C     constraint matrix
% C=======================================================================
% C
      Q(1) = DGX;
      Q(2) = -1.0D0;
      Q(3) = 0.0D0;
      Q(4) = DGX*(RX*(-COT*SIPS - COPS*SIP*SIT) - COP*SIT*XR);
      Q(5) = -(COPS*RX*SIP) - COP*XR + DGX*(COP*COPS*COT*RX -COT*SIP*XR);
%==========================================================================
%==========================================================================
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%                               END OF CONSTM                             %
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%==========================================================================
%==========================================================================
end     




function [Y,FNL,FNR,V,S,ALPHA,OMEGA,...
    RXL, RXR, DRXL, DRXR, D2RXL, D2RXR,...
    DGXL, DGXR, D2GXL, D2GXR,...
    DELTAL, DELTAR, MU, TL, TR, IERR] = CREEP(Y,FNL,FNR,V,S,ALPHA,OMEGA,...
    RXL, RXR, DRXL, DRXR, D2RXL, D2RXR,...
    DGXL, DGXR, D2GXL, D2GXR,...
    DELTAL, DELTAR, MU, TL, TR, IERR )

% C====================================================================
% C
% C     CREEP
% C     =====
% C
% C     PART OF THE TEST SET FOR INITIAL VALUE PROBLEMS, SEE
% C
% C        HTTP://WWW.CWI.NL/CWI/PROJECTS/IVPTESTSET.HTML
% C
% C     VERSION :        NOV. 1995
% C     AUTHORS :        SIMEON, FUEHRER, RENTROP
% C     PURPOSE :        COMPUTATION OF CREEPAGE FORCES FOR
% C                      SIMULATION OF WHEELSET WITH DASSL
% C     DESCRIPTION:     SEE THE ABOVE MENTIONED TESTSET
% C     REFERENCE:       B. SIMEON, C. FUEHRER, P. RENTROP:
% C                      INTRODUCTION TO DIFFERENTIAL-ALGEBRAIC
% C                      EQUATIONS IN VEHICLE SYSTEM DYNAMICS,
% C                      SURV. MATH. IND. I: 1-37 (1991)
% C
% C     APPROXIMATION OF CREEPAGE FORCES DUE TO A. JASCHINSKI (DLR),
% C     IN ORDER TO GUARANTEE SMOOTHNESS OF THE FUNCTIONS.
% C     FOR THIS END SEE:
% C     A. JASCHINSKI, DLR-REPORT DLR-FB 90-06, KÖLN, 1990
% C
% C     FOR THE DESCRIPTION OF THE PARAMETERS SEE RESRAD.
% C
% C ==================================================================
% C
% 
%       DOUBLE PRECISION
%      *           Y(1:23), FNL, FNR, V, S, ALPHA, OMEGA,
%      *           RXL, RXR, DRXL, DRXR, D2RXL, D2RXR,
%      *           DGXL, DGXR, D2GXL, D2GXR,
%      *           DELTAL, DELTAR, MU, TL(1:3), TR(1:3)
%       INTEGER    IERR
%       DOUBLE PRECISION
%      *           CL, CR, RHOG, RHOR, RR, A, B, E, G, SIGMA, GM, PI,
%      *           C11, C22, C23, WVR, WVK1, WVK2, WVK3,WVK4,WABS,
%      *           WNUX, WNUY, WPHIS, TX, TY, XX, YY, BETAP, XRL, XRR,
%      *           ZZ,ZZP, WVK5, WVK6,
%      *           SITL,SITR, COTL,COTR, SIP, COP, SIA, COA, SISL, COSL,
%      *           SISR, COSR,
%      *           SIDL, CODL, SIDR, CODR, WT, SQT,
%      *           XXP, YYP, TETAP, PHIP

E     = 1.3537956D0;
G     = 0.7115218D0;
SIGMA = 0.28D0;
GM    = 7.92D+10;
PI    = 3.1415926D0;
C11   = 4.72772197D0;
C22   = 4.27526987D0;
C23   = 1.97203505D0;

 
      IERR  = 0;
      XX    = Y(1);
      YY    = Y(2);
      ZZ    = Y(3);
 
      SIP   = sin(Y(5));
      COP   = cos(Y(5));
 
      XXP   = Y(6);
      YYP   = Y(7);
      ZZP   = Y(8);
      TETAP = Y(9);
      PHIP  = Y(10);
      BETAP = Y(11);

      SISL  = sin(Y(14));
      COSL  = cos(Y(14));
      XRL   = Y(15);
      SISR  = sin(Y(22));
      COSR  = cos(Y(22));
      XRR   = Y(23);
      SIA   = sin(ALPHA);
      COA   = cos(ALPHA);
      SIDL  = sin(DELTAL);
      CODL  = cos(DELTAL);
      SIDR  = sin(DELTAR);
      CODR  = cos(DELTAR);
      SITL   = sin(Y(4)/CODL);
      COTL   = cos(Y(4)/CODL);
      SITR   = sin(Y(4)/CODR);
      COTR   = cos(Y(4)/CODR);

% C
% C     CONTACT ELLIPSES
% C

      RR    = RXL*sqrt(1.D0 + DRXL*DRXL);
      RHOG  = -D2GXL/(1.D0 + DGXL*DGXL)^1.5d0;
      RHOR  =  D2RXL/(1.D0 + DRXL*DRXL)^1.5d0;
      A     = 0.5D0/RR;
      B     = 0.5D0*(RHOG + RHOR);
      WABS  = abs(FNL)*3.D0;
      CL    = ((WABS*(1.D0-SIGMA)*E)...
          /(2.0D0*PI*(A+B)*GM*sqrt(G)))^(1.0d0/3.0d0);
      
      RR    = RXR*sqrt(1.D0 + DRXR*DRXR);
      RHOG  = -D2GXR/(1.D0 + DGXR*DGXR)^1.5d0;
      RHOR  =  D2RXR/(1.D0 + DRXR*DRXR)^1.5d0;
      A     = 0.5D0/RR;
      B     = 0.5D0*(RHOG + RHOR);
      WABS  = abs(FNR)*3.D0;
      CR    = ((WABS*(1.D0-SIGMA)*E) /...
          (2.0D0*PI*(A+B)*GM*sqrt(G)))^(1.0d0/3.0d0);
      
% C
% C     CREEPAGE LEFT CONTACT POINT
% C
% C     RELATIVE VELOCITY

      WVK1 = -(OMEGA+BETAP)*RXL*(SITL*COSL+COTL*SIP*SISL)...
          + V*S*COA*( RXL*(SITL*SIP*COSL+COTL*SISL)...
          + XRL*SITL*COP   - ZZ             )...
          + XXP - TETAP*(RXL*(SITL*SIP*COSL+COTL*SISL)...
          + XRL*SITL*COP)...
          - PHIP*COTL*(XRL*SIP - RXL*COP*COSL);
      
      WVK2 = (OMEGA+BETAP)*RXL*COP*SISL...
          + V*S*SIA*(ZZ -XRL*SITL*COP...
          - RXL*(SITL*SIP*COSL+COTL*SISL) )...
          + YYP + PHIP*(XRL*COP + RXL*SIP*COSL);
      
      WVK3 = -(OMEGA+BETAP)*RXL*(COTL*COSL-SITL*SIP*SISL) + V + ZZP...
          + V*S * ( XX*COA - YY*SIA...
          + COA*(RXL*(COTL*SIP*COSL-SITL*SISL)+XRL*COTL*COP)...
          + SIA*(RXL*COP*COSL - XRL*SIP)                 )...
          - TETAP*(XRL*COTL*COP + RXL*(COTL*SIP*COSL-SITL*SISL))...
          + PHIP*SITL*(XRL*SIP - RXL*COP*COSL);
      
% C
% C     ROLLING VELOCITY
% C

      WVK4 = WVK1 - 2.D0*XXP + 2.D0*V*S*ZZ*COA;
      WVK5 = WVK2 - 2.D0*YYP - 2.D0*V*S*ZZ*SIA;
      WVK6 = WVK3 - 2.D0*ZZP - 2.D0*V*S*(XX*COA-YY*SIA) - 2.D0*V;
      WVR  = 0.5d0*sqrt( WVK4*WVK4 + WVK5*WVK5 + WVK6*WVK6 );
% C
% C     CREEPAGE
% C

      WNUX = ( SITL*WVK1 + COTL*WVK3 ) / WVR;
      WNUY = ( COTL*CODL*WVK1 + SIDL*WVK2 - SITL*CODL*WVK3 ) / WVR;
      WPHIS= (-SIDL*( OMEGA+BETAP - V*S*SIA )+CODL*( TETAP   - V*S*COA ) ) / WVR;
      
% C
% C     CREEPAGE FORCES
% C

      WT   =  MU*FNL;
      TX   = -WT*tanh( GM*CL*CL*C11*WNUX/WT );
      TY   = -WT*tanh( GM*CL*CL*C22*WNUY/WT + GM*CL*CL*CL*C23*WPHIS/WT);
      SQT  =  sqrt ( TX*TX + TY*TY );
      
      if ( SQT*SQT>WT*WT )
% C         NORMALIZE
          TX   = TX * abs(WT) / SQT;
          TY   = TY * abs(WT) / SQT;
          IERR = 3;
      end
      
% C
% C     TRANSFORMATION TO NOMINAL REFERENCE FRAME
% C

      TL(1) = SITL*TX + COTL*CODL*TY;
      TL(2) = +SIDL*TY;
      TL(3) = COTL*TX - SITL*CODL*TY;
% C
% C     CREEPAGE RIGHT CONTACT POINT
% C
% C     RELATIVE VELOCITY

      WVK1 = -(OMEGA+BETAP)*RXR*(SITR*COSR+COTR*SIP*SISR)...
          + V*S*COA*( RXR*(SITR*SIP*COSR+COTR*SISR)...
          + XRR*SITR*COP -ZZ               )...
          + XXP - TETAP*(RXR*(SITR*SIP*COSR+COTR*SISR)...
          + XRR*SITR*COP)...
          - PHIP*COTR*(XRR*SIP - RXR*COP*COSR);
      
      WVK2 = (OMEGA+BETAP)*RXR*COP*SISR...
          + V*S*SIA*(ZZ -XRR*SITR*COP...
          - RXR*(SITR*SIP*COSR+COTR*SISR) )...
          + YYP + PHIP*(XRR*COP + RXR*SIP*COSR);
      
      WVK3 = -(OMEGA+BETAP)*RXR*(COTR*COSR-SITR*SIP*SISR) + V +ZZP...
          + V*S * ( XX*COA - YY*SIA...
          + COA*(RXR*(COTR*SIP*COSR-SITR*SISR)+XRR*COTR*COP)...
          + SIA*(RXR*COP*COSR - XRR*SIP)                 )...
          - TETAP*(XRR*COTR*COP + RXR*(COTR*SIP*COSR-SITR*SISR))...
          + PHIP*SITR*(XRR*SIP - RXR*COP*COSR);
      
% C
% C     ROLLING VELOCITY
% C

      WVK4 = WVK1 - 2.D0*XXP + 2.D0*V*S*ZZ*COA;
      WVK5 = WVK2 - 2.D0*YYP - 2.D0*V*S*ZZ*SIA;
      WVK6 = WVK3 - 2.D0*ZZP - 2.D0*V*S*(XX*COA-YY*SIA) - 2.D0*V;
      WVR  = 0.5d0*sqrt( WVK4*WVK4 + WVK5*WVK5 + WVK6*WVK6 );
      
% C
% C     CREEPAGE
% C

      WNUX = ( SITR*WVK1 + COTR*WVK3 ) / WVR;
      WNUY = ( COTR*CODR*WVK1 - SIDR*WVK2 - SITR*CODR*WVK3 ) / WVR;
      WPHIS= (+SIDR*( OMEGA+BETAP - V*S*SIA )...
          +CODR*( TETAP   - V*S*COA ) ) / WVR;
% C
% C     CREEPAGE FORCES
% C
      WT   =  MU*FNR;
      TX   = -WT*tanh( GM*CR*CR*C11*WNUX/WT );
      TY   = -WT*tanh( GM*CR*CR*C22*WNUY/WT + GM*CR*CR*CR*C23*WPHIS/WT);
      SQT  =  sqrt ( TX*TX + TY*TY );
      
      if ( SQT*SQT>WT*WT )
% C         NORMALIZE
          TX   = TX * abs(WT) / SQT;
          TY   = TY * abs(WT) / SQT;
          IERR = 4;
      end
      
% C
% C     TRANSFORMATION TO NOMINAL REFERENCE FRAME
% C
      TR(1) = SITR*TX + COTR*CODR*TY;
      TR(2) = -SIDR*TY;
      TR(3) = COTR*TX - SITR*CODR*TY;
%==========================================================================
%==========================================================================
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%                               END OF CREEP                              %
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%==========================================================================
%==========================================================================
end     



function [X,SX,DSX,D2SX,D3SX,IERR] = RAILP (X,SX,DSX,D2SX,D3SX,IERR)

% C====================================================================
% C
% C     RAILP
% C     =====
% C
% C     PART OF THE TEST SET FOR INITIAL VALUE PROBLEMS, SEE
% C
% C        HTTP://WWW.CWI.NL/CWI/PROJECTS/IVPTESTSET.HTML
% C
% C     VERSION :        NOV. 1995
% C     AUTHORS :        SIMEON, FUEHRER, RENTROP
% C     PURPOSE :        EVALUATION OF PROFILE RAIL FOR
% C                      SIMULATION OF WHEELSET WITH DASSL
% C     DESCRIPTION:     SEE THE ABOVE MENTIONED TESTSET
% C     REFERENCE:       B. SIMEON, C. FUEHRER, P. RENTROP:
% C                      INTRODUCTION TO DIFFERENTIAL-ALGEBRAIC
% C                      EQUATIONS IN VEHICLE SYSTEM DYNAMICS,
% C                      SURV. MATH. IND. I: 1-37 (1991)
% C
% C     PARAMETER  (I=INPUT, O=OUTPUT)
% C
% C     X          DISPLACEMENT   (XI (LEFT OR RIGHT))              (I)
% C     SX         VALUE OF RAIL PROFILE FUNKTION   R(XI)                (O)
% C     DSX        ITS DERIVATIVE                                   (O)
% C     D2SX       ITS SECOND DERIVATIVE                            (O)
% C     D3SX       ITS THIRD DERIVATIVE (ONLY USED IN THE INDEX-1 CASE)
% C     IERR       ERROR CODE                                       (O)
% C                IERR =  0 : NO ERROR
% C                     = -2 : OUT OF RANGE
% C
% C
% C     CONSTANTS:
% C
% C     DELTA0     CONE ANGLE
% C     RN0        NOMINAL RADIUS
% C     AR         1/2 WHEEL DISTANCE
% C     RS         RAIL RADIUS
% C     EPS        TOLERANCE
% C
% C====================================================================
% 
%       DOUBLE PRECISION
%      *           X, SX, DSX, D2SX, D3SX
%       INTEGER    IERR
%       DOUBLE PRECISION
%      *           DELTA0, RN0, AR, RS, EPS, T, XABS, SIR
global LOUT

DELTA0 = 0.0262d0;
RN0    = 0.1000d0;
AR     = 0.1506d0;
RS     = 0.06d0;
EPS    = 0.00001d0;


      SIR  = sin(DELTA0)*RS;
      XABS = abs(X);
      if ( ( (AR + SIR + RS - EPS)<=XABS)||...
              ( (AR + SIR - RS + EPS)>=XABS )    )
          
           IERR = -2;
           %WRITE(6,*) 'OUT OF RAIL PROFILE (DERAILMENT)'
           fprintf(LOUT,'OUT OF RAIL PROFILE (DERAILMENT)');
      else
           T    = sqrt( RS*RS - (XABS-AR-SIR)^2 );
           SX   = T - RN0 - cos(DELTA0)*RS;
           
           if (X>=0)    
               S_I_G_N = 1;
           end
           
           if (X<0)
               S_I_G_N = -1;
           end
           
           %DSX  = SIGN(1.D0,X)*(-XABS+AR+SIR)/T;
           DSX  = S_I_G_N*(-XABS+AR+SIR)/T;
           D2SX = -RS*RS/(T*T*T);
           D3SX = 3*RS*RS/(T*T*T*T)*DSX;
           IERR = 0;
      end
      
%==========================================================================
%==========================================================================
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%                               END OF RAILP                             %
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%==========================================================================
%==========================================================================
end     


function [X,RX,DRX,D2RX,D3RX,IERR] = WHEELP(X,RX,DRX,D2RX,D3RX,IERR)

% C====================================================================
% C
% C     WHEELP
% C     ======
% C
% C     PART OF THE TEST SET FOR INITIAL VALUE PROBLEMS, SEE
% C
% C        HTTP://WWW.CWI.NL/CWI/PROJECTS/IVPTESTSET.HTML
% C
% C     VERSION :        NOV. 1995
% C     AUTHORS :        SIMEON, FUEHRER, RENTROP
% C     PURPOSE :        EVALUATION OF PROFILE WHEEL FOR
% C                      SIMULATION OF WHEELSET WITH DASSL
% C     DESCRIPTION:     SEE THE ABOVE MENTIONED TESTSET
% C     REFERENCE:       B. SIMEON, C. FUEHRER, P. RENTROP:
% C                      INTRODUCTION TO DIFFERENTIAL-ALGEBRAIC
% C                      EQUATIONS IN VEHICLE SYSTEM DYNAMICS,
% C                      SURV. MATH. IND. I: 1-37 (1991)
% C
% C ===============================================================
% C
% C     PARAMETER  (I=INPUT, O=OUTPUT)
% C
% C     X          DISPLACEMENT   (XI (LEFT OR RIGHT))              (I)
% C     RX         VALUE OF PROFILE FUNKTION   R(XI)                (O)
% C     DRX        ITS DERIVATIVE                                   (O)
% C     D2RX       ITS SECOND DERIVATIVE                            (O)
% C     D3RX       ITS THIRD DERIVATIVE (ONLY USED IN THE INDEX-1 CASE)
% C     IERR       ERROR CODE                                       (O)
% C                IERR =  0 : NO ERROR
% C                     = -1 : OUT OF RANGE
% C
% C
% C     CONSTANTS:
% C
% C     DELTA0     CONE ANGLE
% C     RN0        NOMINAL RADIUS
% C     AR         1/2 WHEEL DISTANCE
% C     B1         INNER WHEEL LIMIT
% C     B2         OUTER WHEEL LIMIT
% C
% C====================================================================
% C
%       DOUBLE PRECISION
%      *           X, RX, DRX, D2RX, D3RX
%       INTEGER    IERR
%       DOUBLE PRECISION
%      *           DELTA0, RN0, AR, B1, B2, TA, XABS
global LOUT

DELTA0 = 0.0262D0;
RN0    = 0.1000D0;
AR     = 0.1506D0;
B1     = 0.0000D0;
B2     = 4.000D0;

      TA   = tan(DELTA0);
      XABS = abs(X);
      if ( (B1>=XABS)||(XABS>=B2) )
          IERR = -1;
          %WRITE(6,*) 'OUT OF WHEEL PROFILE (DERAILMENT)'
          fprintf(LOUT,'OUT OF WHEEL PROFILE (DERAILMENT)');
      else
          RX   = RN0 + TA*(AR-XABS);
           if (X>=0)    
               S_I_G_N = 1;
           end
           
           if (X<0)
               S_I_G_N = -1;
           end
          
          %DRX  = SIGN(1.D0,X)*(-TA);
          DRX  = S_I_G_N*(-TA);
          D2RX = 0.0d0;
          D3RX = 0.0d0;
          IERR = 0;
      end
%==========================================================================
%==========================================================================
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%                               END OF WHEELP                             %
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%==========================================================================
%==========================================================================
end     




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEQN,T,Y,YPRIME] = INIT(NEQN,T,Y,YPRIME)
            
   
% C
% C     INITIAL CONDITIONS:
% C     POSITION VARIABLES  X, Y, Z, THETA, PHI   =>  P(0)


      Y( 1) =   0.14941D-02;
      Y( 2) =   0.40089D-06;
      Y( 3) =   0.11241D-05;
      Y( 4) =   -.28573D-03;
      Y( 5) =   0.26459D-03;
      
% C     VELOCITIES     =>   v(0) 

%for I=6:10
%        Y(I)=0.0D0;
%end
Y(6:10)=0;

% C     DYNAMIC VARIABLE BETA =>  B(0)
      Y(11) =   0.0d0;
% C     CONTACT COORDINATES PSIL, XIL, PSIR, XIR  => ( q )
% c     (determined numerically by Jacques de Swart, 23 October 1997)
      Y(12)= -7.4122380357667139d-06;
      Y(13)= -0.1521364296121248d0;
      Y(14)=  7.5634406395172940d-06;
      Y(15)=  0.1490635714733819d0;
%c.......... w(0) = v'(0)
      Y(16)= -1.975258894011285d0;
      Y(17)= -1.0898297102811276d-03;
      Y(18)=  7.8855083626142589d-02;
      Y(19)= -5.533362821731549d0;
      Y(20)= -0.3487021489546511d0;
%c..........  m(0) = B'(0)
      Y(21)= -2.132968724380927d0;
%C     LAGRANGE MULTIPLIERS (SCALED BY*1.E-4)  => lambda(0)
      Y(22) =   -.83593d-02;
      Y(23) =   -.74144d-02 ;     

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [neqn,t,y,df,ipar,rpar,ierr] = F(neqn,t,y,df,ipar,rpar,ierr)
%disp('now in F')
   
%  
% c     uses the routines
% c
% c          reswhs
% c          wheelp
% c          railp
% c          creep
% c          constm
% c
      
% c
% c we interchange y(12)-y(22) and y(13)-y(23)
% c
global YPRIME
      y12=y(12);
      y(12)=y(22);
      y(22)=y12;
      y13=y(13);
      y(13)=y(23);
      y(23)=y13;
% c     index-2 formulation: ipar(1)=0
% c
      ipar(1)=0;
      [t,y,YPRIME,df,ires,rpar,ipar] = ...
          RESWHS(t,y,YPRIME,df,rpar,ipar);
% c
% c we interchange y(12)-y(22) and y(13)-y(23)
% c
      y12=y(12);
      y(12)=y(22);
      y(22)=y12;
      y13=y(13);
      y(13)=y(23);
      y(23)=y13;

%disp('letf F')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEQN,AM,LDIM,IPAR,RPAR,IERR] = MAS(NEQN,AM,LDIM,IPAR,RPAR,IERR)
%disp('IN MAS')
  
 

      %DO J=1,11
      %   AM(1,J) = 1.0D0                  
      %ENDDO
AM(1,1:11) = 1;

     
%disp('LEFT MAS')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,y,dfy,NEQN,YPRIME,IPAR,RPAR,IERR]...
    = PDERV(t,y,dfy,NEQN,YPRIME,IPAR,RPAR,IERR)
%disp('now in pderv')

%DUMMY ROUTINE
     
%disp('left pderv')
end
