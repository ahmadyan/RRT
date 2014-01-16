%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE SLIDER CRANK PROBLEM: INDEX = 2.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function CRANK
global MYtagSTIFF FIRST

%FOR USE IN RESMBS%
FIRST = 1;

MYtagSTIFF = 1;

ND=24;
NEQN = ND;
N = ND;

LWORK=(38+3*ND)*ND+3;
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
      
LOUT = fopen('crank.txt','w');
fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE CRANK PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r Number of Equations : %g \r',ND);



%... ENDPOINT OF INTEGRATION
         X=0;
         XEND=0.1;
         XOUT=XEND;
         
%...  REQUIRED TOLERANCE

        %RTOL=.316227766016837972E-06;
        RTOL=1E-07;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
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
        H=1.0E-2*RTOL;
         IWORK(1)=14;
         IWORK(2)=10;
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

function  YEND = SOLN(dummy)         
      YEND(1) = 1.500000000000000E+01;
      YEND(2) = -3.311734987910881E-01;
      YEND(3) = 1.697373326718410E-01;
      YEND(4) = 1.893192460247178E-04;
      YEND(5) = 2.375751865617931E-05;
      YEND(6) = -5.323907988763734E-06;
      YEND(7) = -8.363283141616840E-06;
      YEND(8) = 1.500000000000000E+02;
      YEND(9) = 6.025346682645789E+01;
      YEND(10)= -8.753116989887888E+00;
      YEND(11)= -3.005536801092212E-02;
      YEND(12)= -5.500488291932075E-03;
      YEND(13)= 4.978243404809343E-04;
      YEND(14)= 1.104933470696396E-03;
      YEND(15)= 0.E0;
      YEND(16)= 6.488722210234531E+03;
      YEND(17)= 2.167924253080623E+03;
      YEND(18)= 3.391435115267547E+01;
      YEND(19)= 1.699107480197843E-01;
      YEND(20)= -1.415799354959001E+00;
      YEND(21)= 9.903251655235532E-01;
      YEND(22)= -6.232893262533717E+01;
      YEND(23)= -1.637910131687472E+02;
      YEND(24)= 2.529853213732781E+01;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [ITYP, IEQUA, ICALL,T,X,XD,DELTA,IRES,RPAR,IPAR] = RESMBS(ITYP, IEQUA, ICALL,T,X,XD,DELTA,IRES,RPAR,IPAR)
  
% C
% C     Slider crank - flexible model with sliding block
% C     ------------------------------------------------
% C  ** written by Bernd Simeon, TH Darmstadt, 13/06/95 **
% C  ** extended version with special beam model,
% C  ** for IVPTestset                         11/28/97 **
% C
% C     The flexible slider crank mechanism is described in
% C
% C     Simeon, B.: Modelling a Flexible Slider Crank Mechanism
% C     by a Mixed System of DAEs and PDEs.
% C     Math. Modelling of Systems 2, 1-18, 1995
% C
% C     This version contains all coupling terms (also for 2D)
% C     for discretizations of the connecting rod.
% C     Particular grid used: 2 modal functions lateral;
% C                           2 nodes (quadratic ansatz) longitudinal.
% C
% C     phi1(t) = omega*t prescribed by constraint.
% C
% C     PARAMETERS ON ENTRY:
% C
% C     ITYP     This integer flag determines in which formulation the equations
% C              of motion have to be evaluated:
% C              ITYP = 0: index 3 system;     ITYP = 3: index 1 system;
% C                   = 1: index 2 system;
% C
% C     IEQUA    This integer flag determines whether the complete residual
% C              or only parts of it have to be evaluated.
% C              IEQUA = 0: evaluate complete residual;
% C                    = 2: evaluate only position+velocity constraints
% C                         in DELTA(1:6)
% C
% C     ICALL    This integer flag indicates whether RESMBS has already
% C              been called with the actual parameter set T, X, XD.
% C              ICALL = 0: new values T, X, XD;
% C                    = 1: T, X, XD have not changed since the last call.
% C              (unused)
% C
% C     T        This real variable contains the current value of the
% C              independent variable (time).
% C
% C     X(*)     This array contains the current values of the dependent
% C              variables: The multibody system variables are arranged as
% C              X = [ p                   rigid motion (3 coordinates)
% C                    q                   elastic motion (4 nodes)
% C                    pd                  velocity variables p
% C                    qd                  velocity variables q
% C                    w                   acceleration variables p, q
% C                    lambda              3 Lagrange multipliers
% C
% C     XD(*)    This array contains the derivatives of the solution
% C              components at T.
% C
% C     RPAR,IPAR These are real and integer parameter arrays which
% C              are used for communication between the calling program
% C              and this subroutine.
% C       IPAR(1)  0: only linear stiffness term K*q
% C                1: nonlinear stiffness term included
% C
% C       IPAR(2)  0: no physical damping - purely imaginary EV
% C                1: damping matrix (1 %) included
% C
% C     ON RETURN:
% C
% C     DELTA(NEQ) This array contains the residual of the equations of motion.
% C
% C     IRES     This integer flag indicates a stop condition (unused).
% C
% C     SYSTEM PARAMETERS:
% C
% C       IPAR(1)  0: only linear stiffness term K*q
% C                1: nonlinear stiffness term included
% C
% C       IPAR(2)  0: no physical damping - purely imaginary EV
% C                1: damping matrix (0.5 %) included
% C
      
% C
% C     Data set
% C
    M1 = 0.36E0;
    M2 = 0.151104E0;
    M3 = 0.075552;
    L1 = 0.15E0;     
    L2 = 0.30E0;
    J1 = 0.002727E0;
    J2 = 0.0045339259E0;
    PI = 3.1415927E0;
    EE = .20E12;
    NUE= 0.30E0;
    BB = 0.0080E0;
    HH = 0.0080E0;
    RHO= 7870.0E0;
    GRAV= 0.0E0;
    OMEGA = 150.E0;  
    
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
%                               FROM ORIGINAL
% C
% C     LOCAL VARIABLES:
% C
% C       Q, QD for FE coefficients and time derivatives,
% C       MQ, KQ, DQ, BQ, c1, c2, c12, c21 for FE matrices and vectors,
% C       up to NQMAX = 20 variables.
% C
%         INTEGER  NQMAX, I, J, JJ, NQ, NP, NL, NX
%         PARAMETER( NQMAX = 20 )
%         DOUBLE PRECISION
%      *           Q(NQMAX), QD(NQMAX), MQ(NQMAX,NQMAX),
%      *           KQ(NQMAX,NQMAX), BQ(NQMAX,NQMAX), DQ(NQMAX,NQMAX),
%      *           c1(NQMAX), c2(NQMAX), c12(NQMAX), c21(NQMAX),
%      *           MQQ(NQMAX), KQQ(NQMAX), DQQD(NQMAX),
%      *           QTBQ(NQMAX), BQQD(NQMAX),
%      *           c1TQ, c1TQD, c2TQ, c2TQD, c12TQ, c12TQD,
%      *           QDTBQQD, QTMQQ, QDTMQQ, DDOT, V(2),
%      *           ALC(3), PLC(3), VLC(3),
%      *           AM(NQMAX+3,NQMAX+3), GP(3,NQMAX+3), F(NQMAX+3),
%      *           COSP1, COSP2, SINP1, SINP2, COSP12, SINP12,
%      *           QKU, QKV, QDKU, QDKV, FACM, FACK, FACB
%         SAVE     MQ, KQ, DQ, BQ, c1, c2, c12, c21
% C
% C       FIRST for first call - evaluation of FE matrices.
% C
%         LOGICAL  FIRST
%         DATA     FIRST / .TRUE. /
% C
% C_______________End of declaration part RESMBS____________________________
% C
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^



% C
% C       FIRST for first call - evaluation of FE matrices.
% C
global FIRST
persistent MQ KQ DQ BQ c1 c2 c12 c21

        NQMAX = 20;
        %WILL BE DECLARED OUTSIDE IN DRIVER ROUTINE.
        %FIRST = 1;
        IRES  = 0;
        NQ    = 4;
        NP    = 7;
        NL    = 3;
        NX    = 3*NP + NL;
        KU    = 4;
        KV    = 0;
        AM  = zeros(NQMAX+3);
        
      if (FIRST) 
        
% C
% C       Initialize grid data.
% C


%DECLERATIONS ON FIRST CALL ONLY (SO NOT TO RESET).
%*************************************************%
MQ  = zeros(NQMAX);
KQ  = zeros(NQMAX);
DQ  = zeros(NQMAX);
BQ  = zeros(NQMAX);
c1  = zeros(1,NQMAX);
c2  = zeros(1,NQMAX);
c12 = zeros(1,NQMAX);
c21 = zeros(1,NQMAX);

%*************************************************%

        FACM = RHO*BB*HH*L2;
        FACK = EE*BB*HH/L2;
        FACB = BB*HH*L2;

        for I=1:NQ
           for J=1:NQ
              MQ(J,I) = 0;
              KQ(J,I) = 0;
              BQ(J,I) = 0;
              DQ(J,I) = 0;
           end
           c1(I) = 0;
           c2(I) = 0;
           c12(I)= 0;
           c21(I)= 0;
        end

        MQ(1,1) = FACM*.5;
        MQ(2,2) = FACM*.5;
        MQ(3,3) = FACM*8;
        MQ(3,4) = FACM*1;
        MQ(4,3) = FACM*1;
        MQ(4,4) = FACM*2;

        KQ(1,1) = FACK*PI^4/24*(HH/L2)^2;
        KQ(2,2) = FACK*PI^4*2/3*(HH/L2)^2;
        KQ(3,3) = FACK*16/3;
        KQ(3,4) = -FACK*8/3;
        KQ(4,3) = -FACK*8/3;
        KQ(4,4) = FACK*7/3;

        BQ(1,3) = -FACB*16/PI^3;
        BQ(1,4) = FACB*(8/PI^3-1/PI);
        BQ(2,4) = FACB*0.5/PI;
        BQ(3,1) = FACB*16/PI^3;
        BQ(4,1) = -FACB*(8/PI^3-1/PI);
        BQ(4,2) = -FACB*0.5/PI;

        c1(3)  = FACB*2/3;
        c1(4)  = FACB*1/6;
        c2(1)  = FACB*2/PI;
        c12(3) = L2*FACB*1/3;
        c12(4) = L2*FACB*1/6;
        c21(1) = L2*FACB*1/PI;
        c21(2) = -L2*FACB*0.5/PI;

        if (IPAR(2)==1)
% C
% C       0.5 per cent damping
% C
           DQ(1,1) = 5;
           DQ(2,2) = 25;
           DQ(3,3) = 0.5*2.308375455264791E+02;
           DQ(3,4) = -0.5*2.62688487992052E+02;
           DQ(4,3) = -0.5*2.626884879920526E+02;
           DQ(4,4) = 0.5*4.217421837156818E+02;
        end
        
        FIRST = 0;
      end

      COSP1  = cos(X(1));
      COSP2  = cos(X(2));
      SINP1  = sin(X(1));
      SINP2  = sin(X(2));
      COSP12 = cos(X(1)-X(2));
      SINP12 = sin(X(1)-X(2));
      V(1)   = X(NP+1);
      V(2)   = X(NP+2);

      for I=1:NQ
         Q(I)  = X(3+I);
         QD(I) = X(NP+3+I);
      end
      
% C
% C     Evaluate scalar products and quadratic forms.
% C

%       c1TQ  = DDOT(NQ,c1,1,Q,1)
%       c1TQD = DDOT(NQ,c1,1,QD,1)
%       c2TQ  = DDOT(NQ,c2,1,Q,1)
%       c2TQD = DDOT(NQ,c2,1,QD,1)
%       c12TQ = DDOT(NQ,c12,1,Q,1)
%       c12TQD= DDOT(NQ,c12,1,QD,1)
%       DO 10 I=1,NQ
%          MQQ(I) = DDOT(NQ,MQ(1,I),1,Q,1)
%          KQQ(I) = DDOT(NQ,KQ(1,I),1,Q,1)
%          DQQD(I)= DDOT(NQ,DQ(1,I),1,QD,1)
%          QTBQ(I)= DDOT(NQ,Q,1,BQ(1,I),1)
%          BQQD(I)= DDOT(NQ,BQ(I,1),NQMAX,QD,1)
%   10  CONTINUE
%       QTMQQ   = DDOT(NQ,Q,1,MQQ,1)
%       QDTMQQ  = DDOT(NQ,QD,1,MQQ,1)
%       QDTBQQD = DDOT(NQ,QD,1,BQQD,1)


%THESE ARE NEEDED TO IMPLEMENT THE MATLAB DOT%
        %*************************%
        %Q
        Q_temp(1:NQ)    = Q(1:NQ);
        %c1
        c1_temp(1:NQ)   = c1(1:NQ);
        %c2
        c2_temp(1:NQ)   = c2(1:NQ);
        %c3
        c12_temp(1:NQ)  = c12(1:NQ);
        %QD
        QD_temp(1:NQ)   = QD(1:NQ);
        %**************************%
        
        
        c1TQ  = dot(c1_temp,Q_temp);
        c1TQD = dot(c1_temp,QD_temp);
        c2TQ  = dot(c2_temp,Q_temp);
        c2TQD = dot(c2_temp,QD_temp);
        c12TQ = dot(c12_temp,Q_temp);
        c12TQD= dot(c12_temp,QD_temp);
       
        
        
      for I=1:NQ
          
                %THESE ARE NEEDED TO IMPLEMENT THE MATLAB DOT%
          %*******************************************************%
          %MQ.  
          MQ_temp(1:NQ) = MQ((I-1)*NQMAX+1:(I-1)*NQMAX+NQ);
          %KQ. 
          KQ_temp(1:NQ) = KQ((I-1)*NQMAX+1:(I-1)*NQMAX+NQ);
          %DQ.
          DQ_temp(1:NQ) = DQ((I-1)*NQMAX+1:(I-1)*NQMAX+NQ);
          %BQ. 
          BQ_temp1(1:NQ) = BQ((I-1)*NQMAX+1:(I-1)*NQMAX+NQ);
          %BQ. (HERE WE WANT ROWS see below).
          BQ_temp2A = BQ(I,:);
          BQ_temp2B(1:NQ) = BQ_temp2A(1:NQ);
          %********************************************************%
          
          
          MQQ(I) = dot(MQ_temp,Q_temp);
          
          KQQ(I) = dot(KQ_temp,Q_temp);
          
          DQQD(I)= dot(DQ_temp,QD_temp);
          
          QTBQ(I)= dot(Q_temp,BQ_temp1);
          
          %THIS IS DOTing THE ROWS NOT COLUMNS%
          BQQD(I)= dot(BQ_temp2B,QD_temp);
          %BQQD(I)= DDOT(NQ,BQ(I,1),NQMAX,QD,1)



         
      end
      
%       MQQ
%       KQQ
%       DQQD
%       QTBQ
%       vpa(BQQD)
%       pause
      
      QTMQQ   = dot(Q_temp,MQQ);
      QDTMQQ  = dot(QD_temp,MQQ);
      QDTBQQD = dot(QD_temp,BQQD);

% C
% C     Kinematic and dynamic equations.
% C

      for I=1:NP
         DELTA(I)    = XD(I)    - X(NP+I);
         DELTA(NP+I) = XD(NP+I) - X(2*NP+I);
      end
      
% C
% C     Compute mass matrix.
% C

          AM(1,1) = J1 + M2*L1*L1;
          AM(1,2) = .5*L1*L2*M2*COSP12;
          AM(2,2) = J2;
          AM(1,3) = 0;
          AM(2,3) = 0;
          AM(3,1) = 0;
          AM(3,2) = 0;
          AM(3,3) = M3;
          AM(1,2) = AM(1,2) + RHO*L1*(SINP12*c2TQ+COSP12*c1TQ);
          AM(2,2) = AM(2,2) + QTMQQ + 2*RHO*c12TQ;
          
          
          for I=1:NQ
             AM(1,3+I) = RHO*L1*(-SINP12*c1(I) + COSP12*c2(I));
             AM(2,3+I) = RHO*c21(I) + RHO*QTBQ(I);
             AM(3,3+I) = 0;
          end
          
          for I=1:NQ
             for J=1:I
                AM(3+J,3+I) = MQ(J,I);
             end
          end
          
          for I=1:NP
             for J=(I+1):NP
                AM(J,I) = AM(I,J);
             end
          end
          
% C
% C     Compute constraint matrix.
% C

          if (KU==0) 
             QKU = 0;
          else
             QKU = Q(KU);
          end
          
          if (KV==0)
             QKV = 0;
          else
             QKV = Q(KV);
          end
          
          GP(1,1) = L1*COSP1;
          GP(1,2) = L2*COSP2 + QKU*COSP2 - QKV*SINP2;
          GP(1,3) = 0;
          GP(2,1) = L1*SINP1;
          GP(2,2) = L2*SINP2 + QKU*SINP2 + QKV*COSP2;
          GP(2,3) = 1;
          GP(3,1) = 1;
          GP(3,2) = 0;
          GP(3,3) = 0;
          
          for I=1:NQ
             GP(1,3+I) = 0;
             GP(2,3+I) = 0;
             GP(3,3+I) = 0;
          end
          
          if (KU~=0)
             GP(1,3+KU) = SINP2;
             GP(2,3+KU) = -COSP2;
          end
          
          if (KV~=0)
             GP(1,3+KV) = COSP2;
             GP(2,3+KV) = SINP2;
          end
          
% C
% C     Forces - rigid motion entries.
% C

          F(1) = -.5*L1*GRAV*(M1+2*M2)*COSP1-.5*L1*L2*M2*V(2)*V(2)*SINP12;
          F(2) = -.5*L2*GRAV*M2*COSP2+.5*L1*L2*M2*V(1)*V(1)*SINP12;
          F(3) = 0;
% C
% C     Superposition of flexible motion (term f^e).
% C

          F(1) = F(1)+ RHO*L1*V(2)*V(2)*(-SINP12*c1TQ+COSP12*c2TQ)-...
              2*RHO*L1*V(2)*(COSP12*c1TQD+SINP12*c2TQD);
          F(2) = F(2)+ RHO*L1*V(1)*V(1)*(SINP12*c1TQ-COSP12*c2TQ)-...
              2*RHO*V(2)*c12TQD - 2*V(2)*QDTMQQ-...
              RHO*QDTBQQD - RHO*GRAV*(COSP2*c1TQ-SINP2*c2TQ);
% C
% C     Coriolis and gravity terms flexible motion (Gamma).
% C

          for I=1:NQ
             F(3+I) = V(2)*V(2)*MQQ(I)+ RHO*(V(2)*V(2)*c12(I)+...
                 L1*V(1)*V(1)*(COSP12*c1(I)+SINP12*c2(I))+...
                 2*V(2)*BQQD(I) )- RHO*GRAV*(SINP2*c1(I)+COSP2*c2(I));
          end
          
% C
% C         Stiffness + damping terms - K q - D q'.
% C

          for I=1:NQ
             F(3+I) = F(3+I) - KQQ(I) - DQQD(I);
          end
          
          if (IPAR(1)==1)
% C
% C            Nonlinear stiffness term
% C

             FACK = 0.5*EE*BB*HH/L2^2*PI^2;
             FACB = 80.D0/(PI^2*9);
             F(4) = F(4) -FACK*(Q(1)*Q(4)-FACB*Q(2)*(-4*Q(3)+2*Q(4)));
             F(5) = F(5) -FACK*(4*Q(2)*Q(4)-FACB*Q(1)*(-4*Q(3)+2*Q(4)));
             F(6) = F(6) -FACK*4*FACB*Q(1)*Q(2);
             F(7) = F(7) -FACK*(0.5*Q(1)^2+2*Q(2)^2-2*FACB*Q(1)*Q(2));
         
          end
          
% C
% C     Dynamics part II ( M*w - f + G(T)*lambda ).
% C

 
         for I=1:NP
 
             AM_temp(1:NP) = AM((I-1)*(NQMAX+3)+1:(I-1)*(NQMAX+3)+NP);
            
             X_temp(1:NP) = X((2*NP+1):(2*NP+1)+NP-1);
             
         DELTA(2*NP+I) = dot(AM_temp,X_temp)- F(I) + GP(1,I)*X(NX-2)+GP(2,I)*X(NX-1)+GP(3,I)*X(NX);
         end
 

% C
% C     Acceleration level constraints.
% C


      if (KU==0)
          QDKU = 0;
      else
          QDKU = QD(KU);
      end
      
      if (KV==0)
          QDKV = 0;
      else
          QDKV = QD(KV);
      end
      
      ALC(1) = -L1*SINP1*V(1)*V(1) - (L2+QKU)*SINP2*V(2)*V(2)...
          +2*V(2)*(COSP2*QDKU-SINP2*QDKV) - COSP2*V(2)*V(2)*QKV;
      ALC(2) =  L1*COSP1*V(1)*V(1) + (L2+QKU)*COSP2*V(2)*V(2)...
          +2*V(2)*(SINP2*QDKU+COSP2*QDKV) - SINP2*V(2)*V(2)*QKV;
      ALC(3) = 0;
      
      for I=1:NP
         ALC(1) = ALC(1) + GP(1,I)*X(2*NP+I);
         ALC(2) = ALC(2) + GP(2,I)*X(2*NP+I);
         ALC(3) = ALC(3) + GP(3,I)*X(2*NP+I);
      end
      
% C
% C     Position level constraints.
% C


      PLC(1) = L1*SINP1 + L2*SINP2 + QKU*SINP2 + QKV*COSP2;
      PLC(2) = X(3) - L1*COSP1 - L2*COSP2-QKU*COSP2 + QKV*SINP2;
      PLC(3) = X(1) - OMEGA*T;
      
% C
% C     Velocity level constraints.
% C

      VLC(1) = 0;
      VLC(2) = 0;
      VLC(3) = -OMEGA;
      
      for I=1:NP
         VLC(1) = VLC(1) + GP(1,I)*X(NP+I);
         VLC(2) = VLC(2) + GP(2,I)*X(NP+I);
         VLC(3) = VLC(3) + GP(3,I)*X(NP+I);
      end
      

      if (IEQUA==2)
% C
% C         Evaluate only the constraints.
% C
          DELTA(1) = PLC(1);
          DELTA(2) = PLC(2);
          DELTA(3) = PLC(3);
          DELTA(4) = VLC(1);
          DELTA(5) = VLC(2);
          DELTA(6) = VLC(3);
      else
% C
% C         Select constraints defined by ITYP.
% C
          if (ITYP==0)
% C
% C             Index 3 system.
% C

              DELTA(NX-2) = PLC(1);
              DELTA(NX-1) = PLC(2);
              DELTA(NX)   = PLC(3);
              
          elseif (ITYP==1)
% C
% C             Index 2 system.
% C
              DELTA(NX-2) = VLC(1);
              DELTA(NX-1) = VLC(2);
              DELTA(NX)   = VLC(3);
              
          elseif (ITYP==3)
% C
% C             Index 1 system.
% C

              DELTA(NX-2) = ALC(1);
              DELTA(NX-1) = ALC(2);
              DELTA(NX)   = ALC(3);
              
          end
          
      end
      
% C
% C_______________End of subroutine RESMBS____________________________
% C
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEQN,T,y,YPRIME] = INIT(NEQN,T,y,YPRIME)
            
%  
% C
% C     Initial values: 'Close' to smooth motion,
% C     accelerations and multipliers consistent
% C     for linear stiffness term and no damping
% C     (ipar(1) = 0, ipar(2) = 0).
% C
% C     Position variables
% C     phi1, phi2, x3, q1, q2, q3, q4

      y(1) = 0;
      y(2) = 0;
      y(3) = .450016933;
      y(4) = 0;
      y(5) = 0;
      y(6) = .103339863E-04;
      y(7) = .169327969E-04;
%C     Initial values velocity variables
      y(8) =  .150000000E+03;
      y(9) = -.749957670E+02;
      y(10)= -.268938672E-05;
      y(11)=  .444896105E+00;
      y(12)=  .463434311E-02;
      y(13)= -.178591076E-05;
      y(14)= -.268938672E-05;
%C     Initial values acceleration variables
      y(15)= 0;
      y(16)= -1.344541576008661E-03;
      y(17)= -5.062194923138079E+03;
      y(18)= -6.833142732779555E-05;
      y(19)=  1.449382650173157E-08;
      y(20)= -4.268463211410861E+00;
      y(21)=  2.098334687947376E-01;
%C     Lagrange multipliers
      y(22)= -6.397251492537153E-08;
      y(23)=  3.824589508329281E+02;
      y(24)= -4.376060460948886E-09;
 
%      for i=1:24
%          YPRIME(i) = 0;
%      end
YPRIME(1:24) = 0;
     
    %  for i=1:14
    %      YPRIME(i) = y(7+i);
    % end
YPRIME(1:14) = y(8:21);
     
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEQN,T,Y,DY,IPAR,RPAR,IERR] = F(NEQN,T,Y,DY,IPAR,RPAR,IERR)
%disp('now in F')

        ub=6;
        uf=0.026;
        alpha=0.99;
        beta=1E-6;
        r0=1000;r1=9000;r2=9000;r3=9000;
        r4=9000;r5=9000;r6=9000;r7=9000;
        r8=9000;r9=9000;

      uet   = 0.1*sin(200*pi*T);
      fac1  = beta*(exp((Y(2)-Y(3))/uf)-1);
      fac2  = beta*(exp((Y(5)-Y(6))/uf)-1);
      
      DY(1) = (Y(1)-uet)/r0;
      DY(2) = Y(2)/r1+(Y(2)-ub)/r2+(1-alpha)*fac1;
      DY(3) = Y(3)/r3-fac1;
      DY(4) = (Y(4)-ub)/r4+alpha*fac1;
      DY(5) = Y(5)/r5+(Y(5)-ub)/r6+(1-alpha)*fac2;
      DY(6) = Y(6)/r7-fac2;
      DY(7) = (Y(7)-ub)/r8+alpha*fac2;
      DY(8) = Y(8)/r9;


%disp('letf F')
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEQN,AM,LDIM,IPAR,RPAR,IERR] = MAS(NEQN,AM,LDIM,IPAR,RPAR,IERR)
%disp('IN MAS')
   
     AM(1,1:14)  = 1;                  
 
     AM(1,15:24) = 0;
     
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
