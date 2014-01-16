%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE CHEMICAL AKZO NOBEL PROBLEM: INDEX = 1.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


function CHEMAKZO
global MYtagSTIFF
MYtagSTIFF = 1;
ND=6;
LWORK=(38+3*ND)*ND+3;
LIWORK=ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);
YPRIME = zeros(1,ND);

MBND = zeros(1,4);
TRUE  = zeros(1,ND);
ERROR = zeros(1,ND);
MASBND = zeros(1,4);
MYVAR = 0;
RPAR = 0;
IPAR = 0;
IERR = 0;

 
      NEQN=ND;
      N = ND;
        MASBND(1)=1;
        MASBND(2)=0;
        MASBND(3)=0;
        MASBND(4)=1;
        MBND(1) = N;
        MBND(2) = N;
        MBND(3) = N;
        MBND(4) = N;
% C     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
% C     WORK() HOUSES THE FOLLOWING ARRAYS
% C
% C     Y(N,12)   , YHOLD(N,12) , YNHOLD(N,2) , YMAX(N)
% C     ERRORS(N) , SAVE1(N) , SAVE2(N) , SCALE(N) , ARH(N) , PW(MBND(4)*N)
% C     PWCOPY(MBND(4)*N) ,AM(MASBND(4)*N)
% C     IF THE SPARSE OPTION IS NOT BEING USED THEN MBND(4)=N=MASBND(4).
% C     <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
IWORK_15 = zeros(1,N);
WORK_1   = 0; %WORK(1)
WORK_2   = 0; %WORK(2)
WORK_3   = zeros(N,12); %WORK(3).....ETC
WORK_I1  = zeros(N,12);
WORK_I2  = zeros(N,2);
WORK_I3  = zeros(1,N);
WORK_I4  = zeros(1,N);
WORK_I5  = zeros(1,N);
WORK_I6  = zeros(1,N);
WORK_I7  = zeros(1,N);
WORK_I8  = zeros(1,N);
WORK_I9  = zeros(1,MBND(4)*N);
WORK_I10 = zeros(1,MBND(4)*N);
WORK_I11 = zeros(1,MASBND(4)*N);




      
LOUT = fopen('chemakzo.txt','w');
fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE CHEMAKZO PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r Number of Equations : %g \r',ND);

%... ENDPOINT OF INTEGRATION
    
         XEND=180;
         XOUT=XEND;
         X=0;
%...  REQUIRED TOLERANCE

        RTOL=10.0E-8;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES

        [NEQN,X,Y,YPRIME] = INIT(NEQN,X,Y,YPRIME,IPAR,RPAR,IERR);         
%...  SET DEFAULT VALUES
        MF=21;
        INDEX=1;
        MAXDER=7;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=100000;
        H=1E-8;
        fprintf(LOUT,'\r Hstart = %E',H);
        IWORK(1)=N;
        IWORK(2)=0;
        IWORK(3)=0;         
%...
%... CALL OF THE SUBROUTINE
%...

MYtag = 0;

while (MYtag == 0)
             
          [N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,WORK_3,...
              WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,WORK_I7,...
              WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,IWORK_15,...
              MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR] =MEBDF...
              (N,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,...
              WORK_3,WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,...
              WORK_I7,WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,...
              IWORK_15,MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR);  

           if (INDEX==1)
              INDEX = 0;
              
              continue
           else
             
              MYtag =34;
              
           end
           
end


      TRUE(1)= 0.1150794920661702;
      TRUE(2)= 0.1203831471567715E-2;
      TRUE(3)= 0.1611562887407974;
      TRUE(4)= 0.3656156421249283E-3;
      TRUE(5)= 0.1708010885264404E-1;
      TRUE(6)= 0.4873531310307455E-2;



      
for t = 1:N
         fprintf(LOUT,'\r\r APROX Y%g = %E',t,Y(t));
         fprintf(LOUT,'\rTRUE Y%g = %E\r\r',t,TRUE(t));
end      

      fclose(LOUT)
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEQN,T,Y,YPRIME] = INIT(NEQN,T,Y,YPRIME,IPAR,RPAR,IERR)
      
      ks   =115.83;

      Y(1) = 0.444;
      Y(2) = 0.00123;
      Y(3) = 0;
      Y(4) = 0.007;
      Y(5) = 0;
      Y(6) = ks*Y(1)*Y(4);

      [NEQN,T,Y,YPRIME,IPAR,RPAR,IERR] = F(NEQN,T,Y,YPRIME,IPAR,RPAR,IERR);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEQN,X,Y,DF,ipar,rpar,ierr] = F(NEQN,X,Y,DF,ipar,rpar,ierr)

%disp('now in F')

        k1   = 18.7;
        k2   = 0.58;
        k3   = 0.09;
        k4   = 0.42;
        kbig = 34.4;
        kla  = 3.3;
        ks   = 115.83;
        po2  = 0.9;
        hen  = 737;

      if (Y(2)< 0)
         ierr = -1
         return
      end
      
      r1  = k1*(Y(1)^4)*sqrt(Y(2));
      r2  = k2*Y(3)*Y(4);
      r3  = k2/kbig*Y(1)*Y(5);
      r4  = k3*Y(1)*(Y(4)^2);
      r5  = k4*(Y(6)^2)*sqrt(Y(2));
      fin = kla*(po2/hen-Y(2));
      
      DF(1) =   -2*r1 +r2 -r3     -r4;
      DF(2) = -0.5*r1             -r4     -0.5*r5 + fin;
      DF(3) =        r1 -r2 +r3;
      DF(4) =           -r2 +r3 -2*r4;
      DF(5) =            r2 -r3         +r5;
      DF(6) = ks*Y(1)*Y(4)-Y(6);

%disp('letf F')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NEQN,AM,LDIM,IPAR,RPAR,IERR] = MAS(NEQN,AM,LDIM,IPAR,RPAR,IERR)
%disp('IN MAS')

      for I = 1:(NEQN-1)
          AM(1,I)=1;
      end

      AM(1,NEQN)=0;

%disp('LEFT MAS')

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [T,Y,DFDY,NEQN,MEBAND,IPAR,RPAR,IERR] =  PDERV(T,Y,DFDY,NEQN,MEBAND,IPAR,RPAR,IERR)
%disp('now in pderv')

    
        k1   =18.7;
        k2   =0.58;
        k3   =0.09;
        k4   =0.42;
        kbig =34.4;
        kla  =3.3;
        ks   =115.83;

   
      if (Y(2)<0)
         IERR = -1;
         return
      end
      
      r11  = 4*k1*(Y(1)^3)*sqrt(Y(2));
      r12  = 0.5*k1*(Y(1)^4)/sqrt(Y(2));
      r23  = k2*Y(4);
      r24  = k2*Y(3);
      r31  = (k2/kbig)*Y(5);
      r35  = (k2/kbig)*Y(1);
      r41  = k3*Y(4)^2;
      r44  = 2*k3*Y(1)*Y(4);
      r52  = 0.5*k4*(Y(6)^2)/sqrt(Y(2));
      r56  = 2*k4*Y(6)*sqrt(Y(2));
      fin2 = -kla;
      DFDY(1,1) = -2*r11-r31-r41;
      DFDY(1,2) = -2*r12;
      DFDY(1,3) = r23;
      DFDY(1,4) = r24-r44;
      DFDY(1,5) = -r35;
      DFDY(2,1) = -0.5*r11-r41;
      DFDY(2,2) = -0.5*r12-0.5*r52+fin2;
      DFDY(2,4) = -r44;
      DFDY(2,6) = -0.5*r56;
      DFDY(3,1) = r11+r31;
      DFDY(3,2) = r12;
      DFDY(3,3) = -r23;
      DFDY(3,4) = -r24;
      DFDY(3,5) = r35;
      DFDY(4,1) = r31-2*r41;
      DFDY(4,3) = -r23;
      DFDY(4,4) = -r24-2*r44;
      DFDY(4,5) = r35;
      DFDY(5,1) = -r31;
      DFDY(5,2) = r52;
      DFDY(5,3) = r23;
      DFDY(5,4) = r24;
      DFDY(5,5) = -r35;
      DFDY(5,6) = r56;
      DFDY(6,1) = ks*Y(4);
      DFDY(6,4) = ks*Y(1);
      DFDY(6,6) = -1;

%disp('left pderv')
end

