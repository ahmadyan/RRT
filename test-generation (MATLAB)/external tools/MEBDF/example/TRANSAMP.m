%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE TRANSISTOR AMPLIFIER PROBLEM: INDEX = 1.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


function TRANSAMP
global MYtagSTIFF
MYtagSTIFF = 1;
ND=8;
LWORK=(41+2*ND)*ND+3;
LIWORK=ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);
YPRIME = zeros(1,ND);

MBND = zeros(1,4); 
TRUE  = zeros(1,ND);
ERROR = zeros(1,ND);
MASBND = zeros(1,4);
MASBND(1)=1;
        MASBND(2)=1;
        MASBND(3)=1;
        MASBND(4)=3;
        MBND(1) = 2;
        MBND(2) = 1;
        MBND(3) = 4;
        MBND(4) = 6;
MYVAR = 0;
RPAR = 0;
IPAR = 0;
IERR = 0;

%...  DIMENSION OF THE SYSTEM
      N = ND;
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
      
LOUT = fopen('transamp.txt','w');
fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE TRANSISTOR AMP PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r Number of Equations : %g \r',ND);

%... ENDPOINT OF INTEGRATION
    
         XEND=0.2;
         XOUT=XEND;
         X=0;
%...  REQUIRED TOLERANCE

        RTOL=1E-7;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES
        [N,X,Y] = INIT(N,X,Y);
%...  SET DEFAULT VALUES
        MF=23;
        INDEX=1;
        MAXDER=7;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=200000;
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


      TRUE(1) = -0.556214501E-02;
      TRUE(2) =  0.3006522473125E+01;
      TRUE(3) =  0.284995878984E+01;
      TRUE(4) =  0.29264225362E+01;
      TRUE(5) =  0.27046178656E+01;
      TRUE(6) =  0.2761837776452E+01;
      TRUE(7) =  0.4770927631E+01;
      TRUE(8) =  0.1236995867E+01;




      
for t = 1:N
         fprintf(LOUT,'\r\r APROX Y%g = %E',t,Y(t));
         fprintf(LOUT,'\rTRUE Y%g = %E\r\r',t,TRUE(t));
end      

      fclose(LOUT)
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,T,Y] = INIT(N,T,Y)

        ub=6;
        r1=9000;
        r2=9000;
        r3=9000;
        r5=9000;
        r6=9000;
        r7=9000;
        c2=2E-6;
        c4=4E-6;

      Y(1) = 0;
      Y(2) = ub/(r2/r1+1);
      Y(3) = Y(2);
      Y(4) = ub;
      Y(5) = ub/(r6/r5+1);
      Y(6) = Y(5);
      Y(7) = Y(4);
      Y(8) = 0;
      
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

function [neqn,am,ldim,ipar,rpar,ier] = MAS(neqn,am,ldim,ipar,rpar,ier)
%disp('IN MAS')
  
     c1=1E-6;
     c2=2E-6;
     c3=3E-6;
     c4=4E-6;
     c5=5E-6;

      for i=1:neqn
         am(1,i) = 0;
         am(3,i) = 0;
      end
      
      am(1,2) = c1;
      am(1,5) = c3;
      am(1,8) = c5;
      am(2,1) = -c1;
      am(2,2) = -c1;
      am(2,3) = -c2;
      am(2,4) = -c3;
      am(2,5) = -c3;
      am(2,6) = -c4;
      am(2,7) = -c5;
      am(2,8) = -c5;
      am(3,1) = c1;
      am(3,4) = c3;
      am(3,7) = c5;
      
  
%disp('LEFT MAS')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,y,jac,NEQN,ldim,ipar,rpar,ier] = pderv(t,y,jac,NEQN,ldim,ipar,rpar,ier)
%disp('now in pderv')
          
uf=0.026;alpha=0.99;beta=1E-6;
r0=1000;r1=9000;r2=9000;r3=9000;
r4=9000;r5=9000;r6=9000;r7=9000;
r8=9000;r9=9000;

      fac1p = beta*exp((y(2)-y(3))/uf)/uf;
      fac2p = beta*exp((y(5)-y(6))/uf)/uf;

      for i=1:8
         jac(1,i) = 0;
         jac(3,i) = 0;
         jac(4,i) = 0;
      end

      jac(1,3) = -(1-alpha)*fac1p;
      jac(1,6) = -(1-alpha)*fac2p;
      jac(2,1) = 1/r0;
      jac(2,2) = 1/r1+1/r2+(1-alpha)*fac1p;
      jac(2,3) = 1/r3+fac1p;
      jac(2,4) = 1/r4;
      jac(2,5) = 1/r5+1/r6+(1-alpha)*fac2p;
      jac(2,6) = 1/r7+fac2p;
      jac(2,7) = 1/r8;
      jac(2,8) = 1/r9;
      jac(3,2) = -fac1p;
      jac(3,3) = -alpha*fac1p;
      jac(3,5) = -fac2p;
      jac(3,6) = -alpha*fac2p;
      jac(4,2) = alpha*fac1p;
      jac(4,5) = alpha*fac2p;

%disp('left pderv')
end
