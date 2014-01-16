%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE VAN DER POL EQUATION. WITH F, MAS AND PDERV TO BE 
% SAVED AS SEPARATE M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
function VANDER

global MYtagSTIFF
MYtagSTIFF = 1;

ND=2;
LWORK =(33+3*ND)*ND+3;
LIWORK = ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);

MBND = zeros(1,4);
TRUE  = zeros(1,ND);
ERROR = zeros(1,ND);
MASBND = zeros(1,4);
MYY = zeros(1,20);
MYVAR = 0;
RPAR = 0;
IPAR = 0;
IERR = 0;

N=2;
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
WORK_I9  = zeros(1,N*N);
WORK_I10 = zeros(1,N*N);
WORK_I11 = zeros(1,MASBND(4)*N);
      
LOUT = fopen('VANDER.txt','w');

fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE VAN DER POL PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r Number of Equations : %g \r',ND);

%... ENDPOINT OF INTEGRATION
         X=0;
         XOUT=1;
         XEND=11;
%...  REQUIRED TOLERANCE
        RTOL=1E-7;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);       
%...  INITIAL VALUES
        Y(1)=2;
        Y(2)=0;
%...  SET DEFAULT VALUES
        MF=21;
        INDEX=1;
        MAXDER=7;
%...  MAXIMAL NUMBER OF STEPS
%       IWORK(1) = 2;
        IWORK(14)=200000;
        H=1E-8;
        MASBND(1)=0;
        fprintf(LOUT,'Hstart = %E',H);
%...
%... CALL OF THE SUBROUTINE
%...
 MYtag = 0;
 MYcase = 220;
 for I = 1:11
    
     MYtag = 0;
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
            
                XOUT=XOUT+1;
                MYY1(I) = Y(1);
                MYY2(I) = Y(2);
                XX(I) = X(1);
                MYtag = 34;
          continue
                
        end

     end

 end

for K = 1:11
      MYYTRUE = SOLN(K);
      
      fprintf(LOUT,'\r \r APPROX Y1 =  %E', MYY1(K));
      fprintf(LOUT,'\r TRUE Y1 =  %E', MYYTRUE(1));
      
      fprintf(LOUT,'\r \r APPROX Y2 =  %E', MYY2(K));
      fprintf(LOUT,'\r TRUE Y2 =  %E', MYYTRUE(2));
      fprintf(LOUT,'\r \r AT X =  %E', XX(K));
end
      fclose(LOUT)
     
end

function TRUE = SOLN(I)
      
switch I
    case 1
         TRUE(1)=-0.1863646254808130E+01;
         TRUE(2)= 0.7535430865435460E+00;
    case 2
         TRUE(1)=  0.1706167732170456E+01;
         TRUE(2)= -0.8928097010248257E+00;
    case 3
         TRUE(1)=  -0.1510606936744095E+01;
         TRUE(2)=   0.1178380000730945E+01;
    case 4
         TRUE(1)=   0.1194414677695905E+01;
         TRUE(2)=  -0.2799585996540082E+01;
    case 5
         TRUE(1)=   0.1890428596416747E+01;
         TRUE(2)=  -0.7345118680166940E+00;
    case 6
         TRUE(1)=-0.1737716306805883E+01;
         TRUE(2)=   0.8604008653025923E+00;
    case 7
         TRUE(1)=   0.1551614645548223E+01;
         TRUE(2)=  -0.1102382892533321E+01;
    case 8
         TRUE(1)=-0.1278631984330405E+01;
         TRUE(2)=   0.2013890883009155E+01;
    case 9
         TRUE(1)=  -0.1916552949489830E+01;
         TRUE(2)=   0.7169573003463228E+00;
    case 10
         TRUE(1)=   0.1768163792391936E+01;
         TRUE(2)=-0.8315276407898496E+00;
    case 11
         TRUE(1)=-0.1590150544829062E+01;
         TRUE(2)= 0.1040279389212485E+01;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,X,Y,DF,IPAR,RPAR,ierr] = F(N,X,Y,DF,IPAR,RPAR,ierr)
%F subroutine for the Van Der Pol Problem
%disp('now in F')
      
      EPS=1E-6;
      DF(1)=Y(2);
      PROD=1-Y(1)*Y(1);
      DF(2)=(PROD*Y(2)-Y(1))/EPS;
%disp('letf F')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MAS(n,am,m,IPAR,RPAR,ierr)

%DUMMY ROUTINE

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [X,Y,DFY,N,MEBAND,IPAR,RPAR,ierr] =  PDERV(X,Y,DFY,N,MEBAND,IPAR,RPAR,ierr)
%PDERV subroutine for the Van Der Pol Problem
%disp('now in pderv')

      EPS=1E-6;
      DFY(1,1)=0;
      DFY(1,2)=1;
      DFY(2,1)=(-2*Y(1)*Y(2)-1)/EPS;
      DFY(2,2)=(1-Y(1)^2)/EPS;

%disp('left pderv')
end