%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE HIRES PROBLEM. WITH F, MAS AND PDERV TO BE 
% SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function HIRES
global MYtagSTIFF
MYtagSTIFF = 1;
%MAXN = 400;
ND=8;
LWORK =(31+3*ND)*ND+3;
LIWORK = ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);

MBND = zeros(1,4);
TRUE  = zeros(1,ND);
ERROR = zeros(1,ND);
MASBND = zeros(1,4);
%MYY = zeros(1,20);
MYVAR = 0;
RPAR = 0;
IPAR = 0;
IERR = 0;

%...  DIMENSION OF THE SYSTEM
N=8;
IDID =1;
IERR = 0;



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




      
LOUT = fopen('HIR.txt','w');
PROBLEM = 'HIRES';
SOLVER = 'MEBDFAE';
 fprintf(LOUT,'\r COMPUTATIONAL STATISTICS OF THE HIRES PROBLEM USING MEBDF \r',PROBLEM,SOLVER);
 fprintf(LOUT,'\r Number of Equations : %g \r',ND);

%... ENDPOINT OF INTEGRATION
    
         XEND=421.8122;
         XOUT=XEND-100;
         X=0;
%...  REQUIRED TOLERANCE
        RTOL=1E-7;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES
       Y(1)=1;
      Y(2)=0;
      Y(3)=0;
      Y(4)=0;
      Y(5)=0;
      Y(6)=0;
      Y(7)=0;
      Y(8)=0.0057;
%...  SET DEFAULT VALUES
        MF=21;
        INDEX=1;
        MAXDER=7;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=200000;
        H=1E-8;
        fprintf(LOUT,'\rHstart = %E',H);
        MASBND(1)=0;

%...
%... CALL OF THE SUBROUTINE
%...

MYtag = 0;
 MYcase = 220;
 for I = 1:2
    
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
            
                XOUT=XOUT+100;
                MYY1(I) = Y(1);
                MYY2(I) = Y(2);
                MYY3(I) = Y(3);
                MYY4(I) = Y(4);
                MYY5(I) = Y(5);
                MYY6(I) = Y(6);
                MYY7(I) = Y(7);
                MYY8(I) = Y(8);
                XX(I)   = X(1);
                MYtag = 34;
          continue
                
        end

     end

 end




for K = 1:2
      MYYTRUE = SOLN(K);
      
      fprintf(LOUT,'\r \r APPROX Y1 =  %E', MYY1(K));
      fprintf(LOUT,'\r TRUE Y1 =  %E', MYYTRUE(1));
      
      fprintf(LOUT,'\r \r APPROX Y2 =  %E', MYY2(K));
      fprintf(LOUT,'\r TRUE Y2 =  %E', MYYTRUE(2));
      
      fprintf(LOUT,'\r \r APPROX Y3 =  %E', MYY3(K));
      fprintf(LOUT,'\r TRUE Y3 =  %E', MYYTRUE(3));
      
      fprintf(LOUT,'\r \r APPROX Y4 =  %E', MYY4(K));
      fprintf(LOUT,'\r TRUE Y4 =  %E', MYYTRUE(4));

            fprintf(LOUT,'\r \r APPROX Y5 =  %E', MYY5(K));
      fprintf(LOUT,'\r TRUE Y5 =  %E', MYYTRUE(5));
      
      fprintf(LOUT,'\r \r APPROX Y6 =  %E', MYY6(K));
      fprintf(LOUT,'\r TRUE Y6 =  %E', MYYTRUE(6));
      
      fprintf(LOUT,'\r \r APPROX Y7 =  %E', MYY7(K));
      fprintf(LOUT,'\r TRUE Y7 =  %E', MYYTRUE(7));
      
      fprintf(LOUT,'\r \r APPROX Y8 =  %E', MYY8(K));
      fprintf(LOUT,'\r TRUE Y8 =  %E', MYYTRUE(8));
         
      fprintf(LOUT,'\r \r AT X =  %g', XX(K));
end

      fclose(LOUT)
      
     
end

function TRUE = SOLN(I)
      switch I
          case 1
         TRUE(1)=0.000737131257332567;
         TRUE(2)=0.000144248572631618;
         TRUE(3)=0.000058887297409676;
         TRUE(4)=0.001175651343283149;
         TRUE(5)=0.002386356198831330;
         TRUE(6)=0.006238968252742796;
         TRUE(7)=0.002849998395185769;
         TRUE(8)=0.002850001604814231;
          case 2
         TRUE(1)=0.000670305503581864;
         TRUE(2)=0.000130996846986347;
         TRUE(3)=0.000046862231597733;
         TRUE(4)=0.001044668020551705;
         TRUE(5)=0.000594883830951485;
         TRUE(6)=0.001399628833942774;
         TRUE(7)=0.001014492757718480;
         TRUE(8)=0.004685507242281520;
      end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,X,Y,DY,IPAR,RPAR,ierr] = F(N,X,Y,DY,IPAR,RPAR,ierr)
%F subroutine for the Hires problem
%disp('now in F')
      
      DY (1) = -1.71*Y(1) + 0.43*Y(2) + 8.32*Y(3) + 0.0007;
      DY (2) = 1.71*Y(1) - 8.75*Y(2);
      DY (3) = -10.03*Y(3) + 0.43*Y(4) + 0.035*Y(5);
      DY (4) = 8.32*Y(2) + 1.71*Y(3) - 1.12*Y(4);
      DY (5) = -1.745*Y(5) + 0.43*Y(6) + 0.43*Y(7);
      DY (6) = -280*Y(6)*Y(8) + 0.69*Y(4) + 1.71*Y(5) - 0.43D0*Y(6) + 0.69D0*Y(7);
      DY (7) = 280*Y(6)*Y(8) - 1.81*Y(7);
      DY (8) = -DY (7);

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
%PDERV subroutine for the Hires problem
%disp('now in pderv')
      DFY(1,1)=  -1.71;
      DFY(1,2)=   0.43;
      DFY(1,3)=   + 8.32;

      DFY(2,1)=   1.71;
      DFY(2,2)=   - 8.75;

      DFY(3,3)=   -10.03;
      DFY(3,4)=   0.43;
      DFY(3,5)=   + 0.035;

      DFY(4,2)=    8.32;
      DFY(4,3)=   + 1.71;
      DFY(4,4)=   - 1.12;

      DFY(5,5)=   -1.745;
      DFY(5,6)=   + 0.43;
      DFY(5,7)=   + 0.43;

      DFY(6,4)=    + 0.69;
      DFY(6,5)=    + 1.71;
      DFY(6,6)=    - 0.43 -280*Y(8);
      DFY(6,7)=    + 0.69;
      DFY(6,8)=     -280*Y(6);
 
      DFY(7,6)=     280*Y(8);
      DFY(7,7)=      - 1.81;
      DFY(7,8)=     280*Y(6);

      DFY(8,6)=    - 280*Y(8);
      DFY(8,7)=       1.81;
      DFY(8,8)=    - 280*Y(6);


%disp('left pderv')
end
