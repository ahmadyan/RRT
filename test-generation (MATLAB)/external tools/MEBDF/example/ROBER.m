%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE ROBERTSON PROBLEM. WITH F, MAS AND PDERV TO BE 
% SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function ROBER
global MYtagSTIFF
MYtagSTIFF = 1;
%MAXN = 400;
ND=3;
LWORK =(32+3*ND)*ND+3;
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
N=3;
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
      
LOUT = fopen('ROB.txt','w');

fprintf(LOUT,'\r COMPUTATIONAL STATISTICS OF THE ROBERTSON PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r NUMBER OF EQUAITONS : %g \r',ND);

%... ENDPOINT OF INTEGRATION
         X=0;
         XOUT=1;
         XEND=1E11;
%...  REQUIRED TOLERANCE
        RTOL=1E-13;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES
        Y(1)=1;
        Y(2)=0;
        Y(3)=0;
%...  SET DEFAULT VALUES
        MF=21;
        INDEX=1;
        MAXDER=7;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=200000;
        H=1E-8;
        fprintf(LOUT,'\rHstart = %E',H)
        MASBND(1)=0;

%...
%... CALL OF THE SUBROUTINE
%...

MYtag = 0;
 MYcase = 220;
 for I = 1:12
    
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
            
                XOUT=XOUT*10;
                MYY1(I) = Y(1);
                MYY2(I) = Y(2);
                MYY3(I) = Y(3);
                XX(I)   = X(1);
                
                MYtag = 34;
          continue
                
        end

     end

 end


for K = 1:12
      MYYTRUE = SOLN(K);
      
      fprintf(LOUT,'\r \r APPROX Y1 =  %E', MYY1(K));
      fprintf(LOUT,'\r TRUE Y1 =  %E', MYYTRUE(1));
      
      fprintf(LOUT,'\r \r APPROX Y2 =  %E', MYY2(K));
      fprintf(LOUT,'\r TRUE Y2 =  %E', MYYTRUE(2));
      
      fprintf(LOUT,'\r \r APPROX Y3 =  %E', MYY3(K));
      fprintf(LOUT,'\r TRUE Y3 =  %E', MYYTRUE(3));
      fprintf(LOUT,'\r \r AT X =  %E', XX(K));
end

      fclose(LOUT)
      
     
end

function TRUE = SOLN(I)
      switch I
          case 1
         TRUE(1)= 0.9664597373330035E+00;
         TRUE(2)=0.3074626578578675E-04;
         TRUE(3)=0.3350951640121071E-01;
          case 2
         TRUE(1)=0.8413699238414729E+00;
         TRUE(2)=0.1623390937990473E-04;
         TRUE(3)=0.1586138422491472E+00;
          case 3
         TRUE(1)=0.6172348823960878E+00;
         TRUE(2)=0.6153591274639123E-05;
         TRUE(3)=0.3827589640126376E+00;
          case 4
         TRUE(1)=0.3368745306607069E+00;
         TRUE(2)=0.2013702318261393E-05;
         TRUE(3)=0.6631234556369748E+00;
          case 5
         TRUE(1)=0.1073004285378040E+00;
         TRUE(2)=0.4800166972571660E-06;
         TRUE(3)=0.8926990914454987E+00;
          case 6
         TRUE(1)=0.1786592114209946E-01;
         TRUE(2)=0.7274751468436319E-07;
         TRUE(3)=0.9821340061103859E+00;
          case 7
         TRUE(1)= 0.2031483924973415E-02;
         TRUE(2)=0.8142277783356159E-08;
         TRUE(3)=0.9979685079327488E+00;
          case 8
         TRUE(1)=0.2076093439016395E-03;
         TRUE(2)=0.8306077485067610E-09;
         TRUE(3)=0.9997923898254906E+00;
          case 9
         TRUE(1)=0.2082417512179460E-04;
         TRUE(2)=0.8329841429908955E-10;
         TRUE(3)=0.9999791757415798E+00;
          case 10
         TRUE(1)=0.2083229471647004E-05;
         TRUE(2)=0.8332935037760723E-11;
         TRUE(3)=0.9999979167621954E+00;
          case 11
         TRUE(1)=0.2083328471883087E-06;
         TRUE(2)=0.8333315602809495E-12;
         TRUE(3)=0.9999997916663195E+00;
          case 12
         TRUE(1)=0.2083340149701284E-07;
         TRUE(2)=0.8333360770334744E-13;
         TRUE(3)=0.9999999791665152E+00;
      end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,X,Y,DF,IPAR,RPAR,ierr] = F(N,X,Y,DF,IPAR,RPAR,ierr)
%F subroutine for the Robertson problem
%disp('now in F')
      
      DF(1)=-0.04*Y(1)+1.E4*Y(2)*Y(3);
      DF(3)=3.E7*Y(2)*Y(2);
      DF(2)=-DF(1)-DF(3);
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
%PDERV subroutine for the Robertson problem.
%disp('now in pderv')

         PROD1=1.0E4*Y(2);
      PROD2=1.0E4*Y(3);
      PROD3=6.0E7*Y(2);
      DFY(1,1)=-0.04;
      DFY(1,2)=PROD2;
      DFY(1,3)=PROD1;
      DFY(2,1)=0.04;
      DFY(2,2)=-PROD2-PROD3;
      DFY(2,3)=-PROD1;
      DFY(3,1)=0;
      DFY(3,2)=PROD3;
      DFY(3,3)=0;

%disp('left pderv')
end