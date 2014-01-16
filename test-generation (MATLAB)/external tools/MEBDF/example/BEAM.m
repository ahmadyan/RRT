%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE ELASTIC BEAM PROBLEM PROBLEM. WITH F, MAS AND  
% PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function BEAM

global MYtagSTIFF MY_COUNTER
global NNNN NCOM NNCOM NSQ NQUATR DELTAS N
MYtagSTIFF = 1;
ND=80;
LWORK=(32+3*ND)*ND+3;
LIWORK=ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);
YPRIME = zeros(1,ND);
N = 40;
MASBND = zeros(1,4);
MBND = zeros(1,4);
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
      
LOUT = fopen('BEAM.txt','w');

 fprintf(LOUT,'\r COMPUTATIONAL STATISTICS OF THE BEAM PROBLEM USING MEBDF \r');
 fprintf(LOUT,'\r Number of Equations : %g \r',ND);

%... ENDPOINT OF INTEGRATION
    
         XEND=5;
         XOUT=XEND;
         X=0;
%...  REQUIRED TOLERANCE

        RTOL=1E-7;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %f',RTOL);
        fprintf(LOUT,'\r ATOL = %f',ATOL);
       
%...  INITIAL VALUES
      Y = zeros(1,ND);
         NN=2*N;
         NCOM=N;
         NSQ=N*N;
         NQUATR=NSQ*NSQ;
         NNCOM=NN;
         AN=N;
         DELTAS=1/AN;

%...  SET DEFAULT VALUES
        MF=22;
        INDEX=1;
        MAXDER=3;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=200000;
        H=1E-8;
        fprintf(LOUT,'Hstart = %E',H);
        MASBND(1) = 0;
%...
%... CALL OF THE SUBROUTINE
%...

MYtag = 0;

while (MYtag == 0)                    
           [NN,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,WORK_3,...
               WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,WORK_I7,...
               WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,IWORK_15,...
               MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR] =...
               MEBDF(NN,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,...
               WORK_2,WORK_3,WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,...
               WORK_I6,WORK_I7,WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,...
               IWORK,IWORK_15,MBND,MASBND,MAXDER,ITOL,RTOL,...
               ATOL,RPAR,IPAR,IERR);  

           if (INDEX==1)
              INDEX = 0;
              
              continue
           else
             
              MYtag =34;
              
           end
           
end


      TRUE(1)=-0.005792366591294675;
      TRUE(2)=-0.016952985507199259;
      TRUE(3)=-0.027691033129713322;
      TRUE(4)=-0.038008156558781729;
      TRUE(5)=-0.047906168597422688;
      TRUE(6)=-0.057387104352737008;
      TRUE(7)=-0.066453273134522699;
      TRUE(8)=-0.075107305819780661;
      TRUE(9)=-0.083352197654124544;
      TRUE(10)=-0.091191346546446469;
      TRUE(11)=-0.098628587001297248;
      TRUE(12)=-0.105668220037774708;
      TRUE(13)=-0.112315039540924422;
      TRUE(14)=-0.118574355272698475;
      TRUE(15)=-0.124452012875526880;
      TRUE(16)=-0.129954411326390999;
      TRUE(17)=-0.135088518061004200;
      TRUE(18)=-0.139861881919410397;
      TRUE(19)=-0.144282644101482929;
      TRUE(20)=-0.148359547246256976;
      TRUE(21)=-0.152101942900106414;
      TRUE(22)=-0.155519797806080921;
      TRUE(23)=-0.158623699341992299;
      TRUE(24)=-0.161424860370167541;
      TRUE(25)=-0.163935123819275499;
      TRUE(26)=-0.166166967344037066;
      TRUE(27)=-0.168133508177817718;
      TRUE(28)=-0.169848508060189926;
      TRUE(29)=-0.171326378244038509;
      TRUE(30)=-0.172582184746215274;
      TRUE(31)=-0.173631653797526901;
      TRUE(32)=-0.174491177383960691;
      TRUE(33)=-0.175177818786287100;
      TRUE(34)=-0.175709317871242317;
      TRUE(35)=-0.176104096022807288;
      TRUE(36)=-0.176381260717507812;
      TRUE(37)=-0.176560609756417469;
      TRUE(38)=-0.176662635226010517;
      TRUE(39)=-0.176708527080694206;
      TRUE(40)=-0.176720176107510191;
      TRUE(41)= 0.037473626808570053;
      TRUE(42)=0.109911788012810762;
      TRUE(43)=0.179836047447039129;
      TRUE(44)=0.247242730557127186;
      TRUE(45)=0.312129382035491301;
      TRUE(46)=0.374494737701689822;
      TRUE(47)=0.434338607372647125;
      TRUE(48)=0.491662035432760524;
      TRUE(49)=0.546467785483476383;
      TRUE(50)=0.598760970245279030;
      TRUE(51)=0.648549361126755851;
      TRUE(52)=0.695843516905088648;
      TRUE(53)=0.740657266848912124;
      TRUE(54)=0.783008174791347177;
      TRUE(55)=0.822917665884869456;
      TRUE(56)=0.860411030561688098;
      TRUE(57)=0.895517550233742218;
      TRUE(58)=0.928270826293034365;
      TRUE(59)=0.958708933474210358;
      TRUE(60)=0.986874782150222219;
      TRUE(61)=1.012816579967983789;
      TRUE(62)=1.036587736684594479;
      TRUE(63)=1.058246826485315033;
      TRUE(64)=1.077857811432700289;
      TRUE(65)=1.095490221995530989;
      TRUE(66)=1.111219164319120026;
      TRUE(67)=1.125125269269998022;
      TRUE(68)=1.137294526582397119;
      TRUE(69)=1.147818025203744592;
      TRUE(70)=1.156792131966898566;
      TRUE(71)=1.164318845152484938;
      TRUE(72)=1.170505992580311363;
      TRUE(73)=1.175467424328008220;
      TRUE(74)=1.179323003206967714;
      TRUE(75)=1.182198586301326345;
      TRUE(76)=1.184226111211404704;
      TRUE(77)=1.185543909813440450;
      TRUE(78)=1.186297084230907673;
      TRUE(79)=1.186637618874913665;
      TRUE(80)=1.186724615129383839;



for t = 1:NN
    
         fprintf(LOUT,'\r\r APROX Y%g = %E',t,Y(t));
         fprintf(LOUT,'\rTRUE Y%g = %E\r\r',t,TRUE(t));
end      

      fclose(LOUT)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NN,T,TH,DF,IPAR,RPAR,IERR] = F(NN,T,TH,DF,IPAR,RPAR,IERR)
%disp('now in F')
     
    

global NNNN NCOM NNCOM NSQ NQUATR DELTAS N
      
      
%C...  CALCUL DES TH(I) ET DES SIN ET COS -------------
    
      for I=2:N
         THDIFF=TH(I)-TH(I-1);
         STH(I)=sin(THDIFF);
         CTH(I)=cos(THDIFF);
      end
      %STH(2:N)=sin(THDIFF);
      %CTH(2:N)=cos(THDIFF);
      
%C...  CALCUL DU COTE DROIT DU SYSTEME LINEAIRE -----
 

if (T > pi)
    
% C --------- CASE T GREATER PI ------------
% C ---------- I=1 ------------
  

     TERM1=(-3*TH(1)+TH(2))*NQUATR;
     V(1)=TERM1;
         
%C -------- I=2,..,N-1 -----------

     for I=2:N-1
        TERM1=(TH(I-1)-2*TH(I)+TH(I+1))*NQUATR;
        V(I)=TERM1;
     end
%C ----------- I=N -------------
        
        TERM1=(TH(N-1)-TH(N))*NQUATR;
        V(N)=TERM1;
         
else
    
%C --------- CASE T LESS EQUAL PI ------------

         FABS=1.5*sin(T)*sin(T);
         FX=-FABS;
         FY= FABS;
         
%C ---------- I=1 ------------
  
         TERM1=(-3*TH(1)+TH(2))*NQUATR;
         TERM2=NSQ*(FY*cos(TH(1))-FX*sin(TH(1)));
         V(1)=TERM1+TERM2;
         
%C -------- I=2,..,N-1 -----------

         for I=2:N-1
            TERM1=(TH(I-1)-2*TH(I)+TH(I+1))*NQUATR;
            TERM2=NSQ*(FY*cos(TH(I))-FX*sin(TH(I)));
            V(I)=TERM1+TERM2;
         end
%C ----------- I=N -------------

         TERM1=(TH(N-1)-TH(N))*NQUATR;
         TERM2=NSQ*(FY*cos(TH(N))-FX*sin(TH(N)));
         V(N)=TERM1+TERM2;
         
end

%C -------- COMPUTE PRODUCT DV=W ------------
 
      W(1)=STH(2)*V(2);
      
      for I=2:N-1
        W(I)=-STH(I)*V(I-1)+STH(I+1)*V(I+1);
      end
      %W(2:N-1)=-STH(2:N-1)*V(1:N-2)+STH(3:N)*V(3:N);
      
      W(N)=-STH(N)*V(N-1);
      
%C -------- TERME 3 -----------------
 
       for I=1:N
         W(I)=W(I)+TH(N+I)*TH(N+I);
       end
       %W(1:N)=W(1:N)+TH(N+1:N+N)*TH(N+1:N+N);
%C ------- SOLVE SYSTEM CW=W ---------
      ALPHA(1)=1;
      
      for I=2:N
         ALPHA(I)=2;
         BETA(I-1)=-CTH(I);
      end
      %ALPHA(2:N)=2;
      %BETA(1:N-1)=-CTH(2:N);
      
      ALPHA(N)=3;
      
      for I=N-1:-1:1
          
         Q=BETA(I)/ALPHA(I+1);
         W(I)=W(I)-W(I+1)*Q;
         ALPHA(I)=ALPHA(I)-BETA(I)*Q;
      end
      
      W(1)=W(1)/ALPHA(1);
      
     for I=2:N
         W(I)=(W(I)-BETA(I-1)*W(I-1))/ALPHA(I);
     end
     %W(2:N)=(W(2:N)-BETA(1:N-1)*W(1:N-1))/ALPHA(2:N);
     
%C -------- COMPUTE U=CV+DW ---------

      U(1)=V(1)-CTH(2)*V(2)+STH(2)*W(2);
      
      for I=2:N-1
         U(I)=2*V(I)-CTH(I)*V(I-1)-CTH(I+1)*V(I+1)-STH(I)*W(I-1)+STH(I+1)*W(I+1);
      end
      %U(2:N-1)=2*V(2:N-1)-CTH(2:N-1)*V(1:N-2)-CTH(3:N)*V(3:N)-STH(2:N-1)*W(1:N-2)+STH(3:N)*W(3:N);
      
      U(N)=3*V(N)-CTH(N)*V(N-1)-STH(N)*W(N-1);
      
%C -------- PUT  DERIVATIVES IN RIGHT PLACE -------------

      for I=1:N
         DF(I)=TH(N+I);
         DF(N+I)=U(I);
      end
      %DF(1:N)=TH(N+1:N+N);
      %DF(N+1:N+N)=U(1:N);

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
%disp('now in pderv')

% DUMMY ROUTINE  

%disp('left pderv')
end


