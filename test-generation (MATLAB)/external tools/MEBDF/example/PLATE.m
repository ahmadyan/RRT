%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE MOVEMENT OF A RECTANGULAR PLATE PROBLEM PROBLEM.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function PLATE
global MYtagSTIFF
MYtagSTIFF = 1;
%MAXN = 400;
MX=8;MY=5;MACHS1=2;MACHS2=4;
ND=2*MX*MY;
LWORK =(32+3*ND)*ND+3;
LIWORK = ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);

global  TRANS NX NXM1 NY NYM1 NDEMI NACHS1 NACHS2 NDUMMY OMEGA STIFFN DELX USH4 FAC WEIGHT



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

      N=ND;
      NX=MX;
      NY=MY;
      NACHS1=MACHS1;
      NACHS2=MACHS2;
      NXM1=NX-1;
      NYM1=NY-1;
      NDEMI=NX*NY;
      OMEGA=1000;
      STIFFN=100;
      WEIGHT=200;
      DENOM=NX+1;
      DELX=2/DENOM;
      USH4=1/(DELX^4);
      FAC=STIFFN*USH4;



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




      
LOUT = fopen('plate.txt','w');

 fprintf(LOUT,'\r COMPUTATIONAL STATISTICS OF THE PLATE PROBLEM USING MEBDFAE \r');
 fprintf(LOUT,'\r Number of Equations : %g \r',ND);

%... ENDPOINT OF INTEGRATION
    
         XEND=7;
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
      Y = zeros(1,N);


%...  SET DEFAULT VALUES
        MF=21;
        INDEX=1;
        MAXDER=7;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=100000;
        H=1E-8;
        fprintf(LOUT,'\r Hstart = %E',H);
        MASBND(1)=0;
        

%...
%... CALL OF THE SUBROUTINE
%...

MYtag = 0

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


      TRUE(1)=0.000490143813851336;
      TRUE(2)=0.000980081485560611;
      TRUE(3)=0.001462893811482190;
      TRUE(4)=0.001915822464411935;
      TRUE(5)=0.002285152533727002;
      TRUE(6)=0.002461353376688549;
      TRUE(7)=0.002254597413097122;
      TRUE(8)=0.001438312591933600;
      TRUE(9)=0.000849025149228402;
      TRUE(10)=0.001697885005625757;
      TRUE(11)=0.002535239886068847;
      TRUE(12)=0.003323989552181772;
      TRUE(13)=0.003977902193560667;
      TRUE(14)=0.004320231736082990;
      TRUE(15)=0.004025679955083897;
      TRUE(16)=0.002643206356123840;
      TRUE(17)=0.000980287627702671;
      TRUE(18)=0.001960162971121222;
      TRUE(19)=0.002925787622964379;
      TRUE(20)=0.003831644928823870;
      TRUE(21)=0.004570305067454005;
      TRUE(22)=0.004922706753377098;
      TRUE(23)=0.004509194826194244;
      TRUE(24)=0.002876625183867201;
      TRUE(25)=0.000849025149228402;
      TRUE(26)= 0.001697885005625757;
      TRUE(27)=0.002535239886068847;
      TRUE(28)=0.003323989552181772;
      TRUE(29)=0.003977902193560667;
      TRUE(30)=0.004320231736082990;
      TRUE(31)=0.004025679955083897;
      TRUE(32)=0.002643206356123840;
      TRUE(33)=0.000490143813851336;
      TRUE(34)=0.000980081485560611;
      TRUE(35)=0.001462893811482190;
      TRUE(36)=0.001915822464411935;
      TRUE(37)=0.002285152533727002;
      TRUE(38)=0.002461353376688549;
      TRUE(39)=0.002254597413097122;
      TRUE(40)=0.001438312591933600;
      TRUE(41)=- 0.001177590304545409;
      TRUE(42)=-0.002409005827992214;
      TRUE(43)=-0.003722140831656533;
      TRUE(44)=-0.005078780056048207;
      TRUE(45)=-0.006302661811097914;
      TRUE(46)=-0.006973399942926759;
      TRUE(47)=-0.006394575120415784;
      TRUE(48)=-0.003960464551310118;
      TRUE(49)=-0.002040148244040460;
      TRUE(50)=-0.004174829877953482;
      TRUE(51)=-0.006456510337516159;
      TRUE(52)=-0.008832503276738242;
      TRUE(53)=-0.011029624807177369;
      TRUE(54)=-0.012352389570141255;
      TRUE(55)=-0.011524177328690540;
      TRUE(56)=-0.007253301886026949;
      TRUE(57)=-0.002355180609090818;
      TRUE(58)=-0.004818011655984428;
      TRUE(59)=-0.007444281663313065;
      TRUE(60)=-0.010157560112096413;
      TRUE(61)=-0.012605323622195828;
      TRUE(62)=-0.013946799885853517;
      TRUE(63)=-0.012789150240831569;
      TRUE(64)=-0.007920929102620235;
      TRUE(65)=-0.002040148244040460;
      TRUE(66)=-0.004174829877953482;
      TRUE(67)=-0.006456510337516159;
      TRUE(68)=-0.008832503276738242;
      TRUE(69)=-0.011029624807177369;
      TRUE(70)=-0.012352389570141255;
      TRUE(71)=-0.011524177328690540;
      TRUE(72)=-0.007253301886026949;
      TRUE(73)=-0.001177590304545409;
      TRUE(74)=-0.002409005827992214;
      TRUE(75)=-0.003722140831656533;
      TRUE(76)=-0.005078780056048207;
      TRUE(77)=-0.006302661811097914;
      TRUE(78)=-0.006973399942926759;
      TRUE(79)=-0.006394575120415784;
      TRUE(80)=-0.003960464551310118;

      
for t = 1:N
         fprintf(LOUT,'\r\r APROX Y%g = %E',t,Y(t));
         fprintf(LOUT,'\rTRUE Y%g = %E\r\r',t,TRUE(t));
end      

      fclose(LOUT)
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [N,X,Y,DF,IPAR,RPAR,ierr] = F(N,X,Y,DF,IPAR,RPAR,ierr)

%disp('now in F')
   
global  TRANS NX NXM1 NY NYM1...
    NDEMI NACHS1 NACHS2 NDUMMY ...
    OMEGA STIFFN DELX USH4 FAC WEIGHT
% C
% C -------- LA BOUCLE -------
      for I=1:NX
         for J=1:NY
            K=I+NX*(J-1);
% C -------- DERIVEE DEUXIEME ----
            DF(K)=Y(K+NDEMI);
% C ------ POINT CENTRAL ---
            UC=16*Y(K);
            if(I > 1)
               UC=UC+Y(K);
               UC=UC-8*Y(K-1);
            end
            
            if(I<NX)
               UC=UC+Y(K);
               UC=UC-8*Y(K+1);
            end
            
            if(J>1)
               UC=UC+Y(K);
               UC=UC-8*Y(K-NX);
            end
            
            if(J<NY)
               UC=UC+Y(K);
               UC=UC-8*Y(K+NX);
            end
            
            if(I>1 && J>1 )
            UC=UC+2*Y(K-NX-1);
            end
            
            if( I<NX && J>1 )
                UC=UC+2*Y(K-NX+1);
            end
            
            if(I > 1 && J < NY)
                UC=UC+2*Y(K+NX-1);
            end
            
            if(I<NX && J < NY)
            UC=UC+2*Y(K+NX+1);
            end
            
            if (I>2)
                UC=UC+Y(K-2);
            end
            
            if(I<NXM1)
                UC=UC+Y(K+2);
            end
            
            if(J>2)
                UC=UC+Y(K-2*NX);
            end
            
            if(J<NYM1)
                UC=UC+Y(K+2*NX);
            end
            
            if(J==NACHS1||J==NACHS2)
                
               XI=I*DELX;
               
               FORCE = exp(-5*(X-XI-2)^2)+ exp(-5*(X-XI-5)^2);
            else
               FORCE=0;
               
            end
            
            DF(K+NDEMI)=-OMEGA*Y(K+NDEMI)-FAC*UC+FORCE*WEIGHT;
            
         end
      end

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
 
      global  TRANS NX NXM1 NY NYM1 NDEMI NACHS1 NACHS2 NDUMMY OMEGA STIFFN DELX USH4 FAC WEIGHT
% C
% C------ METTRE A ZERO -------
      for I=1:N
         for J=1:N
           DFY(I,J)=0;
         end
      end
      
% C -------- LA BOUCLE -------

      for I=1:NX
         for J=1:NY
            K=I+NX*(J-1);
% C -------- DERIVEE DEUXIEME ----
            DFY(K,K+NDEMI)=1;
% C ------ POINT CENTRAL ---
            DFY(K+NDEMI,K)=-FAC*16;
            if(I>1)
               DFY(K+NDEMI,K)=DFY(K+NDEMI,K)-FAC;
               DFY(K+NDEMI,K-1)=FAC*8;
            end
            
            if(I<NX)
               DFY(K+NDEMI,K)=DFY(K+NDEMI,K)-FAC;
               DFY(K+NDEMI,K+1)=FAC*8;
            end
            
            if(J>1)
               DFY(K+NDEMI,K)=DFY(K+NDEMI,K)-FAC;
               DFY(K+NDEMI,K-NX)=FAC*8;
            end
            
            if(J<NY)
               DFY(K+NDEMI,K)=DFY(K+NDEMI,K)-FAC;
               DFY(K+NDEMI,K+NX)=FAC*8;
            end
            
            if(I>1 && J>1 )
                DFY(K+NDEMI,K-NX-1)=-FAC*2;
            end
            
            if(I<NX&&J>1 )
                DFY(K+NDEMI,K-NX+1)=-FAC*2;
            end
            
            if(I>1&&J<NY)
                DFY(K+NDEMI,K+NX-1)=-FAC*2;
            end
            
            if(I<NX&&J<NY)
                DFY(K+NDEMI,K+NX+1)=-FAC*2;
            end
            
            if(I>2)
                DFY(K+NDEMI,K-2)=-FAC;
            end
            
            if(I<NXM1)
                DFY(K+NDEMI,K+2)=-FAC;
            end
            
            if(J>2)
                DFY(K+NDEMI,K-2*NX)=-FAC;
            end
            
            if(J<NYM1)
                DFY(K+NDEMI,K+2*NX)=-FAC;
            end
                DFY(K+NDEMI,K+NDEMI)= -OMEGA;
            
         end
         
      end

%disp('left pderv')
end


