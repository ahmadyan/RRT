%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE WATER TUBE PROBLEM: INDEX = 2.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function WATER
global MYtagSTIFF 

MYtagSTIFF = 1;

ND=49;
NEQN = ND;
N = 49;
RTOL = zeros(1,N);
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
      
LOUT = fopen('water.txt','w');
fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE WATER PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r Number of Equations : %g \r',ND);



%... ENDPOINT OF INTEGRATION
         X=0;
         XEND=12*3600;
         XOUT=XEND;
         
%...  REQUIRED TOLERANCE

        %RTOL=.316227766016837972E-06;
        RTOL(1:N)=1E-7;
        ATOL=RTOL;
        ATOL(37:N)=1E6*RTOL(1);
        
        ITOL=5;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL(1));
        fprintf(LOUT,'\r ATOL = %E and %E',ATOL(1),ATOL(37));
       
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
         IWORK(1)=38;
         IWORK(2)=11;
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
   
      y(  1) =  0.2298488296477430d-002;
      y(  2) =  0.1188984650746585d-002;
      y(  3) =  0.1109503645730845d-002;
      y(  4) =  0.1589620100314825d-003;
      y(  5) =  0.1030022640715102d-002;
      y(  6) =  0.8710606306836165d-003;
      y(  7) =  0.3243571480903489d-002;
      y(  8) =  0.1109503645730845d-002;
      y(  9) =  0.7120986206521341d-003;
      y( 10) =  0.6414613963833099d-003;
      y( 11) =  0.9416978549524347d-003;
      y( 12) =  0.3403428519096511d-002;
      y( 13) =  0.2397639310739395d-002;
      y( 14) =  0.2397639310739395d-002;
      y( 15) =  0.3348581430454180d-002;
      y( 16) =  0.1353560017035444d-002;
      y( 17) =  0.1995021413418736d-002;
      y( 18) =  0.5746220741193575d-002;
      y( 19) =  0.4751940452918529d-001;
      y( 20) =  0.4751940452918529d-001;
      y( 21) =  0.4751940452918529d-001;
      y( 22) =  0.4751940452918529d-001;
      y( 23) =  0.4751940452918529d-001;
      y( 24) =  0.4751940452918529d-001;
      y( 25) =  0.4311196778792902d-001;
      y( 26) =  0.4751940452918529d-001;
      y( 27) =  0.4751940452918529d-001;
      y( 28) =  0.4751940452918529d-001;
      y( 29) =  0.4751940452918529d-001;
      y( 30) =  0.4249217433601160d-001;
      y( 31) =  0.4732336439609648d-001;
      y( 32) =  0.4732336439609648d-001;
      y( 33) =  0.4270002118868241d-001;
      y( 34) =  0.4751940452918529d-001;
      y( 35) =  0.4751940452918529d-001;
      y( 36) =  0.3651427026675656d-001;
      y( 37) =  0.1111268591478108d+006;
      y( 38) =  0.1111270045592387d+006;
      y( 39) =  0.1111271078730254d+006;
      y( 40) =  0.1111269851929858d+006;
      y( 41) =  0.1111269255355337d+006;
      y( 42) =  0.1111269322658045d+006;
      y( 43) =  0.1111269221703983d+006;
      y( 44) =  0.1111270121140691d+006;
      y( 45) =  0.1111274419515807d+006;
      y( 46) =  0.1111255158881087d+006;
      y( 47) =  0.1111278793439227d+006;
      y( 48) =  0.1111270995171642d+006;
      y( 49) =  0.1111298338971779d+006;
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [neqn,T,y,yprime] = INIT(neqn,T,y,yprime)
            

     %for i=1:neqn
     %    y(i) = 0d0;         
     %    yprime(i) = 0d0;
     %end
     y(1:neqn) = 0;
     yprime(1:neqn) = 0;
     
     %for i=19:36
     %    y(i) = 0.47519404529185289807d-1;
     %end
     y(19:36) = 0.47519404529185289807d-1;
          
     %for i=37:49
     %    y(i) = 109800d0;
     %end
     
     y(37:49) = 109800d0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [neqn,t,y,df,ipar,rpar,ierr] = F(neqn,t,y,df,ipar,rpar,ierr)
%disp('now in F')
   
nnodes=13;
nu    = 1.31d-6;
g     = 9.8d0;
rho   = 1.0d3;
rcrit = 2.3d3;
length= 1.0d3;
k     = 2.0d-4;
d     = 1.0d0;
b     = 2.0d2;
%pi    = 3.141592653589793238462643383d0;                  

      a = pi*d^2/4d0;
      mu = nu*rho;

      for j=1:nnodes
         ein(j) = 0d0;
         eout(j) = 0d0;
      end
      
      that=t/3600d0;
      that2=that*that;

      ein(1) = (1d0-cos(exp(-that)-1d0))/200d0;
      ein(13) = (1d0-cos(exp(-that)-1d0))/80d0;
      eout(10) = that2*(3d0*that2-92d0*that+720d0)/1d6;

      for j=1:nnodes
         for i=1:nnodes
            phi(i,j) = 0d0;
            lambda(i,j) = 1d0;
         end
      end

      phi( 1, 2) = y( 1);
      phi( 2, 3) = y( 2);
      phi( 2, 6) = y( 3);
      phi( 3, 4) = y( 4);
      phi( 3, 5) = y( 5);
      phi( 4, 5) = y( 6);
      phi( 5,10) = y( 7);
      phi( 6, 5) = y( 8);
      phi( 7, 4) = y( 9);
      phi( 7, 8) = y(10);
      phi( 8, 5) = y(11);
      phi( 8,10) = y(12);
      phi( 9, 8) = y(13);
      phi(11, 9) = y(14);
      phi(11,12) = y(15);
      phi(12, 7) = y(16);
      phi(12, 8) = y(17);
      phi(13,11) = y(18);

      lambda( 1, 2) = y(19);
      lambda( 2, 3) = y(20);
      lambda( 2, 6) = y(21);
      lambda( 3, 4) = y(22);
      lambda( 3, 5) = y(23);
      lambda( 4, 5) = y(24);
      lambda( 5,10) = y(25);
      lambda( 6, 5) = y(26);
      lambda( 7, 4) = y(27);
      lambda( 7, 8) = y(28);
      lambda( 8, 5) = y(29);
      lambda( 8,10) = y(30);
      lambda( 9, 8) = y(31);
      lambda(11, 9) = y(32);
      lambda(11,12) = y(33);
      lambda(12, 7) = y(34);
      lambda(12, 8) = y(35);
      lambda(13,11) = y(36);

      p( 5) = y(37);
      p( 8) = y(38);
      p( 1) = y(39);
      p( 2) = y(40);
      p( 3) = y(41);
      p( 4) = y(42);
      p( 6) = y(43);
      p( 7) = y(44);
      p( 9) = y(45);
      p(10) = y(46);
      p(11) = y(47);
      p(12) = y(48);
      p(13) = y(49);

      for j=1:nnodes
         for i=1:nnodes
            if (lambda(i,j)<0d0)
               ierr = -1;
               return
            end
            
            rtla=sqrt(lambda(i,j));
            r = abs(phi(i,j)*d/(nu*a));
            
            if (r>rcrit)
               rghres(i,j) = 1/rtla - 1.74d0 +2d0*log10(2d0*k/d + 18.7d0/(r*rtla));
               fdba(i,j) = p(i) - p(j) -lambda(i,j)*rho*length*phi(i,j)^2/(a^2*d);
            else
               rghres(i,j) = 1.d0/rtla - 1.74d0 +2d0*log10(2d0*k/d + 18.7d0/(rcrit*rtla));

               fdba(i,j) = p(i) - p(j) -32d0*mu*length*phi(i,j)/(a*d^2);
            end
            
         end
      end

      for nj=1:nnodes
         netflo(nj) = ein(nj)-eout(nj);
         for i=1:nnodes
            netflo(nj) = netflo(nj)+phi(i,nj);
         end
         
         
         for j=1:nnodes
            netflo(nj) = netflo(nj)-phi(nj,j);
         end
      end
      
      df( 1) = fdba( 1, 2);
      df( 2) = fdba( 2, 3);
      df( 3) = fdba( 2, 6);
      df( 4) = fdba( 3, 4);
      df( 5) = fdba( 3, 5);
      df( 6) = fdba( 4, 5);
      df( 7) = fdba( 5,10);
      df( 8) = fdba( 6, 5);
      df( 9) = fdba( 7, 4);
      df(10) = fdba( 7, 8);
      df(11) = fdba( 8, 5);
      df(12) = fdba( 8,10);
      df(13) = fdba( 9, 8);
      df(14) = fdba(11, 9);
      df(15) = fdba(11,12);
      df(16) = fdba(12, 7);
      df(17) = fdba(12, 8);
      df(18) = fdba(13,11);

      df(19) = rghres( 1, 2);
      df(20) = rghres( 2, 3);
      df(21) = rghres( 2, 6);
      df(22) = rghres( 3, 4);
      df(23) = rghres( 3, 5);
      df(24) = rghres( 4, 5);
      df(25) = rghres( 5,10);
      df(26) = rghres( 6, 5);
      df(27) = rghres( 7, 4);
      df(28) = rghres( 7, 8);
      df(29) = rghres( 8, 5);
      df(30) = rghres( 8,10);
      df(31) = rghres( 9, 8);
      df(32) = rghres(11, 9);
      df(33) = rghres(11,12);
      df(34) = rghres(12, 7);
      df(35) = rghres(12, 8);
      df(36) = rghres(13,11);

      df(37) = netflo( 5);
      df(38) = netflo( 8);
      df(39) = netflo( 1);
      df(40) = netflo( 2);
      df(41) = netflo( 3);
      df(42) = netflo( 4);
      df(43) = netflo( 6);
      df(44) = netflo( 7);
      df(45) = netflo( 9);
      df(46) = netflo(10);
      df(47) = netflo(11);
      df(48) = netflo(12);
      df(49) = netflo(13);
      
      

%disp('letf F')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [neqn,am,ldim,ipar,rpar,ierr] = MAS(neqn,am,ldim,ipar,rpar,ierr)
%disp('IN MAS')
  
nu    = 1.31d-6;
g     = 9.8d0;
rho   = 1.0d3;
rcrit = 2.3d3;
length= 1.0d3;
k     = 2.0d-4;
d     = 1.0d0;
b     = 2.0d2;
%pi    = 3.141592653589793238462643383d0 )  
         %3.1415926535897932384626433832795

      a = pi*d^2/4d0;
      c = b/(rho*g);
      v = rho*length/a;

      %for i=1:18
      %   am(1,i) = v;
      %end
      am(1,1:18) = v;
      
      %for i=19:36
      %   am(1,i) = 0d0;
      %end
      am(1,19:36) = 0;
      
      
      
      %for i=37:38
      %    am(1,i) = c;
      %end
      am(1,37:38) = c;
      

      %for i=39:neqn
      %   am(1,i) = 0d0;
      %end
       am(1,39:neqn) = 0d0;
     
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

