%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE FEKETE PROBLEM: INDEX = 2.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


function FEKETE
global MYtagSTIFF 

MYtagSTIFF = 1;

ND=160;
NART = 20;
NEQN = ND;
N = 160;

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
      
LOUT = fopen('fekete.txt','w');


%... ENDPOINT OF INTEGRATION
         X=0;
         XEND=1E+3;
         XOUT=XEND;
         
%...  REQUIRED TOLERANCE

        RTOL=1E-13;
        ATOL=RTOL;
        ITOL=2;
        fprintf(LOUT,'RESULT WITH THE FOLLOWING TOL :');
        fprintf(LOUT,'\r RTOL = %E',RTOL);
        fprintf(LOUT,'\r ATOL = %E',ATOL);
       
%...  INITIAL VALUES
        [N,X,Y,YPRIME] = INIT(N,X,Y,YPRIME);
%...  SET DEFAULT VALUES
        MF=21;
        INDEX=1;
        MAXDER=7;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=200000;
        H=1E-7;
         IWORK(1)=6*NART;
         IWORK(2)=2*NART;
         IWORK(3)=0;

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

function y = SOLN(dummy)
      y(  1) =  -0.4070263380333202E+00;
      y(  2) =   0.3463758772791802E+00;
      y(  3) =   0.8451942450030429E+00;
      y(  4) =   0.7752934752521549E-01;
      y(  5) =  -0.2628662719972299E+00;
      y(  6) =   0.9617122871829146E+00;
      y(  7) =   0.7100577833343567E+00;
      y(  8) =   0.1212948055586120E+00;
      y(  9) =   0.6936177005172217E+00;
      y( 10) =   0.2348267744557627E+00;
      y( 11) =   0.7449277976923311E+00;
      y( 12) =   0.6244509285956391E+00;
      y( 13) =  -0.4341114738782885E+00;
      y( 14) =   0.8785430442262876E+00;
      y( 15) =   0.1992720444237660E+00;
      y( 16) =  -0.9515059600312596E+00;
      y( 17) =   0.2203508762787005E+00;
      y( 18) =   0.2146669498274008E+00;
      y( 19) =  -0.6385191643609878E+00;
      y( 20) =  -0.4310833259390688E+00;
      y( 21) =   0.6375425027722121E+00;
      y( 22) =  -0.1464175087914336E+00;
      y( 23) =  -0.9380871635228862E+00;
      y( 24) =   0.3139337298744690E+00;
      y( 25) =   0.5666974065069942E+00;
      y( 26) =  -0.6739221885076542E+00;
      y( 27) =   0.4740073135462156E+00;
      y( 28) =   0.9843259538440293E+00;
      y( 29) =  -0.1696995357819996E+00;
      y( 30) =  -0.4800504290609090E-01;
      y( 31) =   0.1464175087914331E+00;
      y( 32) =   0.9380871635228875E+00;
      y( 33) =  -0.3139337298744656E+00;
      y( 34) =  -0.7092757549979014E+00;
      y( 35) =   0.5264062637139616E+00;
      y( 36) =  -0.4688542938854929E+00;
      y( 37) =  -0.8665731819284478E+00;
      y( 38) =  -0.4813878059756024E+00;
      y( 39) =  -0.1315929352982178E+00;
      y( 40) =  -0.2347897778700538E+00;
      y( 41) =  -0.8594340408013130E+00;
      y( 42) =  -0.4541441287957579E+00;
      y( 43) =   0.5530976940074118E+00;
      y( 44) =  -0.7674370265615124E+00;
      y( 45) =  -0.3242273140037833E+00;
      y( 46) =   0.7711050969896927E+00;
      y( 47) =   0.6357041816577034E+00;
      y( 48) =   0.3573685519777001E-01;
      y( 49) =   0.7103951209379591E+00;
      y( 50) =   0.2403570431280519E+00;
      y( 51) =  -0.6614886725910596E+00;
      y( 52) =  -0.3038208738735660E-01;
      y( 53) =   0.4501923293640461E+00;
      y( 54) =  -0.8924145871442046E+00;
      y( 55) =  -0.5772996158107093E+00;
      y( 56) =  -0.1766763414971813E+00;
      y( 57) =  -0.7971892020969544E+00;
      y( 58) =   0.2414481766969039E+00;
      y( 59) =  -0.3416456818373135E+00;
      y( 60) =  -0.9082846503446250E+00;
      y( 61) =   0.2409619682166627E-15;
      y( 62) =  -0.1139818460497816E-15;
      y( 63) =   0.1627536276556335E-15;
      y( 64) =   0.1745651819597609E-15;
      y( 65) =  -0.1914278710633076E-15;
      y( 66) =  -0.6639600671806291E-16;
      y( 67) =   0.1708576733899083E-15;
      y( 68) =  -0.2277602521390053E-15;
      y( 69) =  -0.1350782790950654E-15;
      y( 70) =   0.2411941341109454E-15;
      y( 71) =  -0.1438238671800488E-15;
      y( 72) =   0.8087033550666644E-16;
      y( 73) =   0.1618239105233347E-15;
      y( 74) =   0.1837556152070701E-16;
      y( 75) =   0.2715177369929503E-15;
      y( 76) =   0.7930078658689191E-16;
      y( 77) =   0.7482020588342764E-16;
      y( 78) =   0.2746974939098084E-15;
      y( 79) =   0.8849338913035911E-16;
      y( 80) =  -0.5940734725324115E-16;
      y( 81) =   0.4845984056889910E-16;
      y( 82) =  -0.3728835248155620E-16;
      y( 83) =  -0.4600332954062859E-16;
      y( 84) =  -0.1548568884846698E-15;
      y( 85) =   0.2507541692375411E-16;
      y( 86) =  -0.1560155223230823E-15;
      y( 87) =  -0.2517946296860555E-15;
      y( 88) =  -0.3739779361502470E-16;
      y( 89) =  -0.1381663620885020E-15;
      y( 90) =  -0.2784051540342329E-15;
      y( 91) =   0.6624397102887671E-16;
      y( 92) =   0.4226207488883120E-16;
      y( 93) =   0.1571821772296610E-15;
      y( 94) =  -0.4112243677286995E-16;
      y( 95) =   0.1939960344265876E-15;
      y( 96) =   0.2800184977692136E-15;
      y( 97) =  -0.9189023375328813E-16;
      y( 98) =   0.1392943179389155E-15;
      y( 99) =   0.9556003995587458E-16;
      y(100) =  -0.2234188557495892E-15;
      y(101) =   0.1276804778190781E-15;
      y(102) =  -0.1261196211463950E-15;
      y(103) =  -0.1887754149742397E-15;
      y(104) =  -0.2140788698695373E-16;
      y(105) =  -0.2713591291421657E-15;
      y(106) =   0.1107887633060814E-15;
      y(107) =  -0.1318443715631340E-15;
      y(108) =  -0.4521275683078691E-16;
      y(109) =  -0.1277688851278605E-15;
      y(110) =   0.4850914012115388E-16;
      y(111) =  -0.1195891666741192E-15;
      y(112) =  -0.1569641653843750E-15;
      y(113) =   0.1856239009452638E-15;
      y(114) =   0.9898466095646496E-16;
      y(115) =  -0.2068030800303723E-15;
      y(116) =   0.2451470336752085E-15;
      y(117) =   0.9542986459336358E-16;
      y(118) =  -0.2456074075580993E-15;
      y(119) =   0.1532475480661800E-15;
      y(120) =  -0.1229326332276474E-15;
      y(121) =  -0.4750000000000000E+01;
      y(122) =  -0.4750000000000001E+01;
      y(123) =  -0.4750000000000000E+01;
      y(124) =  -0.4750000000000000E+01;
      y(125) =  -0.4750000000000000E+01;
      y(126) =  -0.4750000000000000E+01;
      y(127) =  -0.4750000000000000E+01;
      y(128) =  -0.4750000000000000E+01;
      y(129) =  -0.4750000000000000E+01;
      y(130) =  -0.4750000000000000E+01;
      y(131) =  -0.4750000000000001E+01;
      y(132) =  -0.4750000000000001E+01;
      y(133) =  -0.4750000000000000E+01;
      y(134) =  -0.4750000000000000E+01;
      y(135) =  -0.4750000000000000E+01;
      y(136) =  -0.4750000000000000E+01;
      y(137) =  -0.4749999999999999E+01;
      y(138) =  -0.4750000000000000E+01;
      y(139) =  -0.4750000000000000E+01;
      y(140) =  -0.4750000000000000E+01;
      y(141) =  -0.3537526598492654E-19;
      y(142) =   0.2338193888161182E-18;
      y(143) =  -0.3267771993164953E-18;
      y(144) =   0.2915679914072042E-18;
      y(145) =   0.1965183195887647E-18;
      y(146) =  -0.6224992924096233E-19;
      y(147) =  -0.1715878416756298E-18;
      y(148) =  -0.2704741705248803E-18;
      y(149) =   0.3008700893194513E-18;
      y(150) =  -0.2703121624910402E-18;
      y(151) =   0.4243755291982164E-18;
      y(152) =   0.2862063003949612E-18;
      y(153) =   0.1222125408406218E-19;
      y(154) =  -0.4958862706817728E-18;
      y(155) =  -0.7070673036251212E-18;
      y(156) =  -0.4454983024194383E-18;
      y(157) =  -0.1125384872521777E-18;
      y(158) =   0.1512898724592511E-18;
      y(159) =  -0.6163704221424137E-19;
      y(160) =   0.6255426995473074E-19;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [n,t,y,yprime] = INIT(n,t,y,yprime)
      
nart=20;
neqn=8*nart;

      %pi=3.141592653589793238462643383d0      
      nart= n/8;
      y(1) = 0d0;
      y(2) = 0d0;
      y(3) = 0d0;
      
      for i=1:3
         alpha=2*pi*i/3+pi/13;
         beta=3*pi/8;
         y(3*(i-1)+1)=cos(alpha)*cos(beta);
         y(3*(i-1)+2)=sin(alpha)*cos(beta);
         y(3*(i-1)+3)=sin(beta);
      end
      
      for i=4:10
         alpha=2*pi*(i-3)/7+pi/29;
         beta=pi/8;
         y(3*(i-1)+1)=cos(alpha)*cos(beta);
         y(3*(i-1)+2)=sin(alpha)*cos(beta);
         y(3*(i-1)+3)=sin(beta);
      end
      
      for i=11:16
         alpha=2*pi*(i-10)/(6)+pi/(7);
         beta=-2*pi/(15);
         y(3*(i-1)+1)=cos(alpha)*cos(beta);
         y(3*(i-1)+2)=sin(alpha)*cos(beta);
         y(3*(i-1)+3)=sin(beta);
      end
      
      for i=17:20
         alpha=2*pi*(i-17)/(4)+pi/(17);
         beta=-3*pi/(10);
         y(3*(i-1)+1)=cos(alpha)*cos(beta);
         y(3*(i-1)+2)=sin(alpha)*cos(beta);
         y(3*(i-1)+3)=sin(beta);
      end
      
      for i=(3*nart+1):(6*nart)
         y(i)=0d0;
      end
      
      for i=(6*nart+1):(8*nart)
         y(i)=0d0;
      end
      
      [n,t,y,y,yprime]=F_EVAL(n,t,y,y,yprime);
      
      for i=1:nart
         for j=1:3
            y(6*nart+i)=y(6*nart+i)+y(3*(i-1)+j)*...
                yprime(3*nart+3*(i-1)+j);
         end
         y(6*nart+i)=-y(6*nart+i)/2d0;
      end
      
      [n,t,y,y,yprime] =F_EVAL(n,t,y,y,yprime);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [n,t,y,df,ipar,rpar,ierr] = F(n,t,y,df,ipar,rpar,ierr)
 %disp('now in F')

  [n,t,y,y,df] = F_EVAL(n,t,y,y,df);

 %disp('letf F')
end

function [neqn,t,y,yprime,dy] =F_EVAL(neqn,t,y,yprime,dy)      
 
maxn=150;
alpha=.5d0;
nart=neqn/8;
    
for i=1:nart
    for k=1:3
        p(i,k)=y(3*(i-1)+k);
        q(i,k)=y(3*nart+3*(i-1)+k);
    end
    lam(i)=y(6*nart+i);
    mu(i)=y(7*nart+i);
end

for i=1:nart
    for j=1:nart
        if(i==j)
            for k=1:3
                fa(i,j,k)=0d0;
            end

        else
            rn=0d0;
            for k=1:3
                rn=rn+(p(i,k)-p(j,k))^2;
            end
            
            for k=1:3
                fa(i,j,k)=(p(i,k)-p(j,k))/rn;
            end
            
        end
        
    end
end
    
for i=1:nart
    for k=1:3
        pp(i,k)=q(i,k)+2*mu(i)*p(i,k);
        qp(i,k)=-alpha*q(i,k)+2*lam(i)*p(i,k);
        
        for j=1:nart
            qp(i,k)=qp(i,k)+fa(i,j,k);
        end
    end
end

for i=1:nart
    phi(i)=-1d0;
    gpq(i)=0d0;
    
    for k=1:3
        phi(i)=phi(i)+p(i,k)^2;
        gpq(i)=gpq(i)+2*p(i,k)*q(i,k);
    end
end

for i=1:nart
    for k=1:3
        dy(3*(i-1)+k)=pp(i,k);
        dy(3*nart+3*(i-1)+k)=qp(i,k);
    end
        dy(6*nart+i)=phi(i);
        dy(7*nart+i)=gpq(i);
end
      
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [NEQN,AM,LDIM,IPAR,RPAR,IERR] = MAS(NEQN,AM,LDIM,IPAR,RPAR,IERR)
%disp('IN MAS')

      NART=NEQN/8;      
      for J =1:(6*NART)
         AM(1,J) = 1;         
      end
      
      for J =(6*NART+1):(8*NART)
         AM(1,J) = 0;         
      end
  
%disp('LEFT MAS')

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [t,y,dfy,neqn,yprime,ipar,rpar,ierr] =...
    PDERV(t,y,dfy,neqn,yprime,ipar,rpar,ierr)      
%disp('in pderv')                                            

maxn=150;
alpha=.5d0;
nart=neqn/8;
      
for i=1:nart
    for k=1:3
        p(i,k)=y(3*(i-1)+k);
        q(i,k)=y(3*nart+3*(i-1)+k);
    end
        lam(i)=y(6*nart+i);
        mu(i)=y(7*nart+i);
end

      
for j=1:nart
    for i=1:nart
        rn(i,j)=0d0;
        for k=1:3
            rn(i,j)=rn(i,j)+(p(i,k)-p(j,k))^2;
        end
    end
end
     
for j=1:neqn
    for i=1:neqn
        dfy(i,j) =0d0;
    end
end

%J_pp
 for i=1:nart
    for k=1:3
        dfy(3*(i-1)+k,3*(i-1)+k)=2d0*mu(i);
    end
 end
 
%J_pq
for i=1:nart
    for k=1:3
        dfy(3*(i-1)+k,3*nart+3*(i-1)+k)=1d0;
    end
end

%J_pm
 for i=1:nart
    for k=1:3
        dfy(3*(i-1)+k,7*nart+i)=2d0*p(i,k);
    end
 end
      
%J_qp
for i=1:nart
    for k=1:3
        dfy(3*nart+3*(i-1)+k,3*(i-1)+k)=2d0*lam(i);
        for j=1:nart
            if(j~=i)
            dfy(3*nart+3*(i-1)+k,3*(i-1)+k)=...
                dfy(3*nart+3*(i-1)+k,3*(i-1)+k)+...
                (rn(i,j)-2d0*(p(i,k)-p(j,k))^2)/rn(i,j)^2;
            end
        end
    end
end

%J_qp
for i=1:nart
    for k=1:3
        for m=1:3
            if(m~=k)
                for j=1:nart
                    if(j~=i)
                        dfy(3*nart+3*(i-1)+k,3*(i-1)+m)=...
                            dfy(3*nart+3*(i-1)+k,3*(i-1)+m)-...
                            2d0*(p(i,k)-p(j,k))*(p(i,m)-p(j,m))/...
                            rn(i,j)^2;
                    end
                end
            end
        end
    end
end

%J_qp
for i=1:nart
    for l=1:nart
        if(l~=i)
            for k=1:3
                dfy(3*nart+3*(i-1)+k,3*(l-1)+k)=...
                    (-rn(i,l)+2d0*(p(i,k)-p(l,k))^2)/rn(i,l)^2;
            end
        end
    end
end

%J_qp

for i=1:nart
    for l=1:nart
        if(l~=i)
            for k=1:3
                for m=1:3
                    if(m~=k)
                        dfy(3*nart+3*(i-1)+k,3*(l-1)+m)=...
                            2d0*(p(i,k)-p(l,k))*(p(i,m)-p(l,m))/...
                            rn(i,l)^2;
                    end
                end
            end
        end
    end
end


%J_qq
for i=1:nart
    for k=1:3
        dfy(3*nart+3*(i-1)+k,3*nart+3*(i-1)+k)=-alpha;
    end
end

%J_ql
for i=1:nart
    for k=1:3
        dfy(3*nart+3*(i-1)+k,6*nart+i)=2d0*p(i,k);
    end
end


%J_lp
 for i=1:nart
    for k=1:3
        dfy(6*nart+i,3*(i-1)+k)=2d0*p(i,k);
    end
 end
 

%J_mp
for i=1:nart
    for k=1:3
        dfy(7*nart+i,3*(i-1)+k)=2d0*q(i,k);
    end
end

%J_mq
 for i=1:nart
    for k=1:3
        dfy(7*nart+i,3*nart+3*(i-1)+k)=2d0*p(i,k);
    end
 end
 %disp('left pderv')
 
 end


