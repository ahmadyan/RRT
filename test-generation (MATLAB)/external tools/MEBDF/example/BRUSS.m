%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
% DRIVER ROUTINE FOR THE BRUSSELATOR PROBLEM PROBLEM.   
% WITH F, MAS AND PDERV TO BE SAVED AS !SEPARATE! M-FILES.
% AUGUST 2005.
%^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

function BRUSS
global MYtagSTIFF 
global N N2 GAMMA GAMMA2

MYtagSTIFF = 1;
ND=1000;
LWORK=(32+3*ND)*ND+3;
LIWORK=ND+14;
Y = zeros(1,ND);
IWORK = zeros(1,LIWORK);
YPRIME = zeros(1,ND);
N = 500;
N2 = 2*N;
MASBND = zeros(1,4);
MBND = zeros(1,4);
         MBND(1)=2;
         MBND(2)=2;
         MBND(3)=5;
         MBND(4)=7;

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
      
LOUT = fopen('bruss.txt','w');
fprintf(LOUT,'COMPUTATIONAL STATISTICS OF THE BRUSS PROBLEM USING MEBDF \r');
fprintf(LOUT,'\r Number of Equations : %g \r',ND);

%... ENDPOINT OF INTEGRATION
    
         XEND=10;
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
     USDELQ=(N+1)^2;
     GAMMA=0.02*USDELQ;
         GAMMA2=2*GAMMA;
         for I=1:N
            ANP1=N+1;
            XI=I/ANP1;
            Y(2*I)=3;
            Y(2*I-1)=1+0.5*sin(2*pi*XI);
         end


%...  SET DEFAULT VALUES
        MF=23;
        INDEX=1;
        MAXDER=7;
        IERR = 0;
%...  MAXIMAL NUMBER OF STEPS
%WORK(1) = 2;
        IWORK(14)=100000;
        H=1E-8;
        fprintf(LOUT,'Hstart = %E',H);
        MASBND(1) = 0;
%...
%... CALL OF THE SUBROUTINE
%...

MYtag = 0;

while (MYtag == 0)
          
          [N2,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,WORK_3,...
              WORK_I1,WORK_I2,WORK_I3,WORK_I4,WORK_I5,WORK_I6,WORK_I7,...
              WORK_I8,WORK_I9,WORK_I10,WORK_I11,LIWORK,IWORK,IWORK_15,...
              MBND,MASBND,MAXDER,ITOL,RTOL,ATOL,RPAR,IPAR,IERR] =MEBDF...
              (N2,X,H,Y,XOUT,XEND,MF,INDEX,LOUT,LWORK,WORK_1,WORK_2,...
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


      TRUE(1)=0.9949197002317599;
      TRUE(8)=3.0213845767604077;
      TRUE(15)=0.9594350193986054;
      TRUE(22)=3.0585989778165419;
      TRUE(29)=0.9243010095428502;
      TRUE(36)=3.0952478919989637;
      TRUE(43)=0.8897959106772672;
      TRUE(50)=3.1310118289054731;
      TRUE(57)=0.8561653620284367;
      TRUE(64)=3.1656101198770159;
      TRUE(71)=0.8236197147449046;
      TRUE(78)=3.1988043370624344;
      TRUE(85)=0.7923328094811884;
      TRUE(92)=3.2303999530641514;
      TRUE(99)=0.7624421042573115;
      TRUE(106)=3.2602463873623941;
      TRUE(113)=0.7340499750795348;
      TRUE(120)=3.2882356529108807;
      TRUE(127)=0.7072259700779899;
      TRUE(134)=3.3142998590079271;
      TRUE(141)=0.6820097782458483;
      TRUE(148)=3.3384078449410937;
      TRUE(155)=0.6584146743834650;
      TRUE(162)=3.3605612157873943;
      TRUE(169)=0.6364312187752559;
      TRUE(176)=3.3807900316323134;
      TRUE(183)=0.6160310186921587;
      TRUE(190)=3.3991483695914764;
      TRUE(197)=0.5971703941198909;
      TRUE(204)=3.4157099395342736;
      TRUE(211)=0.5797938277687891;
      TRUE(218)=3.4305638938070224;
      TRUE(225)=0.5638371159206763;
      TRUE(232)=3.4438109320334580;
      TRUE(239)=0.5492301695479158;
      TRUE(246)=3.4555597666485198;
      TRUE(253)=0.5358994429426996;
      TRUE(260)=3.4659239846027008;
      TRUE(267)=0.5237699892215797;
      TRUE(274)=3.4750193162238476;
      TRUE(281)=0.5127671585747183;
      TRUE(288)=3.4829613034792271;
      TRUE(295)=0.5028179665048467;
      TRUE(302)=3.4898633463634923;
      TRUE(309)=0.4938521662914935;
      TRUE(316)=3.4958350971335204;
      TRUE(323)=0.4858030633656755;
      TRUE(330)=3.5009811668111510;
      TRUE(337)=0.4786081100251151;
      TRUE(344)=3.5054001059792705;
      TRUE(351)=0.4722093177200750;
      TRUE(358)=3.5091836216744015;
      TRUE(365)=0.4665535216425440;
      TRUE(372)=3.5124159935026285;
      TRUE(379)=0.4615925290790646;
      TRUE(386)=3.5151736544621075;
      TRUE(393)=0.4572831793403656;
      TRUE(400)=3.5175249049438184;
      TRUE(407)=0.4535873393501199;
      TRUE(414)=3.5195297317024448;
      TRUE(421)=0.4504718553589467;
      TRUE(428)=3.5212397070273984;
      TRUE(435)=0.4479084778719241;
      TRUE(442)=3.5226979467564341;
      TRUE(449)=0.4458737738041973;
      TRUE(456)=3.5239391090719634;
      TRUE(463)=0.4443490371324889;
      TRUE(470)=3.5249894191569453;
      TRUE(477)=0.4433202068820853;
      TRUE(484)=3.5258667077466495;
      TRUE(491)=0.4427777991494095;
      TRUE(498)=3.5265804544017270;
      TRUE(505)=0.4427168579654424;
      TRUE(512)=3.5271318289682063;
      TRUE(519)=0.4431369281018266;
      TRUE(526)=3.5275137272135266;
      TRUE(533)=0.4440420513508381;
      TRUE(540)=3.5277107990730161;
      TRUE(547)=0.4454407863109616;
      TRUE(554)=3.5276994703501980;
      TRUE(561)=0.4473462502188303;
      TRUE(568)=3.5274479611304068;
      TRUE(575)=0.4497761798232572;
      TRUE(582)=3.5269163066324394;
      TRUE(589)=0.4527530066369863;
      TRUE(596)=3.5260563887768472;
      TRUE(603)=0.4563039400688689;
      TRUE(610)=3.5248119894251024;
      TRUE(617)=0.4604610498812091;
      TRUE(624)=3.5231188790654930;
      TRUE(631)=0.4652613370907894;
      TRUE(638)=3.5209049576992761;
      TRUE(645)=0.4707467798082714;
      TRUE(652)=3.5180904678044698;
      TRUE(659)=0.4769643375804777;
      TRUE(666)=3.5145883024867057;
      TRUE(673)=0.4839658945842979;
      TRUE(680)=3.5103044351908528;
      TRUE(687)=0.4918081185812277;
      TRUE(694)=3.5051385005173827;
      TRUE(701)=0.5005522089940899;
      TRUE(708)=3.4989845585737802;
      TRUE(715)=0.5102635039989190;
      TRUE(722)=3.4917320776245013;
      TRUE(729)=0.5210109134090777;
      TRUE(736)=3.4832671712209993;
      TRUE(743)=0.5328661417420937;
      TRUE(750)=3.4734741260299615;
      TRUE(757)=0.5459026646938675;
      TRUE(764)=3.4622372546582585;
      TRUE(771)=0.5601944229089820;
      TRUE(778)=3.4494431032230182;
      TRUE(785)=0.5758142001453760;
      TRUE(792)=3.4349830354873361;
      TRUE(799)=0.5928316594749734;
      TRUE(806)=3.4187562033108012;
      TRUE(813)=0.6113110218368440;
      TRUE(820)=3.4006728962523969;
      TRUE(827)=0.6313083867734524;
      TRUE(834)=3.3806582409098729;
      TRUE(841)=0.6528687160104193;
      TRUE(848)=3.3586561928427350;
      TRUE(855)=0.6760225267555723;
      TRUE(862)=3.3346337311179157;
      TRUE(869)=0.7007823726569721;
      TRUE(876)=3.3085851288057930;
      TRUE(883)=0.7271392249346637;
      TRUE(890)=3.2805361342349380;
      TRUE(897)=0.7550589020044152;
      TRUE(904)=3.2505478606008622;
      TRUE(911)=0.7844787296769868;
      TRUE(918)=3.2187201496972175;
      TRUE(925)=0.8153046416214843;
      TRUE(932)=3.1851941538893653;
      TRUE(939)=0.8474089465959840;
      TRUE(946)=3.1501538739882800;
      TRUE(953)=0.8806289904192589;
      TRUE(960)=3.1138264039027113;
      TRUE(967)=0.9147669230929857;
      TRUE(974)=3.0764806689389470;
      TRUE(981)=0.9495907429372025;
      TRUE(988)=3.0384245041548366;
      TRUE(995)=0.9848367306701233;



for t = 1:7:995
    
         fprintf(LOUT,'\r\r APROX Y%g = %E',t,Y(t));
         fprintf(LOUT,'\rTRUE Y%g = %E\r\r',t,TRUE(t));
end      

      fclose(LOUT)
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [NNN,X,Y,DF,ipar,rpar,ierr] = F(NNN,X,Y,DF,ipar,rpar,ierr)
%disp('now in F')
global N N2 GAMMA GAMMA2     
    
      I=1;
      IU=2*I-1;
      IV=2*I;
      UI=Y(IU);
      VI=Y(IV);
      UIM=1;
      VIM=3;
      UIP=Y(IU+2);
      VIP=Y(IV+2);
      PROD=UI*UI*VI;
      DF(IU)=1+PROD-4*UI+GAMMA*(UIM-2*UI+UIP);
      DF(IV)=3*UI-PROD+GAMMA*(VIM-2*VI+VIP);
      
      for I=2:N-1
         IU=2*I-1;
         IV=2*I;
         UI=Y(IU);
         VI=Y(IV);
         UIM=Y(IU-2);
         VIM=Y(IV-2);
         UIP=Y(IU+2);
         VIP=Y(IV+2);
         PROD=UI*UI*VI;
         DF(IU)=1+PROD-4*UI+GAMMA*(UIM-2*UI+UIP);
         DF(IV)=3*UI-PROD+GAMMA*(VIM-2*VI+VIP);
      end
      
      I=N;
      IU=2*I-1;
      IV=2*I;
      UI=Y(IU);
      VI=Y(IV);
      UIM=Y(IU-2);
      VIM=Y(IV-2);
      UIP=1;
      VIP=3;
      PROD=UI*UI*VI;
      DF(IU)=1+PROD-4*UI+GAMMA*(UIM-2*UI+UIP);
      DF(IV)=3*UI-PROD+GAMMA*(VIM-2*VI+VIP);
      
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

function [X,Y,DFY,NNN,mebnd,ipar,rpar,ierr] = PDERV(X,Y,DFY,NNN,mebnd,ipar,rpar,ierr)
%disp('now in pderv')
 
  
 global N N2 GAMMA GAMMA2   

      for I=1:N
         IU=2*I-1;
         IV=2*I;
         UI=Y(IU);
         VI=Y(IV);
         UIVI=UI*VI;
         UI2=UI*UI;
         DFY(3,IU)=2*UIVI-4-GAMMA2;
         DFY(2,IV)=UI2;
         DFY(4,IU)=3-2*UIVI;
         DFY(3,IV)=-UI2-GAMMA2;
         DFY(2,IU)=0;
         DFY(4,IV)=0;
      end
      
      for I=1:(N2-2)
         DFY(1,I+2)=GAMMA;
         DFY(5,I)=GAMMA;
      end


%disp('left pderv')
end



