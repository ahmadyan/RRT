function [tn,y,pd,n,meband,ipar,rpar,ier] = PDERV(tn,y,pd,n,meband,ipar,rpar,ier)
%disp('now in pderv')
global variation0
global variation1

c=1.6E-8;cs=1E-9;cp=1E-8;r=25E3;rp=50;
lh=4.45;ls1=2E-3;ls2=5E-4;ls3=5E-4;
rg1=36.3;rg2=17.3;rg3=17.3;ri=5E1;rc=6E2;
gamma=40.67286402E-9;delta=17.7493332;


%       for i=1:n
%          do 177 j=1,n
%             pd(i,j)=0.0D+0
%       end      

pd = zeros(n,n);



	  uin2=2.0E0*sin(2.0E4*pi*tn)   + variation1 ;  
	  ud1=+y(3)-y(5)-y(7)-uin2;
	  ud2=-y(4)+y(6)-y(7)-uin2;
	  ud3=+y(4)+y(5)+y(7)+uin2;
	  ud4=-y(3)-y(6)+y(7)+uin2;
	  qpud1=gamma*delta*exp(delta*ud1);
	  qpud2=gamma*delta*exp(delta*ud2);
	  qpud3=gamma*delta*exp(delta*ud3);
	  qpud4=gamma*delta*exp(delta*ud4);
	  pd(1,1)=-1.0E0/(c*r);
	  pd(1,8)=-1.0E0/c;
	  pd(1,10)=-0.5E0/c;
	  pd(1,11)=-pd(1,10);
	  pd(1,14)=pd(1,8);
	  pd(2,2)=pd(1,1);
	  pd(2,9)=pd(1,8);
	  pd(2,12)=pd(1,10);
	  pd(2,13)=pd(1,11);
	  pd(2,15)=pd(1,14);
	  pd(3,3)=(-qpud1-qpud4)/cs;
	  pd(3,5)=qpud1/cs;
	  pd(3,6)=-qpud4/cs;
	  pd(3,7)=(qpud1+qpud4)/cs;
	  pd(3,10)=1.0E0/cs;
	  pd(4,4)=(-qpud2-qpud3)/cs;
	  pd(4,5)=-qpud3/cs;
	  pd(4,6)=qpud2/cs;
	  pd(4,7)=(-qpud2-qpud3)/cs;
	  pd(4,11)=-1.0E0/cs;
	  pd(5,3)=qpud1/cs;
	  pd(5,4)=-qpud3/cs;
	  pd(5,5)=(-qpud1-qpud3)/cs;
	  pd(5,7)=(-qpud1-qpud3)/cs;
	  pd(5,12)=1.0E0/cs;
	  pd(6,3)=-qpud4/cs;
	  pd(6,4)=qpud2/cs;
	  pd(6,6)=(-qpud2-qpud4)/cs;
	  pd(6,7)=(qpud2+qpud4)/cs;
	  pd(6,13)=-1.0E0/cs;
	  pd(7,3)=(qpud1+qpud4)/cp;
	  pd(7,4)=(-qpud2-qpud3)/cp;
	  pd(7,5)=(-qpud1-qpud3)/cp;
	  pd(7,6)=(qpud2+qpud4)/cp;
	  pd(7,7)=(-qpud1-qpud2-qpud3-qpud4-1.0E0/rp)/cp;
	  pd(8,1)=-1.0E0/lh;
	  pd(9,2)=pd(8,1);
	  pd(10,1)=0.5E0/ls2;
	  pd(10,3)=-1.0E0/ls2;
	  pd(10,10)=-rg2/ls2;
	  pd(11,1)=-0.5E0/ls3;
	  pd(11,4)=1.0E0/ls3;
	  pd(11,11)=-rg3/ls3;
	  pd(12,2)=pd(10,1);
	  pd(12,5)=pd(10,3);
	  pd(12,12)=pd(10,10);
	  pd(13,2)=pd(11,1);
	  pd(13,6)=pd(11,4);
	  pd(13,13)=pd(11,11);
	  pd(14,1)=-1.0E0/ls1;
	  pd(14,14)=-(ri+rg1)/ls1;
	  pd(15,2)=pd(14,1);
	  pd(15,15)=-(rc+rg1)/ls1;
%disp('left pderv')
end

