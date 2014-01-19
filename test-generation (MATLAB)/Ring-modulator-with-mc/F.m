function [n,tn,y,dy,ipar,rpar,ier] = F(n,tn,y,dy,ipar,rpar,ier)
%disp('now in F')
global variation0
global variation1

c=1.6E-8;
cs=1E-9;
cp=1E-8;
r=25E3;
rp=50;
lh=4.45;
ls1=2E-3;
ls2=5E-4;
ls3=5E-4;
rg1=36.3;
rg2=17.3;
rg3=17.3;
ri=5E1;
rc=6E2;
gamma=40.67286402E-9;
delta=17.7493332;


	  uin1=0.5E0*sin(2.0E3*pi*tn) + variation0;
	  uin2=2.0E0*sin(2.0E4*pi*tn) + variation1;
	  ud1=+y(3)-y(5)-y(7)-uin2;
	  ud2=-y(4)+y(6)-y(7)-uin2;
	  ud3=+y(4)+y(5)+y(7)+uin2;
	  ud4=-y(3)-y(6)+y(7)+uin2;
	  qud1=gamma*(exp(delta*ud1)-1.0E0);
	  qud2=gamma*(exp(delta*ud2)-1.0E0);
	  qud3=gamma*(exp(delta*ud3)-1.0E0);
	  qud4=gamma*(exp(delta*ud4)-1.0E0);

	  dy(1)=(y(8)-0.5E0*y(10)+0.5E0*y(11)+y(14)-y(1)/r)/c;
	  dy(2)=(y(9)-0.5E0*y(12)+0.5E0*y(13)+y(15)-y(2)/r)/c;
	  dy(3)=(y(10)-qud1+qud4)/cs;
	  dy(4)=(-y(11)+qud2-qud3)/cs;
	  dy(5)=(y(12)+qud1-qud3)/cs;
	  dy(6)=(-y(13)-qud2+qud4)/cs;
	  dy(7)=(-y(7)/rp+qud1+qud2-qud3-qud4)/cp;
	  dy(8)=-y(1)/lh;
	  dy(9)=-y(2)/lh;
	  dy(10)=(0.5E0*y(1)-y(3)-rg2*y(10))/ls2;
	  dy(11)=(-0.5E0*y(1)+y(4)-rg3*y(11))/ls3;
	  dy(12)=(0.5E0*y(2)-y(5)-rg2*y(12))/ls2;
	  dy(13)=(-0.5E0*y(2)+y(6)-rg3*y(13))/ls3;
	  dy(14)=(-y(1)+uin1-(ri+rg1)*y(14))/ls1;
	  dy(15)=(-y(2)-(rc+rg1)*y(15))/ls1;


%disp('letf F')
end