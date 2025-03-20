% s = tf('s');
% H = 5*exp(-120*s)/(180*s+1);
% t = 0:0.1:500;
% u = max(0,min(t-1,10));
% u=t>0*10;
% lsim(H,u,t);
% grid on;


clear; clc;

% col = 'bgmr';
% M = 10; % Step input magnitude
% tau = 3; % Time constant
% strt = 2; % Time when the step started
% K=5;
% y = K*M*(1-exp(-(t-strt)/tau)).*heaviside(t-strt);
% plot(y); 


% G = zpk([-1,-2],[-4, -.5, 2],10);
% step(G); grid on;

% G = zpk(2, [-1, -2], -2.5);
% step(G); grid on;

% kc=500; ti=0.1; td=0.1; f=250; v=10; p=50; cp=1;
% a=p*v/f;
% s=tf('s');
% x=1/(a*s+1);
% step(x); grid on; hold on;
% b=f*cp*ti/kc; c=ti*(p*cp*v/kc+td); d=ti*(f*cp/kc+1);
% y=b*s/(c*s^2+d*s+1);
% step(y); legend('without controller','with controller');

s=tf('s');
% x=(5*exp(-s))/(s^2+2*s+1);
% y=(5*exp(-s))/(s^2+3*s+1);
% step(x); hold on;
% step(y);
% xlim([0 10])

% x=(-3.5-(100000/(s+5.5)))/(s-10000/(s+5.5)-1999.38+.0576/(s+0.88))
% x=(-3.5-(100000/(s+5.5)))/(s+10000/(s+5.5)+2000.62-.0576/(s+0.88))
% 
% v=20; vc=5; deltah=-100000; ko=5*exp(20); er=10000; p=1000; pc=1000; cp=1; cpc=1; a=1200*sqrt(2); b=0.5;
% f=10; fc=2; cain=1.1; tin=430; tcin=300; caout=0.1; tout=500; tcout=450;
% 
% syms caout tout tcout;
% syms v vc deltah ko er p pc cp cpc a b f fc cain tin tcin;
% 
% diff(f/v*cain-f/v*caout-ko*exp(-er/tout)*caout, f)
%     cain/v - caout/v   f/v  - f/v - ko*exp(-er/tout)  -(caout*er*ko*exp(-er/tout))/tout^2
% 
% diff(f/v*tin-f/v*tout+ko/p/cp*exp(-er/tout)*caout*-deltah-a*fc^b/p/cp/v*((tout-tcin+tout-tcout)/2), f)
% 
% diff(fc/vc*tcin-fc/vc*tcout+a*fc^b/pc/cpc/vc*((tout-tcin+tout-tcout)/2), fc)
% (tcin - tcout)/vc  fc/vc  -fc/vc
% -(a*b*fc^(b - 1)*(tcin/2 + tcout/2 - tout))/(cpc*pc*vc)
% (a*fc^b)/(cpc*pc*vc)  -(a*fc^b)/(2*cpc*pc*vc)
% 
% x=(-3.5+(25/(s+5.5)))/(s+10/(s+5.5)-1.38-.0288/(s+0.64))
% 
% factor(z^3 + 4.76*z^2 + 5.018*z + 1.384, z, 'FactorMode', 'complex')
% [z + 3.4060302808987645879239478616422, z + 0.90495495131612061465452247261834, z + 0.44901476778511479742152966573941]
% (s+3.4060303) (s+0.90495495)(s+0.44901477)
% 
% s*(s+5.5)*(s+0.64)+10*(s+0.64)-1.38*(s+5.5)*(s+0.64)-0.0288*(s+5.5)
% 
% 
% x=((-3.5*s+5.75)*(s+0.64))/(s*(s+5.5)*(s+0.64)+10*(s+0.64)-1.38*(s+5.5)*(s+0.64)-0.0288*(s+5.5));
% step(x); hold on; grid on;
% g=-0.107693737172466;f =2.779611932698980; h = 2.658959537572255;
% y=h*(1-g*s)/(f*s+1);
% step(y); legend('Original TF','FOPTD Approx.');
% 
% syms z; factor((z*(z+5.5)*(z+0.64)+10*(z+0.64)-1.38*(z+5.5)*(z+0.64)-0.0288*(z+5.5)), z, 'FactorMode', 'complex')
% 
% 
% t = 0:0.01:10;
% y = -exp(-t)+0.5*exp(-2*t)+0.5-(-exp(-(t-5))+0.5*exp(-2*(t-5))+0.5).*heaviside(t-5);
% plot(t,y)
% y=1/((s+1)*(s+2))
% *(1/s)*(1-exp(-5*s))
% 
% sys1=1/(s^2+3*s+2); sys2=2/(2*s^2+3*s+2); sys3=3/(3*s^2+3*s+2);
% [u,t] = gensig("square",10);
% lsim(sys1, sys2, sys3, u,t)
% grid on