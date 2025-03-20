% s = tf('s');

clear; clc;

syms k;
g1=tf(2,[1 1]);
g2=0.5;
g3=tf(-1,[1 2]);
g4=5;
g5=tf(1,[1 1]);
controlsystemdesigner()

a=linspace(-2,13,10000)
Kc = linspace(-2,13,100);
G = tf(1,[10 7 Kc]); % Process TF
H = tf(1,[1 1]); % Measurement TF
rlocus(G*H); % rlocus(TF,KcValues)


Km=4; Kv=0.1737847223; tauv=1/6; A=3; taui=0.15;
rlocus(tf(1,poly([-3 -0.5 -2]))) 
sys=1+Km*Kc*(1+1/(taui*s))*-Kv/(tauv*s+1)*-1/(A*s);



s = tf('s');
sys=10/(-20*s+1)
step((-0.5*sys)/(1+sys*-0.5)) %no 1
step((-0.5*(1+1/(100*s))*sys)/(1+sys*-0.5*(1+1/(100*s)))) %no 2
syms s; factor((-8000*s^4-79800*s^3-10040*s^2-400.5*s-5), s, 'FactorMode', 'complex') % no 4


s = tf('s'); step((-12050*s^2-500.5*s-5)/(-8000*s^4-79800*s^3+2010*s^2+100*s-12050*s^2-500.5*s-5))
s = tf('s'); 
plot(simulink.Time, simulink.Data); hold on;
step((-12050*s^2-500.5*s-5)/(-8000*s^4-79800*s^3+2010*s^2+100*s-12050*s^2-500.5*s-5),'--r'); 
xlim([0 1])
legend({'Simulink','MATLAB'},'Location','southeast')

gp=10/(-20*s+1); gm=1/(40*s+1); gc=tf([-1205 -50.05 -0.5], [10 100 0]);
step((gm*gp*gc)/(1+gm*gp*gc))


limit(((-12050*x^2-500.5*x-5)/(-8000*x^4-79800*x^3+2010*x^2+100*x-12050*x^2-500.5*x-5)), x, 0)

% 
% K = linspace(-10,10,1000);
% x = arrayfun(@(K) roots([-4 (5*K-16) (15*K-16) (5*K-4)]),K,'Uni',0);
% x = sortrows(cell2mat(x))';
% xR = real(x); xI = imag(x);
% for j = 1:size(x,2) 
%     plot(xR(:,j),xI(:,j),'b','LineWidth',1.1); hold on;
% end
% title('Root-locus diagram');
% hold off; grid on;


% 
% Kc = linspace(-2,13,1000);
% y = arrayfun(@(Kc) roots([10 17 8 (1+Kc)]),Kc,'Uni',0);
% y = sortrows(cell2mat(y))';
% yR = real(y); yI = imag(y);
% for j = 1:size(y,2) 
% plot(yR(:,j),yI(:,j),'b',... 
% 'LineWidth',1.1); hold on;
% end
% hold off; grid on;


syms a; b=0.5;
factor(-4*a^3+a^2*(5*b-16)+a*(15*b-16)+(5*b-4), a, 'FactorMode', 'complex')

plot((-6.1557*exp(-9.8481*t))+(8.8603*exp(-0.0569*t))-(3.6515*exp(-0.0419*t)+((9*10^-4)*exp(-0.0254*t+1))))

syms s;
a1=(-0.1543)/(s+9.8481); b1=ilaplace(a1)
a2=0.3792/(s+0.0596); b2=ilaplace(a2)
a3=-0.2263/(s+0.0419); b3=ilaplace(a3)
a4=0.0014/(s+0.0254); b4=ilaplace(a4)

plot(-0.1543*exp(-9.8481*t)+0.3792*exp(-0.0596*t)-0.2263*exp(-0.0419*t)+0.0014*(-0.0254*t))

t=1:0.1:180;
plot(t, (0.0156665*exp(-9.84808*t) - 6.36113*exp(-0.05961*t) + 5.40083*exp(-0.0419027*t) - 0.0553592*exp(-0.0254079*t) + 1))
hold on;
plot(simulink.Time, simulink.Data, 'or')
legend({'MATLAB','Simulink'},'Location','southeast')

0.0156665 e^(-9.84808 t) - 6.36113 e^(-0.05961 t) + 5.40083 e^(-0.0419027 t) - 0.0553592 e^(-0.0254079 t) + 1

x=[-12050 -500.5 -5]; y=[-8000 -79800 -10040 -400.5 -5 0]; [r,p,k]=residue(x,y)


% K = linspace(-10,10,10000);
% x = arrayfun(@(K) roots([-4 -(5*K-16) -(15*K-16) -(5*K-4)]),K,'Uni',0);
% x = sortrows(cell2mat(x))';
% xR = real(x); xI = imag(x);
% for j = 1:size(x,2) 
%     plot(xR(:,j),xI(:,j),'b','LineWidth',1.1); hold on;
% end
% title('Root-locus diagram');
% hold off; grid on;
% 


x0=0; syms x;
fzero((-0.1543*exp(-9.8481*x)+0.3792*exp(-0.0596*x)-0.2263*exp(-0.0419*x)+0.0014*exp(-0.0254*x)), x0)

syms t;
x=diff (0.0156665*exp(-9.84808*t) - 6.36113*exp(-0.05961*t) + 5.40083*exp(-0.0419027*t) - 0.0553592*exp(-0.0254079*t) + 1);
x0=[0 180];
fzero(x, x0)

syms t;
x=-0.154285*exp(-9.84808*t) + 0.379187*exp(-0.05961*t) - 0.226309*exp(-0.0419027*t) + 0.00140656*exp(-0.0254079*t);
t=0:1:100; plot(t, x)
solve(x,t)

fzero(@(t) (-0.154285*exp(-9.84808*t) + 0.379187*exp(-0.05961*t) - 0.226309*exp(-0.0419027*t) + 0.00140656*exp(-0.0254079*t)), 50)

function x=dt
x=diff (0.0156665*exp(-9.84808*t) - 6.36113*exp(-0.05961*t) + 5.40083*exp(-0.0419027*t) - 0.0553592*exp(-0.0254079*t) + 1)
end
z=@dt; x0=0;
fzero(z, x0)