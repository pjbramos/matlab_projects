clc
clear

%input
xA = [0:0.005:1];
xB = 1 - xA;
xC = cat(1,xA,xB);
xD = xC.';
global gamma;
temp = 308.15;  %in K
pressure = 100; %in kPa

%P-xy data calculations
global b_P;
global y;
[b_P, y] = bubblePcalcs(temp, xD, gamma);

%T-xy data calculations
global b_T;
global y_t;
[b_T, y_t] = bubbleTcalcs(pressure, xD, gamma);

%P-xy data from Part 2 of LPA 2
Pb = [19.625; 37.877; 45.476; 47.969; 49.489; 50.316; 50.969; 51.302; 51.409; 50.196; 45.863];
xl = [0; 0.098; 0.198; 0.283; 0.387; 0.489; 0.598; 0.7; 0.771; 0.895; 1];
yl = [0; 0.482; 0.601; 0.635; 0.665; 0.686; 0.709; 0.732; 0.754; 0.841; 1];

%P-xy plot
%to plot the P-x,y data, enter the following in the command window.
%h = figure;
%set(h,'Color',[1 1 1]);
%plot(xA,b_P,'b',y(:,1),b_P,'r')
%xlabel('x1, y1')
%ylabel('P(in kPa)')
%title('P-x,y diagram of (1)Acetone/(2)Cyclohexane system at T = 308.15K')
%xlim([0 1])
%to overlap the given P-x,y data, enter the succeeding command therafter.
%hold on
%scatter(x1, Pb)
%scatter(y1, Pb)
%hold off

%T-xy plot
%to plot the T-x,y data, enter the following in the command window.
%h = figure;
%set(h,'Color',[1 1 1]);
%plot(xA,b_T,'b',y_t(:,1),b_T,'r')
%xlabel('x1, y1')
%ylabel('T(in K)')
%title('T-x,y diagram of (1)Acetone/(2)Cyclohexane system at P = 1 bar')
%xlim([0 1])

function [b_P,y] = bubblePcalcs(temp, xD, gamma)
%This function computes for the bubble pressure and vapor compositions
%for an (1)Acetone /(2)Cyclohexane system at temperature temp(K) at 
%various liquid fractions of Acetone based on the Modified Raoult's Law.
    %initialize variable and matrix holders
    P1sat = 0;
    P2sat = 0;
    Pbub = zeros(1,numel(xD(:,1)));
    y = zeros(numel(xD(:,1)),2);
    %precalculations
    A1 = 14.3145;
    A2 = 13.6568;
    B1 = 2756.22;
    B2 = 2723.44;
    C1 = 228.060;
    C2 = 220.618;
    %calculate the saturated pressure (Psat(kPa)) of each species at the given
    %temp. using the Antoine Equation
    P1sat = exp(A1-(B1/((temp-273.15)+C1)));
    P2sat = exp(A2-(B2/((temp-273.15)+C2)));
    %calculate bubble pressure and vapor composition of the mixture
        %calculate the bubble pressure, Pbub
        Pbub = gamma(:,1).*xD(:,1).*P1sat + gamma(:,2).*xD(:,2).*P2sat;
        %calculate the vapor phase composition of each species, y(i)
        y1 = (gamma(:,1).*xD(:,1).*P1sat)./(Pbub);
        y2 = (gamma(:,2).*xD(:,2).*P2sat)./(Pbub);
        y = [y1 y2];
        b_P = Pbub;
end

function [b_T, y_t] = bubbleTcalcs(pressure, xD, gamma)
%This function computes for the bubble temperature and vapor compositions
%for an (1)Acetone /(2)Cyclohexane system at pressure pressure(kPa) at 
%various liquid fractions of Acetone based on the Modified Raoult's Law.
    %antoine constants
    A1 = 14.3145;
    A2 = 13.6568;
    B1 = 2756.22;
    B2 = 2723.44;
    C1 = 228.060;
    C2 = 220.618;
    b_T=zeros(size(gamma,1), 1);
    y_t=zeros(size(gamma,1), 1);
    syms T
    for i=1:length(gamma)
        y1=xD(i,1)*exp(antoine(A1,B1,C1,T))*gamma(i,1)/pressure;
        y2=xD(i,2)*exp(antoine(A1,B1,C1,T))*gamma(i,2)/pressure;
        eqn1=y1+y2==1;
        temp=solve(eqn1,T);
        b_T(i,1)=273.15+temp;
        y_t(i,1)=xD(i,1)*exp(antoine(A1,B1,C1,temp))*gamma(i,1)/pressure;
    end
end

function [lnP] = antoine(a,b,c,T)
%This function computes for the ln of the pressure given the antoine
%coefficients and the unknown temperature T
lnP=a-(b/(T+c));
end