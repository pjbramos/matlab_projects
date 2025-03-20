clc;
clear;

%input    
xA = [0:0.005:1];
xB = 1 - xA;
x = cat(1,xA,xB);
x(:,1) = [];
x(:,200) = [];
x;
temp = 308.15;

%activity coefficient calculations
global gamma;
for row1 = 1 : numel(x(1,:))
    gammafinal(row1,:) = gamma_unifac(temp, x(:,row1));
end
gamma01 = [0 1];
gamma02 = [1 0];
gamma = [gamma01; gammafinal; gamma02];

function gamma = gamma_unifac(temp, x)
%This function computes for the activity coefficients for an (1)Acetone
%/(2)Cyclohexane system at temperature temp(K) at various liquid fractions
%of Acetone based on the UNIFAC model.
    %precalculations
    n = 2;
    m = 3;
    vk = [1 2; 0 4; 1 0];                           
    a_mn = [0 0 476.4; 0 0 476.4; 26.7 26.7 0];
    Rk = [0.9011 0.6744 1.6724];
    Qk = [0.8480 0.5400 1.4880];

    %Calculate r and q
    for i = 1 : n
        r(i) = 0;    %Initialize
        q(i) = 0;    %Initialize
        for k = 1 : m
            r(i) = r(i) + vk(k,i)*Rk(k);
            q(i) = q(i) + vk(k,i)*Qk(k);
        end
    end

    %Calculate J and L
    denomJ = 0;     %Initialize
    denomL = 0;     %Initialize
    for a = 1 : n
        denomJ = denomJ + r(a)*x(a);
        denomL = denomL + q(a)*x(a);
    end

    for b = 1 : n
        J(b) = r(b)/denomJ;
        L(b) = q(b)/denomL;
    end

    %Calculate e
    for k = 1 : m
        for j = 1 : n
            e(k,j) = vk(k,j)*Qk(k)/q(j);
        end
    end

    %Calculate tau
    for c = 1 : m
        for d = 1 : m
            tau(c,d) = exp(-(a_mn(c,d))/temp);
        end
    end

    %Calculate beta
    for a = 1 : n
        for b = 1 : m
            beta(a,b) = 0;
            for l = 1 : m
                beta(a,b) = beta(a,b) + e(l,a)*tau(l,b);
            end
        end
    end

    %Calculate theta
    for y = 1 : m
        theta(y) = 0;
        for z = 1 : n
            theta(y) = theta(y) + (x(z)*q(z)*e(y,z));
        end
        theta(y) = theta(y)/denomL;
    end

    %Calculate s
    for t = 1 : m
        s(t) = 0;
        for u = 1 : m
            s(t) = s(t) + (theta(u)*tau(u,t));
        end
    end

    %Calculate gamma
    for i = 1 : n
        lngammaC(i) = 1 - J(i) + log(J(i)) - 5*q(i)*((1 - (J(i)/L(i)) + (log(J(i)/L(i)))));
        lngammaR(i) = 0;
        for j = 1 : m
            lngammaR(i) = lngammaR(i) + (theta(j)*(beta(i,j)/s(j))) - (e(j,i)*(log(beta(i,j)/s(j))));
        end
        lngammaR(i) = q(i)*(1 - lngammaR(i));
        gamma(i) = exp(lngammaC(i) + lngammaR(i));
    end
end

