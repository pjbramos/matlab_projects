clc
clear

% Input the following variables:
temp = input('Enter system temperature (in K): ');
n = input('Enter no. of species in the mixture: ');
m = input('Enter no. of subgroups present in the species of the mixture: ');                 
% Number assignment of subgroups is based on the UNIFAC-VLE Subgroup Parameter list on SVAS. 
% e.g. if the system have species with CH3, CH2, and CH2NH subgroups, CH3 will be assigned with 1, CH2 with 2 and CH2NH with 3.

% Initialize variables
x = zeros(1,n);
vk = zeros(m,n);
a_mn = zeros(m,m);
Rk = zeros(1,m);
Qk = zeros(1,m);

for b = 1 : n
    x(1,b) = input(['Enter liquid composition (x) of species ' num2str(b) ': ']);
end

for i = 1 : m
    for j = 1 : n
        vk(i,j) = input(['Enter no. of subgroup ' num2str(i) ' present in species ' num2str(j) '- '])
    end
end

for f = 1 : m
    for g = 1 : m
        a_mn(f,g) = input(['Enter main group interaction parameter value between subgroup ' num2str(f) ' and subgroup ' num2str(g) '- '])
    end
end

for i = 1 : m
    Rk(1,i) = input(['Enter Rk value for subgroup ' num2str(i) ': ']);
end

for j = 1 : m
    Qk(1,j) = input(['Enter Qk value for subgroup ' num2str(j) ': ']);
end

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

fprintf(['The activity coefficients of the mixture are: ']);
fprintf('\n');
for a = 1 : n
    fprintf('gamma%d:   %.4f\n', a, gamma(1,a));
end




