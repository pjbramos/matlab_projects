% LPA1 Specification 2: Cubic Equation of State
% Peng-Robinson Equation of State (PR EoS)
% Multi-component Gas/Vapor System
%
% Definition:
%   y = vapor composition
%   press = total pressure of the system
%   temp = total temperature of the system
%   pressc = critical pressure of each species in the system
%   tempc = critical temperature of each species in the system
%   w = acentric factor
%   z = mixture compressibility factor
%   phi = fugacity coefficients of each species in the system
%   fhat = fugacity of each species in the system

clc
clear

% Input the following variables
  n = input('Enter no. of species in the mixture: ');
  for b = 1:n
    y(1,b)=input(['Enter vapor composition (y) of species ' num2str(b) ': ']);
  end
  temp = input('Enter Temperature (in K): ');
  press = input('Enter Pressure (in bar): ');
  fprintf('\n');
  tempc = zeros(1,n);
  pressc = zeros(1,n);
  w = zeros(1,n);
  for b = 1:n
      tempc(1,b) = input(['Enter critical temperature (Tc) of species ' num2str(b) ': ']);
      pressc(1,b) = input(['Enter critical pressure (Pc) of species ' num2str(b) ': ']);
      w(1,b) = input(['Enter acentric factor (w) of species ' num2str(b) ': ']);
  end

% Parameter constants for PR EoS
  psi = 0.45724;
  omega = 0.07780;
  sigma = 1 + sqrt(2);
  epsilon = 1 - sqrt(2);

% Determine the reduced temperature (tempr) for each species in the mixture
  tempr = zeros(1,n);
  for b = 1:n
      tempr(1,b) = temp/tempc(1,b);
  end
  
% Determine the parameter alpha in the PR EoS for each species in the system
  alpha_pr = zeros(1,n);
  for b = 1:n
      alpha_pr(1,b) = (1 + (0.37464 + 1.54226*w(1,b) - 0.26992*w(1,b)^2)*(1 - tempr(1,b)^(1/2)))^2;
  end

% Determine the values of the dimensionless parameters for pure species (a and b)
  % Solve for ai
    ai = zeros(1,n);
    for b = 1:n
        ai(1,b) = (psi*alpha_pr(1,b)*(83.14472^2)*(tempc(1,b)^2))/(pressc(1,b)); 
    end
  % Solve for bi
    bi = zeros(1,n);
    for b = 1:b
        bi(1,b) = (omega*83.14472*tempc(1,b))/(pressc(1,b));
    end
    
% Solve for the unlike interaction parameter (aij)
  % Initialize 
    count = 0;
    for b = 1:n
        count = count + b;
    end
  % Solve for aij
    d = 1;
    aij = zeros(1, count);
    for b = 1:n-1
        for c = b+1:n
            aij(1,d) = (ai(1,b)*ai(1,c))^(1/2);
            d = d + 1;
        end
    end
    for e = 1:n
        aij(1,d) = (ai(1,e)*ai(1,e))^(1/2);
        d = d + 1;
    end
  % Transform aij from an array to a matrix
    aijmat = zeros(n);
    d = 1;
    for b = 1:n-1
        for c = b+1:n
            aijmatrix(b,c) = aij(1,d);
            d = d + 1;
        end
    end
    d = 1;
    for b = 1:n-1
        for c = b+1:n
            aijmatrix(c,b) = aij(1,d);
            d = d + 1;
        end
    end
    for e = 1:n
        aijmatrix(e,e) = aij(1,d);
        d = d + 1;
    end
  
% Determine the values of the dimensionless parameters for mixtures (amix, bmix, beta, and q) 
  % Solve for amix
    amix = 0;
    for i = 1:n
        for j = 1:n
            amix = amix + (y(1,i)*y(1,j)*aijmatrix(i,j));
        end
    end
  % Solve for bmix
    bmix = 0;
    for i = 1:n
        bmix = bmix + (y(1,i)*bi(1,i));
    end
  % Solve for beta
    beta = 0;
    beta = (bmix*press)/(83.14472*temp);
  % Solve for q
    q = 0;
    q = (amix)/(bmix*83.14472*temp);
    
% Determine the mixture compressibility factor (z)
  % Use fzero function to find z
    z = fzero(@(z) (1 + beta - (q*beta)*((z - beta)/((z + epsilon*beta)*(z + sigma*beta))))-z, 0.5);
   
    fprintf('\n');
    fprintf('The value of the mixture compressibility factor (Z) is ');
    disp(z);
    fprintf('\n');
    
% Determine I in the PR EoS
  I = 0;
  I = (1/(sigma - epsilon))*log((z + sigma*beta)/(z + epsilon*beta));

% Determine the partial molar EoS properties of the system (aibar, qibar)
  % Solve for aibar
  aibar = zeros(1,n);
  d = 1;
  for k = 1:n
      sum = 0;
      for i = 1:n
          for j = 1:n
              if i==j
                  sum = sum + 0;
              else
                  sum = sum + y(1,j)*aijmatrix(i,j);
              end
          end
      end
      sum = sum + (2*y(1,k)*aijmatrix(k,k)) - amix;
      aibar(1,k)= sum;
  end
  % Solve for qibar
  qibar = zeros(1,n);
  d = 1;
  for b = 1:n
      qibar(1,b) = q*(1 + (aibar(1,b)/amix) - (bi(1,b)/bmix));
  end
  
% Determine the fugacity coefficients and fugacity of each individual species in the mixture (phi, fhat)
  % Solve for lnphi
  lnphi = zeros(1,n);
  d = 1;
  for b = 1:n
      lnphi(1,b)=((bi(1,b)/bmix)*(z - 1) - log(z - beta) - (qibar(1,b)*I));
  end
  % Solve for phi
  phi = zeros(1,n);
  for a = 1:n
      phi(1,a) = exp(lnphi(1,a));
  end
  
  fprintf('The fugacity coefficients of each individual species in the mixture (phi values) are: ');
  fprintf('\n');
  for a = 1:n
      fprintf('phi%d:   %.4f\n', a, phi(1,a));
  end
  fprintf('\n');

  % Solve for fhat
  fhat = zeros(1,n);
  for a = 1:n
    fhat(1,a) = phi(1,a)*press*y(1,a);
  end
  
  fprintf('The fugacity of each individual species in the mixture (fhat values) are: ');
  fprintf('\n');
  for a = 1:n
      fprintf('fhat%d:   %.4f bar \n', a, fhat(1,a)); 
  end 









