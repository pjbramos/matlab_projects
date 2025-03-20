clear;
clc;

%asking for input
components=input('Input number of components: ');
y=zeros(1, components);
temperature=input('\nInput mixture temperature (in Kelvin): ');
pressure=input('Input mixture pressure (in bar): ');
crit_temp=zeros(1, components);
crit_pressure=zeros(1, components);
z_crit=zeros(1, components);
v_crit=zeros(1, components);
w=zeros(1, components);

for i=1:components
    y(1,i)=input(['\nInput gas/vapor composition y' num2str(i) ': ']);
    w(1,i)=input(['Input w' num2str(i) ': ']);
    crit_temp(1,i)=input(['Input critical temperature (in Kelvin) of component ' num2str(i) ': ']);
    crit_pressure(1,i)=input(['Input critical pressure (in bar) of component ' num2str(i) ': ']);
    z_crit(1,i)=input(['Input Zc of component ' num2str(i) ': ']);
    v_crit(1,i)=input(['Input Vc (in cm3/mol) of component ' num2str(i) ': ']);
end

%solving for bij
%Tcij, Trij, Zcij, wcij, Vcij, Pcij are required
crit_temp_final=zeros(components, components);
crit_pressure_final=zeros(components, components);
z_crit_final=zeros(components, components);
v_crit_final=zeros(components, components);
w_final=zeros(components, components);
b=zeros(components, components);
reduced_temp=zeros(components, components);
for i=1:components
    for j=1:components
        crit_temp_final(i,j)=sqrt(crit_temp(1,i)*crit_temp(1,j));
        w_final(i,j)=(w(1,i)+w(1,j))/2;
        z_crit_final(i,j)=(z_crit(1,i)+z_crit(1,j))/2;
        v_crit_final(i,j)=((v_crit(1,i)^(1/3)+v_crit(1,j)^(1/3))/2)^3;
        crit_pressure_final(i,j)=(83.14*z_crit_final(i,j)*crit_temp_final(i,j))/v_crit_final(i,j);
        reduced_temp(i,j)=temperature/crit_temp_final(i,j);
        b(i,j)=(0.083-(0.422/(reduced_temp(i,j)^1.6))+w_final(i,j)*(0.139-(0.172/(reduced_temp(i,j)^4.2))))*83.14*crit_temp_final(i,j)/crit_pressure_final(i,j);
    end
end

%output bij
fprintf('\nThe possible Bij values are: \n');
for i=1:components-1
    for j=i+1:components
        fprintf('B%d%d: %g\n', i, j, b(i,j));
    end
end

%solving for delta
delta_matrix=zeros(components, components);
for i=1:components
    for j=1:components
        if i==j
           delta_matrix=delta_matrix;
        else
            delta_matrix(i,j)=2*b(i,j)-b(i,i)-b(j,j);
        end
    end
end

%solving for fugacity, output phi hat and f hat
for k=1:components
    sum=0;
    for i=1:components
        for j=1:components
            sum=sum+(y(1,i)*y(1,j)*(2*delta_matrix(i,k)-delta_matrix(i,j)));
        end
    end
    ln_phi=pressure/83.14/temperature*(b(k,k)+0.5*sum);
    phi=exp(ln_phi);
    fprintf('\nphi_hat_%d = %g\n', k, phi);
    fugacity=phi*y(1,k)*pressure;
    fprintf('f_hat_%d = %g bar\n', k, fugacity);
end