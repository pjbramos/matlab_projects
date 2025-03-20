%converting volume percentage of ethanol to mole fraction:
    liters_gin=10;
    percent_ethanol_gin=0.4;
    percent_water_gin=1-percent_ethanol_gin;
    mol_ethanol=liters_gin*percent_ethanol_gin*789/46.07;
    mol_water=liters_gin*percent_water_gin*1000/18.016;
    mol_total=mol_ethanol+mol_water;
    ethanol_mol_fraction=mol_ethanol/mol_total;
    water_mol_fraction=mol_water/mol_total;

    
%calculating mole fraction in liquid and vapor phase, looping for the 5 runs
z1=ethanol_mol_fraction; z2=water_mol_fraction;
for i=1:5
    A1=16.8958; B1=3795.17; C1=230.918;
    A2=16.3872; B2=3885.70; C2=230.170;
    pressure=101.325;
    bubble_temp=fzero(@(bubble_temp) pressure-z1*exp(A1-B1/(bubble_temp+C1))-z2*exp(A2-B2/(bubble_temp+C2)), 75);
    dew_temp=fzero(@(dew_temp) 1/pressure-z1/exp(A1-B1/(dew_temp+C1))-z2/exp(A2-B2/(dew_temp+C2)), 75);
    temp=input(['Input temperature of run ' num2str(i) ' (in Celsius, between ' num2str(bubble_temp) ' and ' num2str(dew_temp) '): ']);
    sat_pressure_ethanol=exp(A1-(B1/(temp+C1)));
    sat_pressure_water=exp(A2-(B2/(temp+C2)));
    K_ethanol=sat_pressure_ethanol/pressure;
    K_water=sat_pressure_water/pressure;
    k_ethanol=1/(K_ethanol-1); k_water=1/(K_water-1);
    a=fzero(@(a) z1/(k_ethanol+a) + z2/(k_water+a), 1);
    x1=z1/(1+a*(K_ethanol-1)); x2=z2/(1+a*(K_water-1));
    y1=K_ethanol*x1; y2=K_water*x2;

    %performing mass balance on the system
    eqn1=[1, 1; x1, y1]; eqn2=[mol_total; z1*mol_total];
    sol=linsolve(eqn1,eqn2);
    mol_total=sol(2);
    volume_ethanol=mol_total*y1*46.07/789;
    volume_water=mol_total*y2*18.016/1000;
    volume_total=volume_ethanol+volume_water;
    fprintf('Total volume (liters): %.3f\nVolume of ethanol (liters): %.3f, %%v/v ethanol: %.2f%%\nVolume of water (liters): %.3f, %%v/v water: %.2f%%\n', volume_total, volume_ethanol, volume_ethanol/volume_total*100, volume_water, volume_water/volume_total*100);
    
    z1=y1; z2=y2;
end
