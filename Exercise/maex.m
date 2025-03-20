%mex 1

    %     function maex1
%     % Set simulation time:
%     tspan = 0:0.01:100; clc;     % t = 0 min to 30 min (Time step = 0.01)
%     
%     % Input variables, u
%     u = zeros(4,length(tspan)); % Size: [# of inputs] x [# of time steps]
%     u(1,:) = 200;               % F_in = 200 L/min
%     u(2,:) = 1;               % Ca,in = 1 mol/L
%     u(3,:) = 0;              % Cb,in = 0 mol/L
%     
%     % Output variables, y
%     y0(1) = 4;                  % h = 4 m  (initial steady-state)
%     y0(2) = 0.2;                % Ca,out = 0.2 (initial steady-state)
%     y0(3) = 0.8;                % Cb,out = 0.8 (initial steady-state)
%     
%     u(1,tspan >= 5) = 250;       % Step inc. F from 200 to 250 at 5 min
%     
%     % Solve the ODE
%     options = odeset('RelTol',1e-8);            % Set the solver accuracy
%     [t,y] = ode45(@odefcn,tspan,y0,options);    % Solve for y(t)
%     
%     % Plot the results
%     subplot(411); plot(t,u(1,:),'b','LineWidth',1.2); grid on;
%     ylabel('F,in (L/min)'); axis([0 100 150 300]);
%     
%     subplot(412); plot(t,y(:,1),'b','LineWidth',1.2); grid on;
%     ylabel('h (m)');
%     
%     subplot(413); plot(t,y(:,2),'b','LineWidth',1.2); grid on;
%     ylabel('Ca,out (mol/L)');
% 
%     subplot(414); plot(t,y(:,3),'b','LineWidth',1.2); grid on;
%     xlabel('Time (min)'); ylabel('Cb,out (mol/L)');
% 
%     function dydt = odefcn(t,y)
%         % System of ODEs to be solved
%         
%         % Initialize dydt as a column of states, x
%         dydt = zeros(3,1);
%         
%         % Determine the desired input values at time, t
%         ut = interp1(tspan,u',t);
%         
%         % Constant parameters:
%         A = 0.5;      % m^2, area
%         rho = 1000;  % kg/m^3, density
%         k = 2.0;   % L*mol^-1*min^-1
%                 
%         dydt(1) = 1./(A*rho).*(ut(1) - 100.*sqrt(y(1)));  % Eq. (1)
%         dydt(2) = ut(1)./(A*rho*y(1)).*(ut(2) - y(2)) - k.*y(2).*y(2);  % Eq. (2)
%         dydt(3) = ut(1)./(A*rho*y(1)).*(ut(3) - y(3)) + k.*y(2).*y(2);  % Eq. (3)
%     end
% 
% end

    %  function maex12
%     % Input conditions
%     h=0;            %height, m, initial steady-state
%     u=0;            %flow rate, m^3/s, initial steady-state
%     tin=0;          %time at runtime start, s
%     tout=10;        %time at runtime end,s
%     M=1;            %step increase
%     
%     %time span
%     tbef=0:0.01:tin;    %from runtime start to disturbance start
%     taft=tin:0.01:tout;   %from disturbance start to runtime end
% 
%     %input eqns
%     ubef=u.*(ones(size(tbef)));     %flow rate before the disturbance
%     uaft=(u+M).*(ones(size(taft))); %flow rate after disturbance
% 
%     %plotting u(t)
%     subplot(311); tplot=[tbef taft]; uplot=[ubef uaft];
%     plot(tplot,uplot,'b','LineWidth',1); grid on;
%     ylabel('Flow rate (m^3/s)'); axis([0 10 0 2]);
% 
%     %solving for 2a (output eqns)
%     hbef=h.*(ones(size(tbef)));                 %height before the disturbance
%     haft=M.*(1-((1-(3.*taft)).*exp(-taft)));    %height after disturbance
% 
%     %plotting 2a
%     subplot(312); hplot=[hbef haft];
%     plot(tplot,hplot,'b','LineWidth',1); grid on;
%     xlabel('Time (s)'); ylabel('Height (m)'); axis([0 10 0 2]);
% 
%     %solving for 2b
%     tover=0;
%     hover=0;
%     htank=1.5;      %desired tank height
%     while htank-hover>1e-5
%         tover=tover+0.01;
%         hover=M.*(1-((1-(3.*tover)).*exp(-tover)));
%     end
%     fprintf('overflow height: %.3f m\noverflow time: %.2f s\n', hover, tover);
% 
%     %solving for 2c
%     Mover=0;
%     hmax=0;
%     while htank-hmax>1e-5
%         Mover=Mover+1e-5;
%         a=Mover.*(1-((1-(3.*taft)).*exp(-taft))); abef=h.*(ones(size(tbef)));
%         aplot=[abef a];
%         hmax=max(aplot);
%     end
%     fprintf('max step increase: %.3f m\n %.2f s\n', Mover);
% 
%     %plotting 2c
%     subplot(313); plot(tplot,aplot,'b','LineWidth',1); grid on;
%     xlabel('Time (s)'); ylabel('Height (m)'); axis([0 10 0 2]);
% end

%mex 2

    % a=zpk(-86.2069, [-116.081, -60.2144], 53.5138);
% step(a); grid on;
% 
% [sR,sI] = meshgrid(-150:1:-30,-40:1:40);
% G = abs(polyval([0.0077372 0.667],sR+sI*i)./polyval([0.000143066 0.02522 1],sR+sI*i));
% contourf(sR,sI,G,500); caxis([0 5]); colorbar; grid on;
% 
% 
% b=tf([0.1667], [0.00014317 0.02522 1]);
% step(b); grid on;
% 
% [sR,sI] = meshgrid(-150:1:-30,-40:1:40);
% G = abs(polyval([0.1667],sR+sI*i)./polyval([0.00014317 0.02522 1],sR+sI*i));
% contourf(sR,sI,G,500); caxis([0 5]); colorbar; grid on;


    % s=tf('s');
% x=-0.15/(6*s+1)^5; step(x); hold on; grid on;
% y=(-0.15*exp(-21*s))/(9*s+1); step(y); hold on;
% z=(-0.15*exp(-15*s))/((6*s+1)*(9*s+1)); step(z);
% legend('Original TF','FOPTD Approx.', 'SOPTD Approx.');
% 
% 
% s=tf('s');
% x=(-0.15*exp(-35.584*s))/(6*s+1)^5; step(x); hold on; grid on;
% y=(-0.15*exp(-(35.584+21)*s))/(9*s+1); step(y); hold on;
% z=(-0.15*exp(-(35.584+15)*s))/((6*s+1)*(9*s+1)); step(z);
% legend('Original TF','FOPTD Approx.', 'SOPTD Approx.');


%mex 4
tc=5; yt=1-exp(-t/tc); 
t=0:0.1:20; plot (t, yt); hold on;
plot(out.simulink.Time, out.simulink.Data); hold on;
plot(out.simulink1.Time, out.simulink1.Data);


Kc = linspace(-500,500,50000);
G = tf([51.5 15.1 1], [75 772.5 226.5 15 0]);
rlocus(G,Kc);
