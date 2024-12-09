clear; clc;

% establish characteristic dimensions
D=Dimensions(5e3,50,2.99792e8,10);
newcolors = [0.00,0.45,0.74; 1,0,0; 0.39,0.83,0.07; 1,0.41,0.16];
set(groot, "defaultaxescolororder", newcolors);

%% Construct system (USER INPUT HERE) L1 = 5m, L2 = 0m
% simulation information
tmax=D.noDim_t(500e-9);
%Vp2=@(t) (exp(-(t-1.5).^2/0.2));
Vp = @(t) exp(-(t-35./33.4).^2./0.1291);
Rc=D.noDim_R([50, 50]); % cable impedances
u=[0.67, 0.67]; % cable wave speeds
dx=D.noDim_len([5e-2,5e-2]); % target grid spacing
RL=[10]; % load resistances
CL=[0.01]; % load capacitances
Lc=D.noDim_len([5,1e-3]);
name='$L_1=5m,\ L_2=0m$ with Probe @ x = 2.5m';
% for i=1:2
%     switch i
%         case 1
%             Lc=D.noDim_len([5,1e-3]); % cable lengths
%             name='$L_1=5m,\ L_2=0m$';
%         case 2
%             Lc=D.noDim_len([10, 10]); % cable lengths
%             name='$L_1=10m,\ L_2=10m$';
%         % case 3
%         %     Lc=D.noDim_len([10,2]); % cable lengths
%         %     name='$L_1=10m,\ L_2=2m$';
%         % case 4
%         %     Lc=D.noDim_len([1.5,2]); % cable lengths
%         %     name='$L_1=1.5m,\ L_2=2m$';
%     end

    % Store system info
WS=Wire_System(Lc,Rc,u,RL,CL,dx); % store all system information
N=sum(WS.Nx);

% IC setup
%Vp=load_pulse(pulse_file,D);
%Vp=Vp2;
b=@(t)[Vp(t);
    zeros(2*N-1,1)];
J0=Jacobian(0,b(0),WS);
M0=Mass(0,b(0),WS);
[state0,~]=Consistent_IC(J0,M0,b(0),WS);

time=0.0;
cnt=0;
    
J=Jacobian(time,state0,WS);
M=Mass(time,state0,WS);
[state0,dstate0]=Consistent_IC(J,M,state0,WS);
ode=@(t,y) J*y - b(t);
options=odeset('Mass',M,'MassSingular','yes',...
    'RelTol',1e-8,'AbsTol',1e-6,...
    'InitialSlope',dstate0,...
    'Jacobian',J);
    
sol=ode15s(ode,[time tmax],state0, options);
    % Plot Instantaneous
%f2=figure(1);
Nmid1=round(WS.Nx(1)/2);
%Ip1=sol.y(N+Nmid,:);
Vpr=sol.y(Nmid1,:);
Ipr=sol.y(N+Nmid1,:);
VL=sol.y(sum(WS.Nx(1)),:)-sol.y(sum(WS.Nx(1))+1,:);
IL=sol.y(N+sum(WS.Nx(1)),:);
t=sol.x;
    
subplot(2,2,1); hold on
%plotyyy(t,VL,t, IL, t, cumtrapz(t,VL.*IL), {"$V/V_{ch}$", "$I/I_{ch}$", "$E/E_{ch}$"}, 1);
plot(t,Vpr);
plot(t,Ipr);
tspan = linspace(0, tmax, 1001);
plot(tspan, Vp(tspan), '--k');
legend('$V/V{ch}$', 'I/I_{ch}', 'Incident Pulse');
hold off;
hTxt(1,1)=text(-1,-0.27,{'$V/V_{ch}$'},'Interpreter', 'Latex','Rotation',90,'FontSize',16,'color',[0.00,0.45,0.74]);
hTxt(2,1)=text(-1,0.17,{'$I/I_{ch}$'},'Interpreter', 'Latex','Rotation',90,'FontSize',16,'color','r');
% ylabel('$I/I_0\ \mathrm{or}\ V/V_0$', 'Interpreter','Latex');
title(name, 'Interpreter','Latex');
yyaxis right
plot(t,cumtrapz(t,Vpr.*Ipr));
ylim([0 0.6]);
ylabel('$E/E_{ch}$', 'Interpreter','Latex', 'FontSize', 16);
xlabel('Time ($t/t_{ch}$)','Interpreter','Latex', 'FontSize', 16);

name='$L_1=5m,\ L_2=0m$ @ Load (x = 5m)';
subplot(2,2,3); hold on;
plot(t, VL);
plot(t, IL);
tspan = linspace(0, tmax, 1001);
plot(tspan, Vp(tspan), '--k');
legend('$V/V{ch}$', 'I/I_{ch}', 'Incident Pulse');
hold off;
hTxt(1,1)=text(-1,-0.37,{'$V/V_{ch}$'},'Interpreter', 'Latex','Rotation',90,'FontSize',16,'color',[0.00,0.45,0.74]);
hTxt(2,1)=text(-1,0.27,{'$I/I_{ch}$'},'Interpreter', 'Latex','Rotation',90,'FontSize',16,'color','r');
% ylabel('$I/I_0\ \mathrm{or}\ V/V_0$', 'Interpreter','Latex');
title(name, 'Interpreter','Latex');
yyaxis right
plot(t,cumtrapz(t,VL.*IL));
ylim([0 0.6]);
ylabel('$E/E_{ch}$', 'Interpreter','Latex', 'FontSize', 16);
xlabel('Time ($t/t_{ch}$)','Interpreter','Latex', 'FontSize', 16);

%% Construct system (USER INPUT HERE) L1 = 10m, L2 = 10m
% simulation information
tmax=D.noDim_t(500e-9);
%Vp2=@(t) (exp(-(t-1.5).^2/0.2));
Vp = @(t) exp(-(t-35./33.4).^2./0.1291);
Rc=D.noDim_R([50, 50]); % cable impedances
u=[0.67, 0.67]; % cable wave speeds
dx=D.noDim_len([5e-2,5e-2]); % taget grid spacing
RL=[10]; % load resistances
CL=[0.01]; % load capacitances
Lc=D.noDim_len([10,10]);
name='$L_1=10m,\ L_2=10m$ with Probe @ x = 5m';
% for i=1:2
%     switch i
%         case 1
%             Lc=D.noDim_len([5,1e-3]); % cable lengths
%             name='$L_1=5m,\ L_2=0m$';
%         case 2
%             Lc=D.noDim_len([10, 10]); % cable lengths
%             name='$L_1=10m,\ L_2=10m$';
%         % case 3
%         %     Lc=D.noDim_len([10,2]); % cable lengths
%         %     name='$L_1=10m,\ L_2=2m$';
%         % case 4
%         %     Lc=D.noDim_len([1.5,2]); % cable lengths
%         %     name='$L_1=1.5m,\ L_2=2m$';
%     end

    % Store system info
WS=Wire_System(Lc,Rc,u,RL,CL,dx); % store all system information
N=sum(WS.Nx);

% IC setup
%Vp=load_pulse(pulse_file,D);
%Vp=Vp2;
b=@(t)[Vp(t);
    zeros(2*N-1,1)];
J0=Jacobian(0,b(0),WS);
M0=Mass(0,b(0),WS);
[state0,~]=Consistent_IC(J0,M0,b(0),WS);

time=0.0;
cnt=0;
    
J=Jacobian(time,state0,WS);
M=Mass(time,state0,WS);
[state0,dstate0]=Consistent_IC(J,M,state0,WS);
ode=@(t,y) J*y - b(t);
options=odeset('Mass',M,'MassSingular','yes',...
    'RelTol',1e-8,'AbsTol',1e-6,...
    'InitialSlope',dstate0,...
    'Jacobian',J);
    
sol=ode15s(ode,[time tmax],state0, options);
    % Plot Instantaneous
%f2=figure(1);
Nmid1=round(WS.Nx(1)/2);
Nmid2=round(WS.Nx(1)+WS.Nx(2)/2);
Vp1=sol.y(Nmid1,:);
Ip1=sol.y(N+Nmid1,:);
Vp2=sol.y(Nmid2,:);
Ip2=sol.y(N+Nmid2,:);
t=sol.x;
    
subplot(2,2,2); hold on
%plotyyy(t,VL,t, IL, t, cumtrapz(t,VL.*IL), {"$V/V_{ch}$", "$I/I_{ch}$", "$E/E_{ch}$"}, 1);
plot(t,Vp1);
plot(t,Ip1);
tspan = linspace(0, tmax, 1001);
plot(tspan, Vp(tspan), '--k');
legend('$V/V{ch}$', 'I/I_{ch}', 'Incident Pulse');
hold off;
hTxt(1,1)=text(-1,-0.37,{'$V/V_{ch}$'},'Interpreter', 'Latex','Rotation',90,'FontSize',16,'color',[0.00,0.45,0.74]);
hTxt(2,1)=text(-1,0.07,{'$I/I_{ch}$'},'Interpreter', 'Latex','Rotation',90,'FontSize',16,'color','r');
% ylabel('$I/I_0\ \mathrm{or}\ V/V_0$', 'Interpreter','Latex');
title(name, 'Interpreter','Latex');
yyaxis right;
plot(t,cumtrapz(t,Vp1.*Ip1));
ylim([0 0.6]);
ylabel('$E/E_{ch}$', 'Interpreter','Latex', 'FontSize', 16);
xlabel('Time ($t/t_{ch}$)','Interpreter','Latex', 'FontSize', 16);


name='$L_1=10m,\ L_2=10m$ with Probe @ x = 15m';
subplot(2,2,4); hold on;
plot(t, Vp2);
plot(t, Ip2);
tspan = linspace(0, tmax, 1001);
plot(tspan, Vp(tspan), '--k');
legend('', '', 'Incident Pulse', '');
hold off;
hTxt(1,1)=text(-1,-0.37,{'$V/V_{ch}$'},'Interpreter', 'Latex','Rotation',90,'FontSize',16,'color',[0.00,0.45,0.74]);
hTxt(2,1)=text(-1,0.07,{'$I/I_{ch}$'},'Interpreter', 'Latex','Rotation',90,'FontSize',16,'color','r');
% ylabel('$I/I_0\ \mathrm{or}\ V/V_0$', 'Interpreter','Latex');
title(name, 'Interpreter','Latex');
yyaxis right;
plot(t,cumtrapz(t,Vp2.*Ip2));
ylim([0 0.6]);
ylabel('$E/E_{ch}$', 'Interpreter','Latex', 'FontSize', 16);
xlabel('Time ($t/t_{ch}$)','Interpreter','Latex', 'FontSize', 16);

%% Save figure
%f2.Units='centimeters';
%f2.Position=[2,2,20,20];
%f2.PaperUnits='centimeters';
%f2.PaperSize=[20,20];
%print(f2,'figures/Case_Study4b.pdf','-dpdf');