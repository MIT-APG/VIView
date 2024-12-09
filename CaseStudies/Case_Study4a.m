clear; clc;

% establish characteristic dimensions
D=Dimensions(5e3,50,2.99792e8,10);
newcolors = [0.00,0.45,0.74; 1,0,0; 0.39,0.83,0.07; 1,0.41,0.16];
set(groot, "defaultaxescolororder", newcolors);
%% Construct system (USER INPUT HERE)
% Wire information
Lc=D.noDim_len([10,6,1e-3]); % cable lengths
Rc=D.noDim_R([50,100,50]); % cable impedances
u=[0.67, 0.67, 0.67]; % cable wave speeds
dx=D.noDim_len([5e-2,5e-2, 5e-3]); % taget grid spacing
% simulation information
tmax=D.noDim_t(500e-9);
%tmax=5;
%pulse_file='../../Data/Ref_Waveforms/Pulse_300VDC.mat';
Vp = @(t) exp(-(t-35./33.4).^2./0.1291);
%Vp =@(t) (exp(-(t-1.5).^2/0.2));
for i=1:4
    switch i
        case 1
            RL=[0,10]; % load resistances
            CL=[0,0.1]; % load capacitances
            name='Load A';
        case 2
            RL=[0,10]; % load resistances
            CL=[0,0.01]; % load capacitances
            name='Load B';
        case 3
            RL=[0,0.5]; % load resistances
            CL=[0,0.1]; % load capacitances
            name='Load C';
        case 4
            RL=[0,0.5]; % load resistances
            CL=[0,0.01]; % load capacitances
            name='Load D';
    end
%RL=[0,10]; % load resistances
%CL=[0,0.01]; % load capacitances

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
        'RelTol',1e-10,'AbsTol',1e-8,...
        'InitialSlope',dstate0,...
        'Jacobian',J);

    sol=ode15s(ode,[time tmax],state0, options);
% Plot Instantaneous
f2=figure(2);
if i==1
    clf;
    leg={};
end
% Nmid=round(WS.Nx(1)/2);
% Nmid2=round(WS.Nx(2)/2)+WS.Nx(1);
% Vp1=sol.y(Nmid,:);
% Vp2=sol.y(Nmid2,:);
% Ip1=sol.y(N+Nmid,:);
% Ip2=sol.y(N+Nmid2,:);
VL=sol.y(sum(WS.Nx(1:2)),:)-sol.y(sum(WS.Nx(1:2))+1,:);
IL=sol.y(N+sum(WS.Nx(1:2)),:);
t=sol.x;

subplot(3,1,1); hold on
%plot(t,Vp1,t,Vp2);
plot(t,VL);
ylabel('$V/V_{ch}$', 'Interpreter','Latex');
subplot(3,1,2); hold on
plot(t,IL);
ylabel('$I/I_{ch}$', 'Interpreter','Latex');
subplot(3,1,3); hold on
plot(t,cumtrapz(t, IL.*VL));
ylabel('$E/E_{ch}$', 'Interpreter','Latex');

leg{i}=name;
if i==4
for k=1:3
    subplot(3,1,k)
    xlabel('Time ($t/t_{ch}$)','Interpreter','Latex');
    legend(leg{:});
end
end
end
tspan = linspace(0, tmax, 1001);
subplot(3, 1, 1);
hold on;
plot(tspan, Vp(tspan), '--k');
legend('Load A', 'Load B', 'Load C', 'Load D', 'Incident Pulse');
%for i=1:2
   %subplot(2,1,i);
   %legend('Probe 1','Probe 2','Load');
   %xlabel('$t/t_0$','Interpreter','Latex');
%end
%subplot(2,1,1);
%text([4.5, 7, 8.5, 9.5,11], [5,-2.5,-2.5,1,-4],{'A','B','C','D','E'});
%ylim([-6,6]);

%% Save figure
f2.Units='centimeters';
f2.Position=[2,2,20,20];
f2.PaperUnits='centimeters';
f2.PaperSize=[20,20];
print(f2,'Case_Study4a.pdf','-dpdf');