clear; clc;

% establish characteristic dimensions
D=Dimensions(5e3,50,2.99792e8,10);

%% Construct system (USER INPUT HERE)
% Wire information
Lc=D.noDim_len([5,5,1e-3]); % cable lengths
Rc=D.noDim_R([50, 50, 50]); % cable impedances
u=[0.67,0.67,0.67]; % cable wave speeds
% Load information
RL=D.noDim_R([0,0]); % load resistances
CL=D.noDim_C([0, 6.67e-12]); % load capacitances
dx=D.noDim_len([5e-2,5e-2,5e-3]); % taget grid spacing
% simulation information
tmax=10; %non-dimensional time
%Vp2=@(t) (exp(-(t-1.5).^2/0.2));
Vp = @(t) exp(-(t-35./33.4).^2./0.1291);

%% Store system info
WS=Wire_System(Lc,Rc,u,RL,CL,dx); % store all system information
N=sum(WS.Nx);

%% IC setup
for k=1:2
b=@(t)[Vp(t);
    zeros(2*N-1,1)];
J0=Jacobian(0,b(0),WS);
M0=Mass(0,b(0),WS);
[state0,dstate0]=Consistent_IC(J0,M0,b(0),WS);

time=0.0;
cnt=0;
sol={};

while time<tmax
    if mod(cnt,2)==0
        WS.RL=[0,1e6];
        if k==1
            SE=@(t,y) spark_event(t,y,WS,2,4,1);
        else
        SE=@(t,y) spark_event(t,y,WS,2,1.5,1);
        end
   else
        WS.RL=[0,0.2];
        SE=@(t,y) spark_event(t,y,WS,2,0.0,-1);
    end
    cnt=cnt+1;
    J=Jacobian(time,state0,WS);
    M=Mass(time,state0,WS);
    [state0,dstate0]=Consistent_IC(J,M,state0,WS);
    ode=@(t,y) J*y - b(t);
    options=odeset('Mass',M,'MassSingular','yes',...
        'RelTol',1e-10,'AbsTol',1e-8,...
        'InitialSlope',dstate0,...
        'Jacobian',J,...
        'Events',SE);
    sol{cnt}=ode15s(ode,[time tmax],state0, options);
    time=sol{cnt}.x(end);  
    state0=sol{cnt}.ye;
    if cnt==1
        fprintf('Spark gap closed at t=%.2f\n',time);
    else
        fprintf('Spark gap opened at t=%.2f\n',time);
    end
end

%% Plot Video
f1=figure(1);clf;
Np=101;
tgrid=linspace(0,tmax,Np);
cnt=1;
for i=1:Np
    if sol{cnt}.x(end)<tgrid(i)
        cnt=cnt+1;
    end
    y=deval(sol{cnt},tgrid(i));
    plot(WS.x',y(1:N),WS.x',y(N+1:end));
    %legend('V','I');
    tit=sprintf('t=%u', D.reDim_t(tgrid(i)));
    title(tit);
    ylim([-1,1]);
    pause(0.02)
    drawnow()
end

%% Plot Instantaneous
if k==1
f2=figure(2); clf;
end
Is=[]; IL=[];
VL=[]; Vs=[];
t=[];
for i=1:length(sol)
    Is=[Is sol{i}.y(WS.Nx(1)+N,:)];
    Vs=[Vs sol{i}.y(WS.Nx(1),:)];
    VL=[VL sol{i}.y(sum(WS.Nx(1:2)),:) - sol{i}.y(sum(WS.Nx(1:2))+1,:)];
    IL=[IL sol{i}.y(sum(WS.Nx(1:2))+N,:)];
    t=[t sol{i}.x];
end
end

subplot(3,2,3); hold on
plot(t,Is);
title('I_M')
subplot(3,2,4); hold on
plot(t,IL)
title('I_L');
subplot(3,2,1); hold on
plot(t,Vs);
title('V_M');
subplot(3,2,2); hold on
plot(t,VL)
title('V_L');
subplot(3,2,5); hold on
plot(t,cumtrapz(t,Vs .* Is));
title('E_M');
subplot(3,2,6); hold on
plot(t,cumtrapz(t,VL .* IL));
title('E_L');

%%
for i=1:6
    subplot(3,2,i);
    xlabel('$t/t_0$', 'Interpreter','Latex');
    if i<3
        ylabel('$V/V_0$','Interpreter','Latex');
        ylim([-2,2]);
    elseif i<5
        ylabel('$I/I_0$','Interpreter','Latex');
        ylim([-2,2]);
    else
        ylabel('$E/E_0$','Interpreter','Latex');
        ylim([-0.2,0.6]);
    end
    if i==2
        legend('No discharge', 'Discharge','Location','NorthWest');
    end
    if mod(i,2)~=0
        xlim([0,tmax]);
    else
        xlim([0,tmax]);
    end
end

%% Save figure
f2.Units='centimeters';
f2.Position=[2,2,20,16];
f2.PaperUnits='centimeters';
f2.PaperSize=[20,16];
print(f2,'figures/Case_Study6.pdf','-dpdf');

%%
% plot(t,Is, t,VL);
% for i=1:length(sol)-1
%     plot([sol{i}.xe,sol{i}.xe], [-1,1],'--k');
% end
% 
% subplot(2,1,2); hold on
% plot(D.reDim_t(t)*1e9,D.reDim_E(cumtrapz(t,Is.*Vs))*1e3);
% xlabel('Time (ns)');
% ylabel('Energy (mJ)');


%% Fourier Decomposition
% tgrid=linspace(0,tmax,1001);
% Vgrid=deval(sol,tgrid,sum(WS.Nx(1:2)))-deval(sol,tgrid,sum(WS.Nx(1:3)));
% [freq,amp,phase]=FFT_processed(D.reDim_t(tgrid),Vgrid);
% f3=figure(3); clf;
% subplot(2,1,1);
% plot(freq/1e6,amp);
% xlim([0,400]);
% 
% V_filt=FFT_reconstruct(D.reDim_t(tgrid),freq,amp,phase,100);
% subplot(2,1,2);
% plot(sol.x,V,tgrid,V_filt)