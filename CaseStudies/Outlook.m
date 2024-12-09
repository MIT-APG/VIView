clear; clc;

% establish characteristic dimensions
D=Dimensions(1e4,50,3e8,10);

% choose the "truth" values
RL_CL_Vb=[rand()*25, rand()*1e-10, 1+rand()];


%% Construct system (USER INPUT HERE)
% Wire information
Lc=D.noDim_len([5,5,1e-3]); % cable lengths
Rc=D.noDim_R([50, 50, 50]); % cable impedances
u=[0.6,0.6,0.6]; % cable wave speeds
% Load information
RL=D.noDim_R([0,0]); % load resistances
CL=D.noDim_C([0,RL_CL_Vb(2)]); % load capacitances
dx=D.noDim_len([5e-2,5e-2,5e-2]); % taget grid spacing
% simulation information
tmax=5;
Vp2=@(t) (exp(-(t-1.5).^2/0.2));

%% Store system info
WS=Wire_System(Lc,Rc,u,RL,CL,dx); % store all system information
N=sum(WS.Nx);

%% IC setup
b=@(t)[Vp2(t);
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
        SE=@(t,y) spark_event(t,y,WS,2,RL_CL_Vb(3),1);
    else
        WS.RL=D.noDim_R([0,RL_CL_Vb(1)]);
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

%% Extract the probe measurements
t=linspace(0,tmax, 10001);
Is=0*t; IL=0*t;
VL=0*t; Vs=0*t;
for i=1:length(sol)
    tloc=sol{i}.x;
    rng=logical((t>=tloc(1)).*(t<=tloc(end)));
    Is(rng) = deval(sol{i},t(rng), WS.Nx(1)+N);
    Vs(rng) = deval(sol{i},t(rng), WS.Nx(1));
    VL(rng) = deval(sol{i},t(rng), sum(WS.Nx(1:2)))-...
        deval(sol{i},t(rng), sum(WS.Nx(1:2))+1);
    IL(rng) = deval(sol{i},t(rng), sum(WS.Nx(1:2))+N);
end

%% FFT
xsamp=Lc(1);
[freq,g]=FFT_complex(t,Vs);
[~,h]=FFT_complex(t,Is);
T=t(end);
phi=2*pi*freq/(u(1)) *xsamp;

xtarget=Lc(1)+Lc(2);


len=length(phi);
a=(g(1:len) .* exp(1i*phi)+h(1:len) .* exp(1i*phi))/2;
b=(g(1:len) - a .* exp(-1i*phi)) ./ exp(1i*phi);

ampA=abs(a);
ampA(1)=a(1);
ampA(2:end)=2*ampA(2:end);
ampB=abs(b);
ampB(2:end)=2*ampB(2:end);
ampB(1)=b(1);

phaseA=atan2(imag(a),real(a))-phi*xtarget/xsamp;
phaseB=atan2(imag(b),real(b))+phi*xtarget/xsamp;
Vf=FFT_reconstruct(t, freq, ampA, phaseA, 200);
Vr=FFT_reconstruct(t, freq, ampB, phaseB, 200);

figure(1); clf;
plot(t,VL)
hold on
plot(t,Vf)
plot(t,Vr)
plot(t,Vf+Vr)




%%
aL=a.*exp(-1i.*phi.*xtarget ./ xsamp);

aL(100:end)=0;


daL=aL .* 1i*2*pi .*freq;
IaL=aL .* -1i ./ (2*pi*freq);
aL_full=[aL, conj(aL(end:-1:2))];
daL_full=[daL, conj(daL(end:-1:2))];
IaL_full=[IaL, conj(IaL(end:-1:2))];
VfL=ifft(aL_full,'symmetric')*length(t);
dVfL=ifft(daL_full,'symmetric')*length(t);
IVfL=ifft(IaL_full,'symmetric')*length(t);

bL=b.*exp(1i*phi.*xtarget ./ xsamp);

bL(200:end)=0;

dbL=bL .* 1i*2*pi .*freq;
IbL=bL .* -1i ./ (2*pi*freq);
bL_full=[bL, conj(bL(end:-1:2))];
dbL_full=[dbL, conj(dbL(end:-1:2))];
IbL_full=[IbL, conj(IbL(end:-1:2))];
VrL=ifft(bL_full,'symmetric')*length(t);
dVrL=ifft(dbL_full,'symmetric')*length(t);
IVrL=ifft(IbL_full,'symmetric')*length(t);


rL=bL./aL;
figure(2); clf;
subplot(3,1,1)
plot(t,Vf); hold on;
plot(t,Vr);
plot(t,VfL);
plot(t,VrL);
subplot(3,1,2);
plot(t,dVfL);
hold on
plot(t,dVrL);
ylim([-2,2]);
subplot(3,1,3);
plot(t,IVfL);
hold on
plot(t,IVrL);



Vtot=aL_full+bL_full;
Itot=(aL_full-bL_full);
figure(3); clf;
subplot(2,1,1);
plot(t,VL);
hold on;
plot(t,ifft(Vtot)*length(t));
subplot(2,1,2);
plot(t,IL);
hold on;
plot(t,ifft(Itot)*length(t));

Iguess1=VL/1e6 + CL(2)*(dVfL+dVrL);
Iguess2=VL/D.noDim_R(RL_CL_Vb(1)) + CL(2)*(dVfL+dVrL);

tflip=3.16;
Iguess=Iguess1*0.5 .*(tanh(-(t-tflip)*500)+1) + ...
    Iguess2*0.5 .*(tanh((t-tflip)*500)+1);

%fiteq=fittype(@(r,c,tf,x) ((VL ./ 1e6) .* (tanh(-(x-tf)*100)+1) + ...
%    (VL ./ r) .* (tanh(-(x-tf)*100)+1) + ...
%    c*(dVfL+dVrL))');

dVL=dVfL+dVrL;
VL_func=@(ts) interp1(t,VL, ts);
dVL_func=@(ts) interp1(t,dVL, ts);

%fiteq=@(r,x) VL_func(x)/r;
fiteq=fittype(@(r,c,tf,x) ((VL_func(x) ./ 1e6) .* 0.5 .* (tanh(-(x-tf)*500)+1) + ...
   (VL_func(x) ./ r) .* 0.5 .* (tanh((x-tf)*500)+1) + ...
   c*dVL_func(x)));

ft=fit(t',IL',fiteq, 'StartPoint',[0.2,0.05,3.1]);
    
Iguess=fiteq(D.noDim_R(RL_CL_Vb(1)),D.noDim_C(RL_CL_Vb(2)),3.16,t');
plot(t,Iguess);
plot(t,ft(t'));
ylim([-1,2]);




%% Actual rLomega
nu=D.noDim_R(RL_CL_Vb(1));
eta=-1 ./ (D.noDim_C(RL_CL_Vb(2))*2*pi*freq);
rLomega=1./((nu+1)^2+eta.^2) .* ((nu.^2-1+eta .^2 ) - 1i*(2*eta));
nu2=1e6;
rLomega2=1./((nu2+1)^2+eta.^2) .* ((nu2.^2-1+eta .^2 ) - 1i*(2*eta));

figure(2); clf;
plot(rLomega(1:200),'or'); hold on;
plot(rLomega2(1:200),'ob');

%% Plot Instantaneous
f2=figure(2); clf;
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
        xlim([0,6]);
    else
        xlim([0,5]);
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
% plot(sol.x,V,tgrid,V_filt);