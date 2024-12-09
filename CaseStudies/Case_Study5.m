clear; clc;
newcolors = [0.00,0.45,0.74; 1,0,0; 0.39,0.83,0.07; 1,0.41,0.16];
set(groot, "defaultaxescolororder", newcolors);

% establish characteristic dimensions just for normalization, as in paper
D=Dimensions(5e3,50,2.99792e8,10); %5kV, 50 Ohms, speed of light, and 10 m

%% Construct system (USER INPUT HERE)
% Wire information
Lc=D.noDim_len([10,1e-3]); % cable lengths (m)
Rc=D.noDim_R([50, 50]); % cable impedances (Ohms)
u=[0.67, 0.67]; % cable wave speeds u = v/c, c is speed of light
dx=D.noDim_len([5e-2,5e-3]); % target grid spacing, m/pt
% simulation information
%tmax=D.noDim_t(500e-9);%
tmax=15; %non dimensional time
%pulse_file='../../Data/Ref_Waveforms/Pulse_300VDC.mat'; <-- experimental
%file
Vp = @(t) exp(-(t-35./33.4).^2./0.1291); %sample gaussian used in paper

RL=[10]; % load resistances
CL=[0.1]; % load capacitances

%System looks like: power supply --> cable (10m) --> load --> grounding
%cable (1e-3m)
%We put a small grounding cable at the end to help apply end boundary
%condition as ground V(x = L) = 0

%% Store system info
WS=Wire_System(Lc,Rc,u,RL,CL,dx); % store all system information
N=sum(WS.Nx);

%% IC setup
%Vp=load_pulse(pulse_file,D);
%Vp=Vp2;
b=@(t)[Vp(t);
    zeros(2*N-1,1)];
J0=Jacobian(0,b(0),WS);
M0=Mass(0,b(0),WS);
[state0,dstate0]=Consistent_IC(J0,M0,b(0),WS);

time=0.0;
cnt=0;
sol={};

J=Jacobian(time,state0,WS);
M=Mass(time,state0,WS);
[state0,dstate0]=Consistent_IC(J,M,state0,WS);
ode=@(t,y) J*y - b(t);
options=odeset('Mass',M,'MassSingular','yes',...
    'RelTol',1e-8,'AbsTol',1e-6,...
    'InitialSlope',dstate0,...
    'Jacobian',J);

sol=ode15s(ode,[time tmax],state0, options);
%% Waveform Decomposition
% Vpr=sol.y(101,:);
% Ip=sol.y(N+151,:);
% VL=sol.y(sum(WS.Nx(1)),:)-sol.y(sum(WS.Nx(1))+1,:);
% IL=sol.y(N+sum(WS.Nx(1)),:);
% t=sol.x;

t=linspace(0,sol.x(end),10001);
VA=deval(sol,t,101);
IA=deval(sol,t,N+101);
VB=deval(sol,t, 101);
IB=deval(sol,t, N+101);

Iloc2=deval(sol,t,N+151);
Vloc3=deval(sol,t,201);

xsamp=[0.5,0.5];
[freq,g]=FFT_complex(t,VA);
[~,h]=FFT_complex(t,IB);
T=t(end);
phiA=2*pi*freq/(u(1)) *xsamp(1);
phiB=2*pi*freq/(u(1)) *xsamp(2);

len=length(phiA);
a=(g(1:len) .* exp(1i*phiB)+h(1:len) .* exp(1i*phiA)) ./ ...
    (2 * cos(phiB-phiA));
b=(g(1:len) - a .* exp(-1i*phiA)) ./ exp(1i*phiA);

ampA=abs(a);
ampA(1)=a(1);
ampA(2:end)=2*ampA(2:end);
ampB=abs(b);
ampB(2:end)=2*ampB(2:end);
ampB(1)=b(1);

xtarget=xsamp(1);
phaseA=atan2(imag(a),real(a))-phiA*xtarget/xsamp(1);
phaseB=atan2(imag(b),real(b))+phiB*xtarget/xsamp(2);
Vf_loc1=FFT_reconstruct(t, freq, ampA, phaseA, 200);
Vr_loc1=FFT_reconstruct(t, freq, ampB, phaseB, 200);

xtarget=0.75;
phaseA=atan2(imag(a),real(a))-phiA*xtarget/xsamp(1);
phaseB=atan2(imag(b),real(b))+phiB*xtarget/xsamp(2);
Vf_loc2=FFT_reconstruct(t, freq, ampA, phaseA, 200);
Vr_loc2=FFT_reconstruct(t, freq, ampB, phaseB, 200);

xtarget=1.0;
phaseA=atan2(imag(a),real(a))-phiA*xtarget/xsamp(1);
phaseB=atan2(imag(b),real(b))+phiB*xtarget/xsamp(2);
Vf_loc3=FFT_reconstruct(t, freq, ampA, phaseA, 200);
Vr_loc3=FFT_reconstruct(t, freq, ampB, phaseB, 200);
%%
f2=figure(1); clf;
subplot(2,2,1); hold on
box on;
plot(t,VA);
plot(t,IB);
ylabel('$V/V_{ch},\ I/I_{ch}$','Interpreter','Latex');
legend('Voltage','Current', 'Interpreter','Latex');
subplot(2,2,2);
plot(t,Vf_loc1); hold on;
plot(t,Vr_loc1);
ylabel('$V/V_{ch}$','Interpreter','Latex');
legend('$V_f(x_A)$','$V_r(x_A)$', 'Interpreter','Latex');
subplot(2,2,3);
plot(t,Iloc2); hold on
plot(t,Vf_loc2-Vr_loc2,'--'); hold on
legend('Actual','Reconstructed', 'Interpreter','Latex');
title('$I/I_{ch}$ Reconstructed @ x = 7.5m', 'Interpreter','Latex');
ylabel('$I/I_{ch}$','Interpreter','Latex');
subplot(2,2,4);
plot(t,Vloc3); hold on
plot(t,Vf_loc3+Vr_loc3,'--'); hold on
legend('Actual','Reconstructed', 'Interpreter','Latex');
title('$V/V_{ch}$ Reconstructed @ Load (x = 10m)', 'Interpreter','Latex');
ylabel('$V/V_{ch}$','Interpreter','Latex');

for k=1:4
    subplot(2,2,k);
    xlabel('$\mathrm{Time} (t/t_{ch})$','Interpreter','Latex');
end


%% Save figure
% f2.Units='centimeters';
% f2.Position=[2,2,16,16];
% f2.PaperUnits='centimeters';
% f2.PaperSize=[16,16];
% print(f2,'figures/test_case_5A.pdf','-dpdf');