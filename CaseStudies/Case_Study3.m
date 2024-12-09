clear; clc;

% establish characteristic dimensions
D=Dimensions(5e3,50,2.99792e8,10);
newcolors = [0.00,0.45,0.74; 1,0,0; 0.39,0.83,0.07; 1,0.41,0.16];
set(groot, "defaultaxescolororder", newcolors);

%% Construct system (USER INPUT HERE)
% FOR COLIN: Everything in this user input section has been modified to
% what we want. Idea is that we want L1 = 20m, L2 = 100m and probes only at x =
% 10, 30 m (equidistant to load on HV and LV side, nothing at the load).
% With the current waveforms, we want to have the same graphs you did
% before with the incident, transmitted, reflected, and the difference powers. 

% Wire information
Lc=D.noDim_len([20,100]); % cable lengths
Rc=D.noDim_R([50 50]); % cable impedances
u=[0.67 0.67]; % cable wave speeds
dx=D.noDim_len([5e-2,5e-2]); % taget grid spacing
% simulation information
tmax=D.noDim_t(500e-9);
%tmax=5;
%pulse_file='../../Data/Ref_Waveforms/Pulse_300VDC.mat';
%Vp =@(t) (exp(-(t-1.5).^2/0.2));
Vp = @(t) exp(-(t-35./33.4).^2./0.1291); % new gaussian pulse

for i=1:4
    switch i
        case 1
            RL=[10]; % load resistances
            CL=[0.1]; % load capacitances
            name='Load A';
        case 2
            RL=[10]; % load resistances
            CL=[0.01]; % load capacitances
            name='Load B';
        case 3
            RL=[0.5]; % load resistances
            CL=[0.1]; % load capacitances
            name='Load C';
        case 4
            RL=[0.5]; % load resistances
            CL=[0.01]; % load capacitances
            name='Load D';
    end
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
    f2=figure(2);
    if i==1
        clf;
        leg={};
    end
    Nmid1=round(WS.Nx(1)/2); % should be x = 10m
    Nprobe2 = WS.Nx(1) + Nmid1; % should be 20m + 10m, so 30m
    Ip1=sol.y(N+Nmid1,:);
    Ip2=sol.y(N+Nprobe2,:);
    VL=sol.y(sum(WS.Nx(1)),:)-sol.y(sum(WS.Nx(1))+1,:);
    IL=sol.y(N+sum(WS.Nx(1)),:);
    t=sol.x;
    
    subplot(3,2,1); hold on
    p1=plot(t,Ip1);
    title('BCS HV');
    subplot(3,2,2); hold on
    p2=plot(t,Ip2);
    title('BCS LV');
    subplot(3,2,2+i); hold on

    % FOR COLIN: Here's where we got very confused. We understand you're
    % doing some windowing in non-dimensional time (different ND time to t_ch)
    % to separate the incident and reflected pulses but we couldn't find good 
    % numbers to use to replace 1.5, 2.5 and so on. That is what we'd like
    % help with! I can fine tune the plots at the end, but we just want to
    % resolve the incident and reflected well.
    
    % For Sankarsh: the 1.5 is just to place it in the middle of the domain
    % for plotting. The only value that needs to be calibrated is "tflip"
    % symmetric placement of the probes means its the same for both current
    % shunts. If they were not the same distance from the load, we would
    % need different values for each shunt
    % if you want to move the entire graph left/right, change the 1.5
    tflip=Lc(1)/u(1);
    
    rng=1:length(t);
    rng1=rng(t<=1.5.*tflip);
    rng2=rng(logical((t>tflip) .* (t<=2.5*tflip)));
    
    tgrid=linspace(0,1.5*tflip,201);
    Ip_inc=interp1(t(rng1),Ip1(rng1),tgrid);
    Ip_trans=interp1(t(rng2)-tflip,Ip2(rng2),tgrid);
    Ip_refl=interp1(t(rng2)-tflip,Ip1(rng2),tgrid);
    PL=cumtrapz(t,VL.^2/RL(1));
    
    
    plot(tgrid-1.5,Ip_inc.^2);
    plot(tgrid-1.5,Ip_trans.^2);
    plot(tgrid-1.5,Ip_refl.^2, "Color", [1,0.41,0.16]);
    Pnet=Ip_inc.^2-Ip_trans.^2-Ip_refl.^2;
    plot(tgrid-1.5,Pnet, '--k');
    yyaxis right
    rng3=logical(tgrid>=1.5);
    plot(tgrid(rng3)-1.5, cumtrapz(tgrid(rng3),Pnet(rng3)), 'g');
    plot(t-1.5-tflip/2,PL, '--g');
    %title(sprintf('Load %d',i));
    xlim([0 2]);
    legend('Incident','Transmitted','Reflected', 'Difference', ...
    'location','NorthWest');
    %plot(t,cumtrapz(t,Ip.*Vpr));
    %ylabel('$E/E_0$', 'Interpreter','Latex');
    leg{i}=name;
    %
    if i==4
        for k=1:6
            subplot(3,2,k)
            xlabel('Time ($t/t_{ch}$)','Interpreter','Latex');
            if k<3
                ylabel('$I/I_{ch}$', 'Interpreter','Latex');
                legend(leg{:});
            else
                yyaxis left
                ylabel('$P/P_{ch}$', 'Interpreter','Latex');
                ylim([-1,1.5])
                yyaxis right
                ylabel('$E/E_{ch}$', 'Interpreter','Latex');
                ylim([-0.1,0.3]);
                
            end
        end
    end
end
subplot(3,2,3); title('Load A'); subplot(3,2,4); title('Load B');
subplot(3,2,5); title('Load C'); subplot(3,2,6); title('Load D');
% 
% 
%% Save figure
%f2.Units='centimeters';
%f2.Position=[2,2,26,16];
%f2.PaperUnits='centimeters';
%f2.PaperSize=[26,16];
%print(f2,'figures/Case_Study3.pdf','-dpdf');