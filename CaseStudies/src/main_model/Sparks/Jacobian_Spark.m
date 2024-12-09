function [J] = Jacobian_Spark(t,state,WS)
%JACOBIAN Summary of this function goes here
%   Detailed explanation goes here

%% diagonal constructors for 2nd order FD with first order on boundaries
Ntot=sum(WS.Nx);
is=1;
main=zeros(Ntot,1);
up=main;
down=main;
nil=main;
for i=1:length(WS.Nx)
    ie=is+WS.Nx(i)-1;   
    dx=WS.dx(i);
    main(is:ie)=zeros(WS.Nx(i),1);
    up(is:ie)=1./([1;dx;2*dx.*ones(WS.Nx(i)-2,1)]);
    down(is:ie)=-1./([2*dx.*ones(WS.Nx(i)-2,1);dx;1]);
    main(is)=-1/dx;
    main(ie)=1/dx;
    up(is)=0;
    down(ie)=0;
    is=ie+1;
end
% differential operator matrix
Ddiag=[[main;main],[up;up],[down;down]];
MD=spdiags(Ddiag,[0,1,-1],2*Ntot,2*Ntot);

%% Coefficient Matrix
% coefficient matrix
Aup=zeros(Ntot,1);
Adown=zeros(Ntot,1);
is=1;
for i=1:length(WS.Nx)
    ie=is+WS.Nx(i)-1;   
    u=WS.u(i);
    R=WS.Rc(i);
    Aup(is:ie)=u*R*ones(WS.Nx(i),1);
    Adown(is:ie)=u/R*ones(WS.Nx(i),1);
    is=ie+1;
end
Adiag=-[[nil;Aup],...
    [Adown;nil]];
MA=spdiags(Adiag,[Ntot,-Ntot],2*Ntot,2*Ntot);

%% Construct Full Matrix (no BCs)
J=MA*MD;

%% Now add in the BCs
% these require a mass matrix formulation
% Voltage at x=0 is prescibred (MAY UPDATE THIS FOR SPARK GAP)
J(1,:)=0;
J(1,1)=1;

% Contact with ground is 0 Voltage prescribed
J(Ntot,:)=0;
J(Ntot,Ntot)=1;

for i=1:WS.Num_Load
    % always have current continuity
    row=Ntot+sum(WS.Nx(1:i))+1;
    J(row,:)=0; % over-write anything here
    J(row,row-1)=1;
    J(row,row)=-1;
    
    % now apply load condition
    row=sum(WS.Nx(1:i));
    % calculate voltage across load element
    dV=state(row)-state(row+1);
    RL=WS.RL(dV);
    if RL(end)<50
        fprintf('Load Sparks at t=%.2f with dV=%.2f\n', t,dV(end));
    end
    J(row,:)=0; % over-write anything here
    J(row,row+Ntot)=RL(i); % coefficient on current
    J(row,row)=-1; % coefficient on upstream voltage
    J(row,row+1)=1; % coefficient on downstream voltage
end

end
