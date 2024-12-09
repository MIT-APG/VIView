function [M] = Mass(t,state,WS)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    Ntot=sum(WS.Nx);
    M=speye(2*Ntot);
    % row 1 - specified voltage
    M(1,1)=0;
    % Ground connection
    M(Ntot,:)=0;
    
    for i=1:WS.Num_Load
        % Current Continuity
        row=sum(WS.Nx(1:i))+Ntot+1;
        M(row,row)=0;
        % Load Condition
        row=sum(WS.Nx(1:i));
        M(row,row)=WS.RL(i)*WS.CL(i); % over-write anything here
        M(row,row+1)=-WS.RL(i)*WS.CL(i); % coefficient on current
    end
end

