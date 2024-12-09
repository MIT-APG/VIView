function [M] = Mass_Spark_Cheat(t,state,WS,spark_param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    Rspark=spark_param(1);
    tspark=spark_param(2);

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
        if i==2
            R=10^(log10(WS.RL(i))*0.5*(tanh(-(t-tspark)*20)+1)+...
                log10(Rspark)*0.5*(tanh((t-tspark)*20)+1));% coefficient on current
            if R<1
                flag=1;
            end
        else
            R=WS.RL(i);
        end
        row=sum(WS.Nx(1:i));
        M(row,row)=R*WS.CL(i); % over-write anything here
        M(row,row+1)=-R*WS.CL(i); % coefficient on current
    end
end