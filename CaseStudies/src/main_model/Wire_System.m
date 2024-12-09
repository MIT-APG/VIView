classdef Wire_System
    %SYSTEM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Num_Cab % number of cable elements
        Num_Load % number of load elements = Num_Cab-1
        Lc % length of each cable element
        Rc % impedance of each cable element
        u % wave speed of each cable element
        RL % resistance of each load element
        CL % capacitance of each load element   
        dx % grid spacing
        x % grid points
        Nx % number of grid points per cable
    end
    
    methods
        function obj = Wire_System(Lc, Rc, u, RL, CL, dxt)
            %SYSTEM Construct an instance of this class
            %   Detailed explanation goes here
            obj.Num_Cab=length(Lc);
            obj.Num_Load=obj.Num_Cab-1;
            if length(Rc)==1
                Rc=Rc*ones(1,obj.Num_Cab);
            end
            if length(u)==1
                u=u*ones(1,obj.Num_Cab);
            end
            obj.Lc=Lc;
            obj.Rc=Rc;
            obj.u=u;
            obj.RL=RL;
            obj.CL=CL;
            [obj.dx,obj.x, obj.Nx]=make_grid(obj,dxt);            
        end
        
        function [dx,x, N]= make_grid(obj,dxt)
            if length(dxt)==1
                dxt=dxt*ones(1,obj.Num_Cab);
            end
            dx=zeros(1,length(dxt));
            N=zeros(1,length(dxt));
            for i=1:length(dx)
                N(i)=ceil(obj.Lc(i)/dxt(i));
                dx(i)=obj.Lc(i)/N(i);
            end
            x=zeros(1,sum(N));
            cnt=1;
            for i=1:length(dx)
                if cnt==1
                    offset=0;
                else
                    offset=x(cnt-1);
                end
                x(cnt:cnt+N(i))=(0:dx(i):obj.Lc(i))+offset;
                cnt=cnt+N(i)+1;
            end 
            N=N + 1;
        end

              
    end
end

