classdef Dimensions
    %DIMENSIONS Stores information for non-dimensionalizing problem
    %   Detailed explanation goes here
    
    properties
        V0
        Z0
        u0
        Len0
        C0_L
        L0_L
        C0
        L0
        t0
        I0
    end
    
    methods
        function obj = Dimensions(V0,Z0,u0,L)
            %DIMENSIONS Construct an instance of this class
            %   Detailed explanation goes here
            if nargin==0
                V0=1e4; Z0=50; u0=3e8; L0=1;
            end
            obj.V0=V0;
            obj.Z0=Z0;
            obj.u0=u0;
            obj.Len0=L;
            obj.C0_L=1/(u0*Z0);
            obj.L0_L=Z0/u0;
            obj.C0=obj.C0_L*L;
            obj.L0=obj.L0_L*L;
            obj.t0=L/u0;
            obj.I0=V0/Z0;
        end
              
        function Rs = noDim_R(obj,R)
            Rs=R /obj.Z0;
        end
        function Cs = noDim_C(obj,C)
            Cs=C /obj.C0;
        end
        function Ls = noDim_len(obj,L)
            Ls=L /obj.Len0;
        end
        function ts = noDim_t(obj,t)
            ts=t /obj.t0;
        end
        function t = reDim_t(obj,ts)
            t=ts * obj.t0;
        end
        function E = reDim_E(obj,E)
            E=E*obj.V0*obj.I0*obj.t0;
        end
        
        
    end
end

