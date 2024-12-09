function [Vp] = load_pulse(filename, D)
%LOAD_PULSE Summary of this function goes here
%   Detailed explanation goes here

% load the experimental pulse
load(filename)
if exist('tp','var')
    t=tp; V=Vp;
end

% ensure continuity for extrapolation
V(logical((t>0).*(V<0)))=0;

% make signal start at t=0
tp=t-t(1);

% non-dimensionalize
tp=D.noDim_t(tp); V=V/D.V0;

% create continuous function
Vp=@(tq) interp1(tp,V,tq,'linear',0);

end

