function [value, isterminal, direction] = spark_event(t,y,WS, load_index, thresh, dir)
%SPARK_EVENT Summary of this function goes here
%   Detailed explanation goes here
dV=y(sum(WS.Nx(1:load_index)))-y(sum(WS.Nx(1:load_index))+1);
value=dV-thresh;
isterminal = 1;
direction = dir;

end