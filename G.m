function [ G_val ] = G( DT, density, Autonomous_state, H_args_xv, assume_fixed_DT )
%G Summary of this function goes here
%   Detailed explanation goes here

    G_val = DT.Points(:, 2);
end

