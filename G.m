function [ G_val ] = G( evolving_mesh_obj, user_defined_coupler_obj, H_args_xv )
%G Summary of this function goes here
%   Detailed explanation goes here

    DT = evolving_mesh_obj.DT;

    G_val = DT.Points(:, 2);
end

