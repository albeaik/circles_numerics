function [ H_val ] = H_IDM( evolving_mesh_obj, user_defined_coupler_obj, H_args_xv )
% DT ~ density surface mesh (currently implemented as delaunay triangulation)
% density ~ density value at inside each mesh element (triangle)
% H_args_xv ~ points at which to evaluate H function ~ each row is a 2d vector [x, v]
% assume_fixed_DT ~= 0 ==> triangulation is assumed not changing between
%                           iterations ==> this permits pre-computation of
%                           integral domains and relevant matrices
%                      ==> NOT IMPLEMENTED
    
    assume_fixed_DT = 0;
    
    DT = evolving_mesh_obj.DT;
    density = evolving_mesh_obj.density;
    
    Autonomous_state = user_defined_coupler_obj.Autonomous_state;
    activate_coupling = user_defined_coupler_obj.activate_coupling;

    %------------- | parameters
    epsilon_nut = 1;
    d_nut = 2.5;
    V_max = 30;         %same as v_o for IDM??!
    alpha = 1;
    bet
    20;
    
    %H_IDM parameters       
    v_o = V_max;        %the velocity the vehicle would drive at in free traffic
    little_delta = 4;   
    s_o = 2;            %minimum spacings0:  a minimum desired net distance.  A car can not move if the distance from thecar in front is not at leasts0
    T = 1;              %desired time headwayT:  the minimum possible time to the vehicle in front (headway is the distancebetween vehicles in a transit system measured in time or space; the minimum headway is the shortestsuch distance or time achieved by a system without a reduction in the speed of vehicles)
    a = 1.3;              %accelerationa:  the maximum vehicle acceleration
    b = 2;              %comfortable braking decelerationb:  a positive number
    l = 4;              %gives the length of the vehicle
    
    
    %------------- | get mesh surface geometry information
    [DT_centroids] = evolving_mesh_obj.GetMeshCentroids(DT);  %sorted by triangle id --> can pre-compute if DT is fixed
    [DT_areas] = evolving_mesh_obj.GetMeshAreas(DT);          %sorted by triangle id --> can pre-compute if DT is fixed
    
    
%     %------------- | process and evaluate IDM effects
%     num_H_evaluation_points_t = size(H_args_xv, 1);
%     num_DT_trigs_t = size(DT_centroids, 1);
% 
%     areas_t_mtx = repelem(DT_areas, 1, num_H_evaluation_points_t);                          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
% 
%     centroids_x_t_mtx = repelem(DT_centroids(:, 1), 1, num_H_evaluation_points_t);          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
%     centroids_v_t_mtx = repelem(DT_centroids(:, 2), 1, num_H_evaluation_points_t);          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
% 
%     DT_x_t_mtx = repelem(H_args_xv(:, 1)', num_DT_trigs_t, 1);                 %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
%     DT_v_t_mtx = repelem(H_args_xv(:, 2)', num_DT_trigs_t, 1);                 %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
% 
% 
%     H1_integral_domain = ( ((DT_x_t_mtx < centroids_x_t_mtx) & (centroids_x_t_mtx < DT_x_t_mtx+epsilon_nut)) ) & ( ((0 < centroids_v_t_mtx) & (centroids_v_t_mtx < max(DT.Points(:, 2)))) );
% 
% 
%     y_mtx = centroids_x_t_mtx; %(H1_integral_domain);       %note: already incorporated into H1_integral_domain
%     v_telda_mtx = centroids_v_t_mtx; %(H1_integral_domain); %note: already incorporated into H1_integral_domain
%     integral_area_mtx = areas_t_mtx; %(H1_integral_domain);
% 
%     pre_computed_h_mtx = h(DT_x_t_mtx - y_mtx, epsilon_nut);
%     pre_computed_Vterm_mtx = (V(y_mtx - DT_x_t_mtx, d_nut, V_max) - DT_v_t_mtx);
% 
%     sparse_h_V_area_terms_H1 = sparse(H1_integral_domain .* pre_computed_h_mtx .* pre_computed_Vterm_mtx .* integral_area_mtx);
%     sparse_h_V_area_terms_H2 = sparse(pre_computed_h_mtx .* (DT_v_t_mtx - v_telda_mtx) .* integral_area_mtx);
%     %-------------
%     
%     %-------------
%     H1_val = alpha .* sparse_h_V_area_terms_H1' * density;
%     H2_val = - beta .* sparse_h_V_area_terms_H2' * density; %confirm this is supposed to be negative beta
%     %-------------
    
	%------------- | process and evaluate IDM effects
    num_H_evaluation_points_t = size(H_args_xv, 1);
    num_DT_trigs_t = size(DT_centroids, 1);

    areas_t_mtx = repelem(DT_areas, 1, num_H_evaluation_points_t);                          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)

    centroids_x_t_mtx = repelem(DT_centroids(:, 1), 1, num_H_evaluation_points_t);          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
    centroids_v_t_mtx = repelem(DT_centroids(:, 2), 1, num_H_evaluation_points_t);          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)

    DT_x_t_mtx = repelem(H_args_xv(:, 1)', num_DT_trigs_t, 1);                 %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
    DT_v_t_mtx = repelem(H_args_xv(:, 2)', num_DT_trigs_t, 1);                 %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)


    %H1_integral_domain = ( ((DT_x_t_mtx < centroids_x_t_mtx) & (centroids_x_t_mtx < DT_x_t_mtx+epsilon_nut)) ) & ( ((0 < centroids_v_t_mtx) & (centroids_v_t_mtx < max(DT.Points(:, 2)))) );


    y_mtx = centroids_x_t_mtx; %(H1_integral_domain);       %note: already incorporated into H1_integral_domain
    v_telda_mtx = centroids_v_t_mtx; %(H1_integral_domain); %note: already incorporated into H1_integral_domain
    integral_area_mtx = areas_t_mtx; %(H1_integral_domain);

    pre_computed_h_mtx = h(DT_x_t_mtx - y_mtx, epsilon_nut);
    %pre_computed_Vterm_mtx = (V(y_mtx - DT_x_t_mtx, d_nut, V_max) - DT_v_t_mtx);

    %sparse_h_V_area_terms_H1 = sparse(H1_integral_domain .* pre_computed_h_mtx .* pre_computed_Vterm_mtx .* integral_area_mtx);
    %sparse_h_V_area_terms_H2 = sparse(pre_computed_h_mtx .* (DT_v_t_mtx - v_telda_mtx) .* integral_area_mtx);
    %-------------
    
    %-------------
    H_IDM_evalpts_density_mtx = repelem(density, 1, num_H_evaluation_points_t);       %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
    
	H_IDM_effect = a - a.* (H_args_xv(:, 2)./v_o).^little_delta - a.* sum(pre_computed_h_mtx .* H_IDM_evalpts_density_mtx .* ( (s_o + DT_v_t_mtx .* T + DT_v_t_mtx .* (DT_v_t_mtx - v_telda_mtx) ./ sqrt(2.*a.*b) ) ./ (y_mtx - DT_x_t_mtx - l) ).^2 .* integral_area_mtx)';
    
    H_IDM_effect(isnan(H_IDM_effect)) = 0;  %WARNING!! hard coded fix to numerical explosions caused by div by zero.. it should be ok for now, BUT fix later
    %-------------
    
    %------------- | process and evaluate autonomous vehicle coupling effects
    num_aut_cars = size(Autonomous_state, 1);
    aut_evalpts_argx_mtx = repelem(H_args_xv(:, 1)', num_aut_cars, 1);                      %size = (integral domain ~ num autonomous cars, evaluation points ~ num characteristic samples)
    aut_evalpts_argv_mtx = repelem(H_args_xv(:, 2)', num_aut_cars, 1);                      %size = (integral domain ~ num autonomous cars, evaluation points ~ num characteristic samples)
    aut_evalpts_auty_mtx = repelem(Autonomous_state(:, 1), 1, num_H_evaluation_points_t);   %size = (integral domain ~ num autonomous cars, evaluation points ~ num characteristic samples)
    aut_evalpts_autw_mtx = repelem(Autonomous_state(:, 2), 1, num_H_evaluation_points_t);   %size = (integral domain ~ num autonomous cars, evaluation points ~ num characteristic samples)
    %-------------
    
    %-------------
    H1a_val = alpha .* mean(h(aut_evalpts_argx_mtx - aut_evalpts_auty_mtx, epsilon_nut) .* (V(aut_evalpts_auty_mtx - aut_evalpts_argx_mtx, d_nut, V_max) - aut_evalpts_argv_mtx))';
    H2a_val = beta .* mean(h(aut_evalpts_argx_mtx - aut_evalpts_auty_mtx, epsilon_nut) .* (aut_evalpts_autw_mtx - aut_evalpts_argv_mtx))';
    %-------------
    
    
    %------------- | evaluate total effect
    %H_val = H1_val + H2_val + H1a_val + H2a_val;
    H_val = H_IDM_effect + (H1a_val + H2a_val) .* activate_coupling;
    %-------------



function [h_val] = h(x, epsilon_nut)
    h_val = exp(-1./( (epsilon_nut/2)^2 - (-x-epsilon_nut/2).^2 )) .* ((-epsilon_nut < x') & (x' < 0))';
    h_val(isnan(h_val)) = 0;


function [V_val] = V(x, d_nut, V_max)
    V_val = V_max * ((tanh(x/d_nut - 2) + tanh(2)) / (1+tanh(2)));


