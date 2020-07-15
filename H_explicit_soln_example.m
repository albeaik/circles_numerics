function [ H_val ] = H_explicit_soln_example( evolving_mesh_obj, user_defined_coupler_obj, H_args_xv )
% DT ~ density surface mesh (currently implemented as delaunay triangulation)
% density ~ density value at inside each mesh element (triangle)
% H_args_xv ~ points at which to evaluate H function ~ each row is a 2d vector [x, v]

    
    DT = evolving_mesh_obj.DT;
    density = evolving_mesh_obj.density;
    
    %------------- | get mesh surface geometry information
    [DT_centroids] = evolving_mesh_obj.GetMeshCentroids(DT);  %sorted by triangle id --> can pre-compute if DT is fixed
    [DT_areas] = evolving_mesh_obj.GetMeshAreas(DT);          %sorted by triangle id --> can pre-compute if DT is fixed
    
    %------------- | 
    num_H_evaluation_points_t = size(H_args_xv, 1);
    num_DT_trigs_t = size(DT_centroids, 1);
    
    evalpts_argx_mtx = repelem(H_args_xv(:, 1)', num_DT_trigs_t, 1);                      %size = (integral domain ~ num autonomous cars, evaluation points ~ num characteristic samples)
    evalpts_argv_mtx = repelem(H_args_xv(:, 2)', num_DT_trigs_t, 1);                      %size = (integral domain ~ num autonomous cars, evaluation points ~ num characteristic samples)
    
    centroids_x_t_mtx = repelem(DT_centroids(:, 1), 1, num_H_evaluation_points_t);          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
    centroids_v_t_mtx = repelem(DT_centroids(:, 2), 1, num_H_evaluation_points_t);          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
    
    areas_t_mtx = repelem(DT_areas, 1, num_H_evaluation_points_t);                          %size = (integral domain ~ num trigs, evaluation points ~ num characteristic samples)
    
    active_nonlocal_domain = ((DT_centroids(:, 1) >= 0) & (DT_centroids(:, 1) <= 2)) & ((DT_centroids(:, 2) >= 0) & (DT_centroids(:, 2) <= 2));
    
    H_val = (active_nonlocal_domain' * (DT_areas .* density)) * ones([num_H_evaluation_points_t, 1]);
    
    