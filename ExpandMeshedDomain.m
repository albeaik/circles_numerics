function [ DT_new ] = ExpandMeshedDomain( DT_old, mesh_domain_limits, mesh_resolution )
%   expands the xy domain covered by DT_old to include mesh_domain_limits
%   NOTE: rough implementation! 
%   Suggested improvements:
%       1. don't remesh old domain
%       2. with current implementation, convex hull of DT_old is assumed well sampled

    DT_old = delaunayTriangulation(DT_old.Points);
    %figure, trimesh(DT_old.ConnectivityList, DT_old.Points(:, 1), DT_old.Points(:, 2))


    x=linspace(mesh_domain_limits(1, 1), mesh_domain_limits(1, 2), mesh_resolution(1));
    v=linspace(mesh_domain_limits(2, 1), mesh_domain_limits(2, 2), mesh_resolution(2));
    [x_grid, v_grid]=meshgrid(x, v);

    init_domain_points = [x_grid(:), v_grid(:)];

    old_boundary = freeBoundary(DT_old);
    %hold on, plot(DT_old.Points(old_boundary, 1), DT_old.Points(old_boundary, 2))
    bd = [old_boundary(:, 1); old_boundary(1, 1)];
    [in,on] = inpolygon(init_domain_points(:, 1), init_domain_points(:, 2), DT_old.Points(bd, 1), DT_old.Points(bd, 2));
    complement_domain_points = init_domain_points(~in, :);

    %bd = boundary_nodes(q_solution_mesh.DT.ConnectivityList);
    %[in,on] = inpolygon(init_domain_points(:, 1), init_domain_points(:, 2), DT_old.Points(bd, 1), DT_old.Points(bd, 2));

    complement_domain_points = init_domain_points(~in, :);
    %figure, plot(complement_domain_points(:, 1), complement_domain_points(:, 2), '.')

    %ID = pointLocation(DT_old, init_domain_points);
    %complement_domain_points = init_domain_points(isnan(ID), :);

    DT_new = DT_old;
    DT_new.Points = [DT_new.Points; complement_domain_points];
    %figure, trimesh(DT_new.ConnectivityList, DT_new.Points(:, 1), DT_new.Points(:, 2))



    %[bd(:, 1:en-1)
    %DT = delaunayTriangulation(P, C)


end

