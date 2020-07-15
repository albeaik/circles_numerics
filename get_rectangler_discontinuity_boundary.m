function [ points, edges ] = get_rectangler_discontinuity_boundary( lower_left_corner, upper_right_corner, res_step )

    if(length(res_step) == 1)
        resx = res_step;
        resy = res_step;
    else
        resx = res_step(1);
        resy = res_step(2);
    end

    lower_right_corner = [upper_right_corner(1), lower_left_corner(2)];
    upper_left_corner = [lower_left_corner(1), upper_right_corner(2)];
    
    
    tmp1 = [lower_left_corner(1):resx:lower_right_corner(1)]';
    lower_edge_points = [tmp1, lower_left_corner(2)*ones(size(tmp1))];
    
    tmp1 = [lower_right_corner(2):resy:upper_right_corner(2)]';
    right_edge_points = [lower_right_corner(1)*ones(size(tmp1)), tmp1];
    
    tmp1 = [upper_right_corner(1):-resx:upper_left_corner(1)]';
    upper_edge_points = [tmp1, upper_left_corner(2)*ones(size(tmp1))];
    
    tmp1 = [upper_left_corner(2):-resy:lower_left_corner(2)]';
    left_edge_points = [lower_left_corner(1)*ones(size(tmp1)), tmp1];
    
    points = [lower_edge_points; right_edge_points; upper_edge_points; left_edge_points];
    edges = [[1:size(points, 1)]', [2:size(points, 1), 1]'];
end

