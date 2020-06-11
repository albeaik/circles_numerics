function [ edge_lengths ] = GetDelaunayTriangleEdgeLengths( DT )
%GETDELAUNAYTRIANGLEEDGELENGTHS Summary of this function goes here
%   Detailed explanation goes here

    all_edges = [DT.ConnectivityList(:, [1, 2]); DT.ConnectivityList(:, [2, 3]); DT.ConnectivityList(:, [1, 3])];
	unique_edges = unique(sort(all_edges, 2), 'rows');

    edge_end_1 = DT.Points(unique_edges(:, 1), :);
    edge_end_2 = DT.Points(unique_edges(:, 2), :);
    
    edge_lengths = sqrt(sum((edge_end_2 - edge_end_1).^2, 2));

end

