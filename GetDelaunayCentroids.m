function [ centroids, areas ] = GetDelaunayCentroids( DT )
%GETDELAUNAYCENTROIDS Summary of this function goes here
%   Detailed explanation goes here

    trigs_x = reshape(DT.Points(DT.ConnectivityList', 1), [3, size(DT.ConnectivityList, 1)])';   %size = (num trigs, num vert per trig)
    trigs_y = reshape(DT.Points(DT.ConnectivityList', 2), [3, size(DT.ConnectivityList, 1)])';   %size = (num trigs, num vert per trig)
    
    centroids = [mean(trigs_x, 2), mean(trigs_y, 2)];
    areas = polyarea(trigs_x', trigs_y')';

end

