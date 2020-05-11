function [ centroids, areas ] = GetDelaunayCentroids( DT )
%GETDELAUNAYCENTROIDS Summary of this function goes here
%   Detailed explanation goes here

    centroids = zeros([size(DT.ConnectivityList, 1), 2]);
    areas = zeros([size(DT.ConnectivityList, 1), 1]);
    parfor i = 1:size(DT.ConnectivityList, 1)
        trig = DT.Points(DT.ConnectivityList(i, :), :);
        trig_centroid = mean(trig);
        centroids(i, :) = trig_centroid;
        areas(i) = polyarea(trig(:, 1), trig(:, 2));
    end


end

