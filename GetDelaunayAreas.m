function [ areas ] = GetDelaunayAreas( DT )
%GETDELAUNAYAREAS Summary of this function goes here
%   Detailed explanation goes here

    trigs_x = reshape(DT.Points(DT.ConnectivityList', 1), [3, size(DT.ConnectivityList, 1)])';   %size = (num trigs, num vert per trig)
    trigs_y = reshape(DT.Points(DT.ConnectivityList', 2), [3, size(DT.ConnectivityList, 1)])';   %size = (num trigs, num vert per trig)
    
    areas = polyarea(trigs_x', trigs_y')';

end

