function [ total_surface_integral ] = GetTotalDensitySurfaceIntegral( DT, density )
%GETTOTALDENSITYSURFACEINTEGRAL Summary of this function goes here
%   Detailed explanation goes here

    total_surface_integral = density' * GetDelaunayAreas( DT );

end

