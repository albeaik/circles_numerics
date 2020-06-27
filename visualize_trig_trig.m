function visualize_trig_trig( DT_plot, density_plot )
%VISUALIZE_DOT_TRIG Summary of this function goes here
%   Detailed explanation goes here

    %ref: https://www.mathworks.com/matlabcentral/answers/165624-how-to-color-trisurf-faces
    
    trianglation_edge_alpha = 0;

    hh = trisurf(DT_plot.ConnectivityList, DT_plot.Points(:, 1), DT_plot.Points(:, 2), DT_plot.Points(:, 2)*0+1);
    set(gca,'CLim',[min(density_plot), max(density_plot)]);
    set(hh,'FaceColor','flat', 'FaceVertexCData',density_plot, 'CDataMapping','scaled');
    set(hh,'edgecolor', ([1 1 1])*0.3, 'EdgeAlpha', trianglation_edge_alpha)
    colorbar
    
    view(0,90)

end
