function visualize_dot_trig( DT_plot, density_plot, do_plot_trig )
%VISUALIZE_DOT_TRIG Summary of this function goes here
%   Detailed explanation goes here

    %figure
    hold on
    %i = 200;
    %DT_plot = DT_history{i};
    %density_plot = density_history{i};
    
    if(do_plot_trig)
        triplot(DT_plot.ConnectivityList, DT_plot.Points(:, 1), DT_plot.Points(:, 2))
    end

    color_resolution = 100;
    quantization_min = nanmin(density_plot);
    quantization_max = nanmax(density_plot);

    %cmap = jet(length(density_tau)); % Make 1000 colors.
    %cmap = gray(length(density_tau)); % Make 1000 colors.
    cmap = jet(color_resolution); % Make 1000 colors.

    quantizied_densities = floor(((density_plot - quantization_min) ./ quantization_max) * (color_resolution-1));
    plot_cmap = cmap(int32(quantizied_densities)+1, :);

    %scatter(xs, ys, 10, cmap, 'filled')

    centroids = GetDelaunayCentroids(DT_plot);
    scatter(centroids(:, 1), centroids(:, 2), 15, plot_cmap, 'filled')

    colormap(cmap)
    colorbar
    caxis([quantization_min, quantization_max])
    
    %set(gca,'Color','k')
    set(gca,'Color', ([1 1 1])*0.3)

end

