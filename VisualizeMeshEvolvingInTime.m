function VisualizeMeshEvolvingInTime(solution, coupling, plt_times)
% plt_times is list of time instances to plot
% ex: plt_times = [0, 1];
    
    dt = solution.time_history{2}-solution.time_history{1};     %assume constant dt
    fighandle = figure
    %plt_times = [0, 1];
    for t_ind = 1:size(plt_times, 2)
        subplot(2, size(plt_times, 2), t_ind)
        time = plt_times(t_ind);
        n = plt_times(t_ind)/dt+1; %find(cell2mat(solutions{index}.q.time_history) == time);
        solution.VisualizeStep(n);
        coupling.VisualizeStep(n);

        xlimvals = [0, 18];
        ylimvals = [5, 16];
        xlim(xlimvals)
        ylim(ylimvals)
        caxis([0, 2]);
        colorbar off;

        if(t_ind == 1)
            ylabel('Density')
        end
        set(gca,'FontSize',15);

        %---------------
        subplot(2, size(plt_times, 2), size(plt_times, 2)+t_ind)
        time = plt_times(t_ind);
        n = plt_times(t_ind)/dt+1; %find(cell2mat(solutions{index}.q.time_history) == time);
        %solution.VisualizeStep(n);
        %autonomous_car_coupling.VisualizeStep(n);
        triplot(solution.DT_history{n}.ConnectivityList, solution.DT_history{n}.Points(:, 1), solution.DT_history{n}.Points(:, 2));

        xlim(xlimvals)
        ylim(ylimvals)
        caxis([0, 2]);
        colorbar off;

        if(t_ind == 1)
            ylabel('Mesh')
        end
        set(gca,'FontSize',15);
    end
    set(gcf, 'Position',  [74         614        1562         334])
    set(gcf, 'Position',  [1         465        1680         483])
end