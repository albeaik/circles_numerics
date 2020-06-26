function [ output_args ] = animate_solution( input_args )
%ANIMATE_SOLUTION Summary of this function goes here
%   Detailed explanation goes here


    % raw visualizations ---------------------------------------------------------


    i = 1;
    while(true)
        pause(0.1)
        i

        visualize_trig_trig( DT_history{i}, density_history{i} );
        title(['time - ', num2str(i/100)]);
        xlim(x_domain_lim)
        ylim(v_domain_lim)
        %caxis([0, max(max(cell2mat(density_history)))])
        caxis([0, 10])

        hold on
        plot3(Autonomous_state_history{i}(:, 1), Autonomous_state_history{i}(:, 2), 10*ones(size(Autonomous_state_history{i}(:, 2))), 'r*')
        hold off

        i = mod(i+1, size(DT_history, 2) -1 ) + 1;
    end


end

