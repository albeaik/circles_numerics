
%save('../simulations/june_12_presentation/sim_res80_u_n2_n2_long_run_longer', 'DT_history', 'density_history', 'Autonomous_state_history')
%save('../simulations/june_12_presentation/sim_res80_u_n10_n10', 'DT_history', 'density_history', 'Autonomous_state_history')

clear all
close all
clc

x_domain_lim = [0, 20];
%v_domain_lim = [0, 20];
v_domain_lim = [4, 16];

load('../simulations/june_12_presentation/sim_res80_noaut')
density_hist_res80_noaut = density_history;
DT_hist_res80_noaut = DT_history;

load('../simulations/june_12_presentation/sim_res80_u_n10_n10')
density_hist_res80_u_n10_n10 = density_history;
DT_hist_res80_u_n10_n10 = DT_history;
Autonomous_state_history_res80_u_n10_n10 = Autonomous_state_history;

load('../simulations/june_12_presentation/sim_res40_u_n10_n10')
density_hist_res40_u_n10_n10 = density_history;
DT_hist_res40_u_n10_n10 = DT_history;
Autonomous_state_history_res40_u_n10_n10 = Autonomous_state_history;

load('../simulations/june_12_presentation/sim_res20_u_n10_n10')
density_hist_res20_u_n10_n10 = density_history;
DT_hist_res20_u_n10_n10 = DT_history;
Autonomous_state_history_res20_u_n10_n10 = Autonomous_state_history;

load('../simulations/june_12_presentation/sim_res80_u_n2_n2_long_run')
%load('../simulations/june_12_presentation/sim_res80_u_n2_n2_long_run_longer')
density_hist_res80_u_n2_n2_long = density_history;
DT_hist_res80_u_n2_n2_long = DT_history;
Autonomous_state_history_res80_u_n2_n2_long = Autonomous_state_history;


h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = '../plots/june_12_presentation/char_sim.gif';
for n = 1:size(DT_hist_res80_u_n2_n2_long, 2)
    % Draw plot for y = x.^n
    
    %~~~~~~~~ | draw --> shouldn't create side effects after this line
    subplot(3, 2, [1,2])
    density_history = density_hist_res80_u_n2_n2_long;
    DT_history = DT_hist_res80_u_n2_n2_long;
    Autonomous_state_history = Autonomous_state_history_res80_u_n2_n2_long;
    i = mod(n - 1, size(DT_history, 2) ) + 1;
    visualize_trig_trig( DT_history{i}, density_history{i} );
    title(['u = [-2, -2] - res80 - time - ', num2str(i/100)]);
    xlim(x_domain_lim*2)
    ylim(v_domain_lim)
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 5])
    colormap('jet')
    
    hold on
    plot3(Autonomous_state_history{i}(:, 1), Autonomous_state_history{i}(:, 2), 10*ones(size(Autonomous_state_history{i}(:, 2))), 'r*')
    hold off
    
    %~~~~~~~~ | draw --> shouldn't create side effects after this line
    subplot(3, 2, 3)
    density_history = density_hist_res80_noaut;
    DT_history = DT_hist_res80_noaut;
    i = mod(n - 1, size(DT_history, 2) ) + 1;
    visualize_trig_trig( DT_history{i}, density_history{i} );
    title(['no autonomous - res80 - time - ', num2str(i/100)]);
    xlim(x_domain_lim)
    ylim(v_domain_lim)
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 2])
    colormap('jet')
    
    hold on
    %plot3(Autonomous_state_history{i}(:, 1), Autonomous_state_history{i}(:, 2), 10*ones(size(Autonomous_state_history{i}(:, 2))), 'r*')
    hold off
    
    %~~~~~~~~ | draw --> shouldn't create side effects after this line
    subplot(3, 2, 4)
    density_history = density_hist_res80_u_n10_n10;
    DT_history = DT_hist_res80_u_n10_n10;
    Autonomous_state_history = Autonomous_state_history_res80_u_n10_n10;
    i = mod(n - 1, size(DT_history, 2) ) + 1;
    visualize_trig_trig( DT_history{i}, density_history{i} );
    title(['u = [-10, -10] - res80 - time - ', num2str(i/100)]);
    xlim(x_domain_lim)
    ylim(v_domain_lim)
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 2])
    colormap('jet')
    
    hold on
    plot3(Autonomous_state_history{i}(:, 1), Autonomous_state_history{i}(:, 2), 10*ones(size(Autonomous_state_history{i}(:, 2))), 'r*')
    hold off
    
    %~~~~~~~~ | draw --> shouldn't create side effects after this line
    subplot(3, 2, 5)
    density_history = density_hist_res40_u_n10_n10;
    Autonomous_state_history = Autonomous_state_history_res40_u_n10_n10;
    DT_history = DT_hist_res40_u_n10_n10;
    i = mod(n - 1, size(DT_history, 2) ) + 1;
    visualize_trig_trig( DT_history{i}, density_history{i} );
    title(['u = [-10, -10] - res40 - time - ', num2str(i/100)]);
    xlim(x_domain_lim)
    ylim(v_domain_lim)
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 2])
    colormap('jet')
    
    hold on
    plot3(Autonomous_state_history{i}(:, 1), Autonomous_state_history{i}(:, 2), 10*ones(size(Autonomous_state_history{i}(:, 2))), 'r*')
    hold off
    
    %~~~~~~~~ | draw --> shouldn't create side effects after this line
    subplot(3, 2, 6)
    density_history = density_hist_res20_u_n10_n10;
    Autonomous_state_history = Autonomous_state_history_res20_u_n10_n10;
    DT_history = DT_hist_res20_u_n10_n10;
    i = mod(n - 1, size(DT_history, 2) ) + 1;
    visualize_trig_trig( DT_history{i}, density_history{i} );
    title(['u = [-10, -10] - res20 - time - ', num2str(i/100)]);
    xlim(x_domain_lim)
    ylim(v_domain_lim)
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 2])
    colormap('jet')
    
    hold on
    plot3(Autonomous_state_history{i}(:, 1), Autonomous_state_history{i}(:, 2), 10*ones(size(Autonomous_state_history{i}(:, 2))), 'r*')
    hold off
    
    %~~~~~~~~
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
  
h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = '../plots/june_12_presentation/char_sim_res80_n2_n2_long.gif';
for n = 1:size(DT_hist_res80_u_n2_n2_long, 2)
    % Draw plot for y = x.^n
    
    %~~~~~~~~ | draw --> shouldn't create side effects after this line
    density_history = density_hist_res80_u_n2_n2_long;
    DT_history = DT_hist_res80_u_n2_n2_long;
    Autonomous_state_history = Autonomous_state_history_res80_u_n2_n2_long;
    i = mod(n - 1, size(DT_history, 2) ) + 1;
    visualize_trig_trig( DT_history{i}, density_history{i} );
    title(['u = [-2, -2] - res80 - time - ', num2str(i/100)]);
    xlim(x_domain_lim*2)
    ylim(v_domain_lim)
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 5])
    colormap('jet')
    
    hold on
    plot3(Autonomous_state_history{i}(:, 1), Autonomous_state_history{i}(:, 2), 10*ones(size(Autonomous_state_history{i}(:, 2))), 'r*')
    hold off
    
    %~~~~~~~~
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
  



h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = '../plots/june_12_presentation/char_sim_res80_noaut.gif';
for n = 1:size(DT_hist_res80_noaut, 2)
    % Draw plot for y = x.^n
    
    %~~~~~~~~ | draw --> shouldn't create side effects after this line
    density_history = density_hist_res80_noaut;
    DT_history = DT_hist_res80_noaut;
    i = mod(n - 1, size(DT_history, 2) ) + 1;
    visualize_trig_trig( DT_history{i}, density_history{i} );
    title(['no autonomous - res80 - time - ', num2str(i/100)]);
    xlim(x_domain_lim)
    ylim(v_domain_lim)
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 2])
    colormap('jet')
    
    hold on
    %plot3(Autonomous_state_history{i}(:, 1), Autonomous_state_history{i}(:, 2), 10*ones(size(Autonomous_state_history{i}(:, 2))), 'r*')
    hold off
    
    %~~~~~~~~
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
  




h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = '../plots/june_12_presentation/char_sim_res80_u_n10_n10.gif';
for n = 1:size(DT_hist_res80_u_n10_n10, 2)
    % Draw plot for y = x.^n
    
    %~~~~~~~~ | draw --> shouldn't create side effects after this line
    density_history = density_hist_res80_u_n10_n10;
    DT_history = DT_hist_res80_u_n10_n10;
    Autonomous_state_history = Autonomous_state_history_res80_u_n10_n10;
    i = mod(n - 1, size(DT_history, 2) ) + 1;
    visualize_trig_trig( DT_history{i}, density_history{i} );
    title(['u = [-10, -10] - res80 - time - ', num2str(i/100)]);
    xlim(x_domain_lim)
    ylim(v_domain_lim)
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 2])
    colormap('jet')
    
    hold on
    plot3(Autonomous_state_history{i}(:, 1), Autonomous_state_history{i}(:, 2), 10*ones(size(Autonomous_state_history{i}(:, 2))), 'r*')
    hold off
    
    %~~~~~~~~
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
  



h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = '../plots/june_12_presentation/char_sim_res40_u_n10_n10.gif';
for n = 1:size(DT_hist_res40_u_n10_n10, 2)
    % Draw plot for y = x.^n
    
    %~~~~~~~~ | draw --> shouldn't create side effects after this line
    density_history = density_hist_res40_u_n10_n10;
    DT_history = DT_hist_res40_u_n10_n10;
    Autonomous_state_history = Autonomous_state_history_res40_u_n10_n10;
    i = mod(n - 1, size(DT_history, 2) ) + 1;
    visualize_trig_trig( DT_history{i}, density_history{i} );
    title(['u = [-10, -10] - res40 - time - ', num2str(i/100)]);
    xlim(x_domain_lim)
    ylim(v_domain_lim)
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 2])
    colormap('jet')
    
    hold on
    plot3(Autonomous_state_history{i}(:, 1), Autonomous_state_history{i}(:, 2), 10*ones(size(Autonomous_state_history{i}(:, 2))), 'r*')
    hold off
    
    %~~~~~~~~
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end
  



h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = '../plots/june_12_presentation/char_sim_res20_u_n10_n10.gif';
for n = 1:size(DT_hist_res20_u_n10_n10, 2)
    % Draw plot for y = x.^n
    
    %~~~~~~~~ | draw --> shouldn't create side effects after this line
    density_history = density_hist_res20_u_n10_n10;
    DT_history = DT_hist_res20_u_n10_n10;
    Autonomous_state_history = Autonomous_state_history_res20_u_n10_n10;
    i = mod(n - 1, size(DT_history, 2) ) + 1;
    visualize_trig_trig( DT_history{i}, density_history{i} );
    title(['u = [-10, -10] - res20 - time - ', num2str(i/100)]);
    xlim(x_domain_lim)
    ylim(v_domain_lim)
    %caxis([0, max(max(cell2mat(density_history)))])
    caxis([0, 2])
    colormap('jet')
    
    hold on
    plot3(Autonomous_state_history{i}(:, 1), Autonomous_state_history{i}(:, 2), 10*ones(size(Autonomous_state_history{i}(:, 2))), 'r*')
    hold off
    
    %~~~~~~~~
    drawnow 
      % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if n == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
  end