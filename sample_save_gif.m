



h = figure;
axis tight manual % this ensures that getframe() returns a consistent size
filename = '../plots/may_31/sample_sim.gif';
for n = 1:size(DT_history, 2)
    % Draw plot for y = x.^n
    
    visualize_trig_trig( DT_history{n}, density_history{n} );
    title(['time - ', num2str(n/100)]);
    xlim([0, 20])
    ylim([0, 20])
    caxis([0, 10])
    
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