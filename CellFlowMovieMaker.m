%%%% Written by DIPJYOTI DAS (dipjyoti.das@yale.edu,
%%%% dipjyoti.das85@gmail.com)


clear all; close all;

filename = 'Cell_Flow_noise-0.9.dat'; %% Give input file name here
history = dlmread(filename); 

%% Define parameters and variables

 T = history(:,1); % 1st column is the Time Step.
 X = history(:,2); % 2nd column is the x position of the particles.
 Y = history(:,3); % 3rd column is the y position of the particles.
 
 VX = history(:,4);    % 4th column is the x-component of particle velocity.
 VY = history(:,5);    % 5th column is the y-component of particle velocity.
 Index = history(:,6); % 6th column is the particle index.
 
ymax=max(Y);
ymin=min(Y);
tmin = min(T); % Time at beginning.
tmax = max(T); % Time at end
dt = 1;  % How often would you like to capture the frame.
 %% Define movie name etc. 
 
 figure(1), % Shall open an empty figure window. Adjust it to the size you want. 
 % Ideally, you should have a fixed figure with specified position and size, but 
 % I find that very constrained. So, I just open an empty figure window, and adjust 
 % it according to my need. Don't touch it until the code completes. It'll interrupt 
 % the program, and the movie won't be produced. I sometime use this method to kill a 
 % process and redo the movie to my like.   
 
% clear movObj
% movObj = VideoWriter('Particle_Evolution.avi'); % Give whatever name you want to give. Can be automated. 
%movObj.FrameRate = 10; % Change to whatever value you want. 
%movObj.Quality = 75; % Any number between 0-100. 
 

 
 %% The best part: Frame capture 
 
 for t = tmin:dt:tmax %290
 
	id = T == t; % All data index for time t. 
	
    
	x = X(id); 
	y = Y(id); 
	vx = VX(id); 
	vy = VY(id); 
	
    
    plot(x,y,'ob','LineWidth',1.2,'MarkerSize',6,'MarkerFaceColor','c');
    hold on;
	quiver(x,y,vx,vy,0,'k','LineWidth',1.4) % Plot the velocity vectors
	%axis square  % add other axis properties if you feel like. 
	xlabel('X'); 
	ylabel('Y'); 
	title(['Time = ' num2str(t)])
    axis equal
    xlim([-11.5 11.5]);
    %ylim([-1 20]);
    ylim([-.1 ymax+1]);
	hold off;
	% Movie stuff	
    M(t)=getframe;
	frame = getframe(gcf); % Capture the frame
	%writeVideo(movObj,frame); % Add the frame to the movie.
 end
 numtimes=1;
 fps=35;
 movie(M,numtimes,fps)
 %close(movObj); % Close the movie variable and write the file. 
%movie2avi(M, 'particle-history.avi'); %, 'compression', 'none');
