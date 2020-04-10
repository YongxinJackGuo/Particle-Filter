% This codes perform a particle filter to a 2D mobile robot to solve a
% localization problem given a preset environmental and operational noise
% following a given task.
clc
clear
close all

%%
N = 5000; % the number of particles;
dof = 2; % degree of freedom
P_var_x = 1; % Noise variance in process for acceleration along x 
P_var_y = 1; % Noise variance in process for acceleration along y
R_var_x = 5; % Noise variance in process for measurement at x-coordinate
R_var_y = 5; % Noise variance in process for measurement at y-coordinate
totalTime = 200; % total time elapsed
dt = 1; % constant time step
vx0 = 0; % initial x-axis speed;
vy0 = 0; % initial y-axis speed;
v0 = [vx0 vy0];
ax = @(t)2*sin(1.2 * pi * t^2) + 0.8; % non-linear acceleration along x;
ay = @(t)1.1 * cos(0.45 * pi * t^2) + 0.3; % non-linear acceleration along y;
ax1 = @(t)1.2 * cos(2 * pi * (t/totalTime));
ay1 = @(t)1.2 * sin(2 * pi * (t/totalTime));
a = @(t)[ax1(t) ay1(t)];
% put noise into matrix form
P_var = [P_var_x P_var_y];
R_var = [R_var_x R_var_y];

% define the initial actual state for x and y. which serves as the mean of the
% initial guassian distribution. we use initial actual state for simulated
% measurements. We also defined a "wrong" initial state.
x_initial = [0 0]; 
x_wrong_initial = [20 20];
% define the covariance for the initial state.
x_i_var = eye(2) * 200^2;

%% create particle array
guassian = 1;
uniform = 2;
current_initial = uniform;
switch current_initial
    case guassian
        x_particles = mvnrnd(x_wrong_initial, sqrt(x_i_var), N); % 2D initial guassian distribution
    case uniform
        x_particles = 400 * rand(N, 2) - 200; % 2D uniform distribution
    otherwise
        % do nothing
end
% plot initial guassian distribution
plot(x_particles(:,1), x_particles(:,2), 'ro');
xlabel('X-axis (m)');
ylabel('Y-axis (m)');
title(['Initial Particle Distribution at N = ', num2str(N)]);
legend('particles at');
% define a true trajectory array
x_true = [];
x_true_prev = x_initial;
% define a measurement array
%m_true = sqrt((x_initial * 20).^2 + 100) + sqrt(R_var) .* randn(1, 2);  % initial measurement
m_true = x_initial.^2/20 + sqrt(R_var) .* randn(1, 2);
% define the best estimate arrary
x_best_estimate = [];
%% 
% action along x: x(i+1) = x(i) + Vx(i)t + (1/2)axt^2;
% action along y: y(i+1) = y(i) + Vy(i)t + (1/2)ayt^2;
% put it in state-spaced format.
% x = Ax + Bu + noise
% A = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
% B = [0.5 * dt^2; 0.5 * dt^2; dt; dt];
count = 1;
figure
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
tic;
for i = 0: dt: totalTime
    x_true_update = x_true_prev + v0 * dt + 0.5 * dt^2 .* (a(i) + sqrt(P_var) .* randn(1, 2)); % update actual path
    x_true_prev = x_true_update; % store previous x and y position
    x_true = [x_true; x_true_update]; % store path into an array
    %m_true_update = sqrt((x_true_update * 20).^2 + 100) + sqrt(R_var) .* randn(1, 2); % update measurement of actual path with noise
    m_true_update = x_true_update.^2/20 + sqrt(R_var) .* randn(1, 2);
    m_true = [m_true; m_true_update]; % store measurement into an array
    for j = 1: N
        % apply action to the particles.
        x_particles(j, :) = x_particles(j, :) + v0 * dt + 0.5 * dt^2 .* (a(i) + sqrt(P_var) .* randn(1, 2));
        % compute the measurements of the particles without the noise added 
        %m_particles(j, :) = sqrt((x_particles(j, :) * 20).^2 + 100);
        m_particles(j, :) = x_particles(j, :).^2/20 + sqrt(R_var) .* randn(1, 2);
    end
    v0 = a(i) * dt; % update the speed
 
    % store the best estimate for every update
    x_best_estimate = [x_best_estimate; [mean(x_particles(:, 1)), mean(x_particles(:, 2))] ];
    hold on;
    grid on;
    plot(x_true(count, 1), x_true(count, 2), 'bo', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(x_true(:, 1), x_true(:, 2), '-g', 'LineWidth', 2);
    plot(x_best_estimate(:, 1), x_best_estimate(:, 2), '-r', 'LineWidth', 2);
    plot(x_particles(:, 1), x_particles(:, 2), 'ro');
    xlim([-60 80]);
    ylim([-20 120]);
    xlabel('X coordinate of the map (m)');
    ylabel('Y coordinate of the map (m)');
    legend('Current Actual Position', 'Historic Actual Path', 'Best PF Estimated Path', 'Particle Samples');
    title(['The Particle Filter Simulation when N = ', num2str(N)]);
    pause(0.05);
    count = count + 1;
    hold off;
    clf;
    
    % get particle weight
    x_particles_weight = getParticleWeight(m_particles, m_true_update, R_var);
    % use 2-norm to merge the weights of x and y into one single weight.
    x_particles_weight = sqrt(x_particles_weight(:,1).^2 + x_particles_weight(:,2).^2);
    % normalize the weight
    x_particles_weight = x_particles_weight./ sum(x_particles_weight);
    % start resampling
    x_particles = getResample(x_particles_weight, x_particles);
end

% plot section
set(gca,'fontsize',14);
figure
hold on
subplot(121);
plot([0: dt: totalTime], x_best_estimate(:, 1), 'r-', [0: dt: totalTime], x_true(:, 1), 'g-', 'LineWidth', 2);
legend('Best Estimated Path along X-axis from Particle Filter', 'Actual X-axis Path with Noise');
xlabel('Time (s)');
ylabel('X coordinate');
title('X-axis Performance of Particle Filter');
subplot(122);
plot([0: dt: totalTime], x_best_estimate(:, 2), 'r-', [0: dt: totalTime], x_true(:, 2), 'g-', 'LineWidth', 2);
legend('Best Estimated Path along Y-axis from Particle Filter', 'Actual Y-axis Path with Noise');
xlabel('Time (s)');
ylabel('Y coordinate');
title('Y-axis Performance of Particle Filter');

toc;


