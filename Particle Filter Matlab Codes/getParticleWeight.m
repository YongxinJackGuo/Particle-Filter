% a function that takes particle measurement and current true measurement
% as input and output the weight matrix for current time step.
% the input format is N X dof for P_measure, N is the number of particles
% and 1 X dof for X_measure. X_R is noise measurement
function [weight] = getParticleWeight(P_measure, X_measure, X_R)
    [N, dof] = size(P_measure);
    weight = zeros(N, dof); % initialize the weight
    for i = 1: N
        weight(i, :) = (1 ./ sqrt(2 .* pi .* X_R)) .* exp(-(X_measure - P_measure(i, :)).^2 ./ (2 .* X_R)); % compute the weight based on the difference between particles' measurments and actual measurement.
    end
end

