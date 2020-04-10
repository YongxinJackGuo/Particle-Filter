% This function takes particle weight and particle states as the input and
% resample the particles based upon the particle weight. 
% P_weight will be N X 1 with N equals to the number of particles
% X_particles will be N X dof with dof equals to degree of freedom
% The function will output the resampled particles using the specific
% resample method.
function [X_resampled] = getResample(P_weight, X_particles)
    [N, dof] = size(X_particles);
    % initialize the resample matrix
    X_resampled = X_particles;
    % initialize the effective sample size ESSt and coefficient of variation Cvt^2;
    ESS = 0;
    cvt2 = 0;
    % initialize the particle population depletion coefficient
    beta = 0.45;
    % compute cvt2
    for i = 1:N
        cvt2 = cvt2 + (1/N) * (N * P_weight(i) - 1)^ 2;
    end
    % compute ESS
    ESS = N / (1 + cvt2);
    
    if ESS < beta * N
        for i = 1: N
            X_resampled(i, :) = X_particles(find(rand <= cumsum(P_weight), 1), :);
            %X_resampled(i, :) = X_particles(find(rand <= mvncdf(P_weight), 1), :);
        end
    end
end

