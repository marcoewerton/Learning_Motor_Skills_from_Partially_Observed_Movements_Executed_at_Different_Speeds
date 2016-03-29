% complete trajectory given observed part

function [matrix_u, matrix_variance] = complete_observation_with_best_alpha(observation, PSIs_matrix, mu, covariance, alpha, inv_obs_noise)

    length_obs_t = size(observation,1);
    length_obs_z = size( alpha:alpha:alpha*length_obs_t, 2);
    
    % define the observation as a function of z
    original_timeStamps = linspace(0,1,length_obs_t);
    reference_timeStamps = linspace(0,1,length_obs_z);
    observation = interp1(original_timeStamps, observation, reference_timeStamps);
    
    observation = [observation; NaN(size(alpha:alpha:100, 2)-length_obs_z, size(observation,2))];
        
    observation = observation(:);
    
    indexOfObservedValues = find(~isnan(observation));
    
    data = [indexOfObservedValues observation(indexOfObservedValues)];  % 1st column: index of observed values; 2nd column: observed values

    A = blkdiag(PSIs_matrix', PSIs_matrix');
    
    if max(indexOfObservedValues) > size(A,1)
       disp('here');
    end
    
    obsA = A(indexOfObservedValues,:); % I must take exactly the lines of A that correspond to observed values
    a = inv_obs_noise;
    
    covariance = covariance + 10^-7*eye(size(covariance,1));
    
    b = covariance^-1;
    
    lambda = a*(obsA'*obsA) + b;
    
    w_mu_posterior = lambda\(a*obsA'*data(:,2) + b*mu);
            
    u = A*w_mu_posterior;
                    
    variance = (1/a)*eye(size(A,1)) + A*(lambda\A');
    
    length_trajectory_on_phase = size(alpha:alpha:100, 2);
    
    matrix_u = NaN(length_trajectory_on_phase,2);
    for i = 1:length_trajectory_on_phase
        for j = 1:2
           matrix_u(i,j) = u( (j-1)*length_trajectory_on_phase+i ); 
        end       
    end
    
    matrix_variance = NaN(length_trajectory_on_phase,2);
    for i = 1:2
        matrix_variance(:,i) = diag(variance((i-1)*length_trajectory_on_phase+1:(i-1)*length_trajectory_on_phase+length_trajectory_on_phase, (i-1)*length_trajectory_on_phase+1:(i-1)*length_trajectory_on_phase+length_trajectory_on_phase));
    end
    
    
end