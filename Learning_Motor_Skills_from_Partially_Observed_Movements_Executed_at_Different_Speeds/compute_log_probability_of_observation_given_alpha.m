% compute probability of observation given alpha

function log_prob_of_observation = compute_log_probability_of_observation_given_alpha(observation, PSIs_matrix, mu, covariance, alpha, inv_obs_noise)

    length_obs_t = size(observation,1);
    length_obs_z = size(alpha:alpha:alpha*length_obs_t, 2);
    
    if alpha*length_obs_t < 100
    % define the observation as a function of z
    original_timeStamps = linspace(0,1,length_obs_t);
    reference_timeStamps = linspace(0,1,length_obs_z);
    observation = interp1(original_timeStamps, observation, reference_timeStamps);
    
    observation = [observation; NaN(size(alpha:alpha:100, 2)-length_obs_z, size(observation,2))];
        
    observation = observation(:);
    
    indexOfObservedValues = find(~isnan(observation));
    
    A = blkdiag(PSIs_matrix', PSIs_matrix'); % the basis functions are defined on z
    obsA = A(indexOfObservedValues,:);
    
    a = inv_obs_noise;
    
    u = obsA*mu;
                    
    sigma = (1/a)*eye(size(obsA,1)) + obsA*covariance*obsA';
    
    log_prob_of_observation = my_log_mvnpdf(observation(indexOfObservedValues)', u', sigma);
    else
        log_prob_of_observation = -Inf;
    end
end