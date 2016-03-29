% compute probability of training trajectories with current parameters
% estimated by EM

function log_like = compute_log_probability_of_training_trajectories(n_trainingTrajectories, trainingTrajectories, PSIs_cell_training, mu_weights, Sigma_weights, inv_obs_noise)

    log_like = 0;
    
    for i = 1:n_trainingTrajectories
        trajectory = trainingTrajectories{i}(:);
        
        indexOfObservedValues = find(~isnan(trajectory));
        
        PSIs_matrix = PSIs_cell_training{i};
        
        A = blkdiag(PSIs_matrix', PSIs_matrix'); % the basis functions are defined on z
        obsA = A(indexOfObservedValues,:);
        
        a = inv_obs_noise;
        
        u = obsA*mu_weights;
        
        sigma = (1/a)*eye(size(obsA,1)) + obsA*Sigma_weights*obsA';
        
        log_like = log_like + my_log_mvnpdf(trajectory(indexOfObservedValues)', u', sigma);
    end
    
end