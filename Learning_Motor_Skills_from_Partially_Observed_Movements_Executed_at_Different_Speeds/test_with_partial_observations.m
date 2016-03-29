% test with partial observations

n_testTrajectories = 20;

beginningFileName = './Drawings/Shifting_to_the_same_start/As/As_with_Gaussian_alphas_and_Gaussian_weights/Test/A';

%% load test trajectories
testTrajectories = cell(n_testTrajectories,2); % First column of this cell has matrices with the x,y coordinates of the trajectories.
                                                       % Second column of
                                                       % this cell has the
                                                       % alpha of each
                                                       % trajectory.

for i = 1:n_testTrajectories
   fileName = [beginningFileName int2str(i) '.txt'];
   testTrajectories{i,1} = dlmread(fileName);
   testTrajectories{i,2} = findAlpha(testTrajectories{i,1});
end

%% compute log_p(alpha) for a number of alphas
n_alpha_samples = 50;

%% defining deterministically a number of alphas between the minimum and the maximum observed
sampled_alphas = linspace(min(alphas), max(alphas), n_alpha_samples)';

%%
log_p_alphas = log(mvnpdf(sampled_alphas,mu_alpha,Sigma_alpha));

log_p_alphas_observations = NaN(n_alpha_samples, n_testTrajectories);

%% compute the PSIs matrix for each sampled alpha
PSIs_cell_test = cell(n_alpha_samples,1);
for i = 1:n_alpha_samples
    PSIs_cell_test{i} = define_basis_functions(100, N, sampled_alphas(i));
end


figure(5);
set_fig_position([0.161 0.275 0.297 0.468]);
pause;
for test_index = 1:n_testTrajectories
    
    alpha_distances = abs(bsxfun(@minus, sampled_alphas, testTrajectories{test_index,2}));
    [~,index_of_closest_alpha] = min(alpha_distances);
    
    percentage = 0.5;
    for length_obs = percentage*size(testTrajectories{test_index,1},1)
    
        % define which part of the test trajectories is actually observed
        observation = testTrajectories{test_index,1}(1:floor(length_obs),:);
        
        figure(5);
        clf;
        subplot(1,2,1);
        hold on;
        axis([-1 1 -1 1]);
        title(['observed part = ' num2str(floor(length_obs)/size(testTrajectories{test_index,1},1))]);
        plot(observation(:,1), observation(:,2), '-b');
        
        %% compute log_p(observation|alpha) for each of the sampled alphas
        for i = 1:n_alpha_samples
            log_prob_of_observation = compute_log_probability_of_observation_given_alpha(observation, PSIs_cell_test{i}, mu_weights, Sigma_weights, sampled_alphas(i), inv_obs_noise);
            
            %% compute log_p(alpha|observation) \propto log_p(observation|alpha) + log_p(alpha) for each of the sampled alphas
            log_p_alphas_observations(i,test_index) = log_prob_of_observation + log_p_alphas(i);
        end
               
        subplot(1,2,2);
        hold on;
        ylim([0 1]);
        xlabel('index of alpha');
        ylabel('probability of alpha given observation');
        test_index
        
        % mark the sampled alpha that is the closest to the actual one
        plot([index_of_closest_alpha index_of_closest_alpha], [0 1], '-r');
        
        p_alphas_observation = exp(log_p_alphas_observations(:,test_index) - max(log_p_alphas_observations(:,test_index))) / sum(exp(log_p_alphas_observations(:,test_index) - max(log_p_alphas_observations(:,test_index))));
        plot(p_alphas_observation);
        
        [max_p_alphas_observation, index_max_p_alphas_observation] = max(p_alphas_observation);
        
        title(['max probability = ' num2str(max_p_alphas_observation)]);        
                       
    end
    
    %% complete observation using the best alpha
    subplot(1,2,1);
    
    best_alpha = sampled_alphas(index_max_p_alphas_observation);
    
    PSIs_matrix = define_basis_functions(100, N, best_alpha);
    [matrix_u, matrix_variance] = complete_observation_with_best_alpha(observation, PSIs_matrix, mu_weights, Sigma_weights, best_alpha, inv_obs_noise);
    
    original_timeStamps = linspace(0,1,size(best_alpha:best_alpha:100,2));
    reference_timeStamps = linspace(0,1,floor(100/best_alpha));
    matrix_u_on_time = interp1(original_timeStamps, matrix_u, reference_timeStamps);
   
    t = size(observation,1):size(matrix_u_on_time,1);
    plot(matrix_u_on_time(t,1),matrix_u_on_time(t,2), '-r');
    
    plot(testTrajectories{test_index,1}(floor(length_obs)+1:end,1), testTrajectories{test_index,1}(floor(length_obs)+1:end,2), '-k');
    
    legend('observation', 'prediction', 'ground truth');
    xlabel('X');
    ylabel('Y');
    
    figure(6);
    set_fig_position([0.47 0.276 0.292 0.466]);
    clf;
    subplot(1,2,1);
    hold on;
    plot(testTrajectories{test_index,1}(:,1), '-b');
    plot(observation(:,1), '*b');

    shadedErrorBar(1:size(matrix_u_on_time,1), matrix_u_on_time(:,1), 2*sqrt(matrix_variance(:,1)), '-r', 1);
    legend('ground truth', 'observations', 'prediction');
    xlabel('time');
    ylabel('x trajectory');
    
    subplot(1,2,2);
    hold on;
    plot(testTrajectories{test_index,1}(:,2), '-b');
    plot(observation(:,2), '*b');

    shadedErrorBar(1:size(matrix_u_on_time,1), matrix_u_on_time(:,2), 2*sqrt(matrix_variance(:,2)), '-r', 1);
    legend('ground truth', 'observations', 'prediction');
    xlabel('time');
    ylabel('y trajectory');
            
    pause;
end

