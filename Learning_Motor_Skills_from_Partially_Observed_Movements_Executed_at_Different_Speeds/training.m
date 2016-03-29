% script to train phase estimator using EM

clear variables;
close;
clc;

%%%%%%%%%%%%%%%%%%%%%%%% user defined variables %%%%%%%%%%%%%%%%%%%%%%%%%%%

n_trainingTrajectories = 20;

% number of basis functions
N = 20;

% number of DoFs
n_DoFs = 2;

inv_obs_noise = 10000 ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load training trajectories
load('incompleteTrajectories_missing180steps');

PSIs_cell_training = cell(n_trainingTrajectories,1);
alphas = NaN(n_trainingTrajectories,1);

for i = 1:n_trainingTrajectories

   alphas(i) = findAlpha(trainingTrajectories{i,1});
   
   %% define basis functions for each demonstrated alpha
   PSIs_matrix = define_basis_functions(100, N, alphas(i));
   PSIs_cell_training{i} = PSIs_matrix;
   
    
end

%% estimate the parameters of p(alpha), assuming this is a Gaussian
mu_alpha = mean(alphas);
Sigma_alpha = cov(alphas);
std_alpha = std(alphas);

%% initialize mu_weights and Sigma_weights
mu_weights = 5*ones(n_DoFs*N, 1);
Sigma_weights = eye(n_DoFs*N, n_DoFs*N);

%% estimate the parameters of p(w), assuming this is a Gaussian, using Expectation-Maximization (EM)
mu_weights_for_each_trajectory = cell(n_trainingTrajectories,1);
Sigma_weights_for_each_trajectory = cell(n_trainingTrajectories,1);

n_em_iterations = 100;

loglikes = NaN(n_em_iterations,1);

for i = 1:n_em_iterations
    i
    sum_mus = zeros(size(mu_weights,1), size(mu_weights,2));
    sum_Sigmas = zeros(size(Sigma_weights,1), size(Sigma_weights,2));
    for j = 1:n_trainingTrajectories

        PSIs_matrix = PSIs_cell_training{j};
        
        B = blkdiag(PSIs_matrix', PSIs_matrix');
        
        alpha = alphas(j);
        z = alpha:alpha:100;

        trajectory = trainingTrajectories{j,1}(:);
        obs_trajectory = trajectory(~isnan(trajectory));
        
        obsB = B(~isnan(trajectory),:);
        
        Sigma_weights_for_each_trajectory{j} = (obsB'*inv_obs_noise*obsB + Sigma_weights^-1)^-1;
        mu_weights_for_each_trajectory{j} = Sigma_weights_for_each_trajectory{j}*(obsB'*inv_obs_noise*obs_trajectory + Sigma_weights^-1*mu_weights);
    
        sum_mus = sum_mus + mu_weights_for_each_trajectory{j};
        sum_Sigmas = sum_Sigmas + Sigma_weights_for_each_trajectory{j};
    end
    
    mu_weights = sum_mus/n_trainingTrajectories;
    E = NaN(n_trainingTrajectories,N*2);
    for k = 1:n_trainingTrajectories
        E(k,:) = (mu_weights_for_each_trajectory{k} - mu_weights)';
    end
    
    Sigma_weights = (E'*E + sum_Sigmas)/n_trainingTrajectories;
    Sigma_weights = Sigma_weights + 10^-5*eye(size(Sigma_weights,1), size(Sigma_weights,2));
    
    % Compute probability of all training trajectories with the current
    % parameters estimated by EM
    loglikes(i) = compute_log_probability_of_training_trajectories(n_trainingTrajectories, trainingTrajectories, PSIs_cell_training, mu_weights, Sigma_weights, inv_obs_noise);
    
end


% see if the weights and basis functions approximate the real trajectories
% well
close;
figure(3);
for i = 1:n_trainingTrajectories
    clf;
    hold on;
    weights = [mu_weights_for_each_trajectory{i}(1:N,1)'; mu_weights_for_each_trajectory{i}(N+1:end,1)'];
    approxTrajectory = PSIs_cell_training{i}'*weights';
    plot(trainingTrajectories{i,1}(:,1), trainingTrajectories{i,1}(:,2), '-*b');
    plot(approxTrajectory(:,1), approxTrajectory(:,2), '-r'); 
    pause;
end