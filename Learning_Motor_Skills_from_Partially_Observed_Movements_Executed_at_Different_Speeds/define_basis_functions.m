% define basis functions

function PSIs_matrix = define_basis_functions(z_T, N, alpha)

    % z_T: phase function z evaluated at the last time step
    % N: number of basis functions
        
    h = z_T;
    
    %% Define generic basis function
    psi = @(n,z_t) exp(-0.5*(z_t-(n-1)*z_T/(N-1))^2/h);
    
    %% Define PSIs_matrix
    PSIs_matrix = NaN(N,z_T);
    
    z = alpha:alpha:z_T;
    
    for n = 1:N
        for phase_index = 1:size(z,2)
            z_t = z(phase_index);
            PSIs_matrix(n,phase_index) = psi(n,z_t);                        
        end
    end

end