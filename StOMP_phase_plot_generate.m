clear all;
close all;
clc;

% Written in Matlab 2017a

%{

This code was used to generate the plots mentioned in the report for the
compared with StOMP vs OMP algorithm. Note: This code takes around 10hrs to run on a i7-6700HQ gig
with 16GB RAM, it is advised to check the cosamp_report_plots.m file to look at
pre calculated plots, or look at cosamp.m, if you need to run StOMP for a particular
instance of y, x and phi and noise(StOMP.m takes a maximum of a couple of seconds)

%}


% Meghasyam Tummalacherla
% ECE 285: Special Topics in Linear Algebra
% Orthogonal Matching Pursuit (OMP) and variants

%% Setting up the Problem
% s0_max = 9;
rel_error_cell_omp = zeros(50,80,6);
pes_error_cell_omp = zeros(50,80,6);

rel_error_cell_stomp = zeros(50,80,6);
pes_error_cell_stomp = zeros(50,80,6);

s0_vec = [4 4 6 6 8 8]; % Sparsity values
n_eta_norm_vec = [0.3 0.6 0.3 0.6 0.3 0.6]; % Corresponding noise norms
max_stages = 10; % Maximum no of stages for StOMP

for iter_index = 1: length(s0_vec)
    s0 = s0_vec(iter_index);
%     s0
    n_eta_norm = n_eta_norm_vec(iter_index);
    for monte_carlo = 1:1000
        for n=s0:80 % n<=80
            for m=1:min(n-1, 50) % m< n, and m <=50
                
                % The entries of phi are i.i.d entries draw from a standard normal
                % distribution
                phi_mat = randn(m,n);
                
                % Making phi_mat a sample from the Uniform Spherical Ensemble (USE)
                for col_no = 1:size(phi_mat,2)
                    phi_mat(:, col_no) = phi_mat(:, col_no)/norm(phi_mat(:, col_no));
                end
                
                % cardinality of the sparse vector
                % No of non zero elements in x
                
                % Range of x is from [-10 -1] U [1 10] and it is a uniform distribution
                % So we will construct a variable int_ind (interval indicator) which is a
                % uniform variable which can take values <0.5 for the negative interval and
                % >=0.5 for the positive interval.
                
                % Based on the int_ind, we will now pick up a uniform sample from the
                % selected interval
                
                x = zeros(n,1); % Initializing x
                
                % The intervals from which we pick up x
                int_1 = [-10 -1];
                int_2 = [10 1];
                
                % Deciding the random positions of x, which contain the non zero values
                x_pos = [];
                while length(x_pos)<s0
                    x_pos = unique(randi(n, [1 s0]));
                end
                
                % Filling up the non zero values from a uniform distribution
                for j=1:s0
                    
                    % Choosing the interval
                    int_ind = rand(1);
                    
                    if int_ind >=0.5
                        int_curr = int_1;
                    else
                        int_curr = int_2;
                    end
                    
                    % Drawing the uniform random value from the chosen interval
                    curr_val = int_curr(1)*rand(1) + int_curr(2);
                    
                    % Filling up the non zero positions
                    x(x_pos(j)) = curr_val;
                    
                end
                %
                % figure(1);
                % subplot(2,1,1);
                % stem(x);
                % title('The sparse vector x');
                
                % Generating the noise n_eta
                sigma_eta = 1; % The variance of the noise signal
                mean_eta = 0; % Mean of noise is zero (given)
                n_eta = sigma_eta*randn([m 1]) + mean_eta;
                n_eta = n_eta/norm(n_eta)*n_eta_norm;
                
                %% OMP
                % The measured signal y
                y = phi_mat*x + n_eta;
                
                gamma_curr = y; % gamma_0, initial residual
                I_set = []; % set of columns that contribute to y
                col_set = []; % index of the columns
                max_iter = 100; % Maximum number of iterations
                iter_no = 1; % initializing the iteration variable
                error_iter = Inf; % initial error
                thresh = n_eta_norm; % Stopping threshold for error
                norm_gamma = [norm(gamma_curr)]; % Initializing the norm of residual
                
                while((iter_no < max_iter+1) && error_iter > thresh)
                    max_ip = 0;
                    for col_no = 1:size(phi_mat,2)
                        curr_ip = abs(gamma_curr'*phi_mat(:,col_no));
                        if curr_ip > max_ip
                            max_ip = curr_ip;
                            curr_col = col_no;
                        end
                    end
                    col_set = [col_set curr_col];
                    I_set = [I_set phi_mat(:,curr_col)];
                    
                    z_curr = projection_meg(y, I_set);
                    
                    gamma_curr = y - z_curr;
                    
                    if norm(gamma_curr) > norm_gamma(end)
                        break
                    else
                        norm_gamma = [norm_gamma norm(gamma_curr)];
                    end
                    
                    error_iter = norm(gamma_curr);
                    iter_no = iter_no+1;
                end
                
                col_set = unique(col_set);
                
                I_set = [];
                
                % Finding x from I_set;
                for col=1:length(col_set)
                    I_set = [I_set phi_mat(:, col_set(col))];
                end
                
                y_reduced_dim =zeros(length(col_set),1);
                phi_reduced_dim = zeros(length(col_set),length(col_set));
                x_reduced_dim = zeros(length(col_set),1);
                
                for col = 1:length(col_set)
                    y_reduced_dim(col) = y'*I_set(:, col);
                end
                
                for row = 1:length(col_set)
                    for col = 1:length(col_set)
                        phi_reduced_dim(row, col) = I_set(:, row)'*I_set(:, col);
                    end
                end
                
                x_reduced_dim = pinv(phi_reduced_dim)*y_reduced_dim;
                x_est = zeros(size(x));
                for pos = 1:length(col_set)
                    x_est(col_set(pos)) = x_reduced_dim(pos);
                end
                
                
                rel_error_cell_omp(m,n,iter_index) = rel_error_cell_omp(m,n,iter_index)+norm(x_est - x)/norm(x);
                pes_error_cell_omp(m,n,iter_index) = pes_error_cell_omp(m,n,iter_index)+(max(length(x_pos), length(col_set))-length(intersect(x_pos, col_set)))/max(length(x_pos), length(col_set));
                
                %% StOMP
                y = phi_mat*x + n_eta;
                x_s = 0*x; % initial estimate of x
                residual_curr = y; % initial estimate of the residual
                I_s = []; % The support of the signal x_s, initialized to an empty set
                
                for stage_no = 1:max_stages
                    
                    % Estimate of the residual correlations for the current iteration
                    c_s = (phi_mat')*residual_curr;
                    
                    % The formal noise level
                    sigma_s = norm(residual_curr)/sqrt(n);
                    
                    % The allowable fraction of false positives
                    q = 0.5;
                    
                    % Finding the threshold parameter t_s using False Discovery rate (FDR)
                    % or False Dicovery control
                    t_s = fdrthresh(c_s/sigma_s, q);
                    
                    % Find support of c_s that is greater than the set threshold
                    J_s = find(abs(c_s) > sigma_s*t_s);
                    
                    I_snew = [I_s; J_s]; % The merged support
                    I_snew = unique(I_snew); % Disregarding repeated columns
                    
                    % Stopping criteria, if no new elements for the support are found, then
                    % halt
                    if length(I_snew) ==  length(I_s)
                        break
                    else
                        I_s = I_snew;
                    end
                    
                    % Updating the current I_s to the new I_snew
                    I_s = I_snew;
                    
                    % Finding the basis columns of the merged support
                    phi_mat_I_s = zeros(m, length(I_s));
                    for col_no = 1: length(I_s)
                        phi_mat_I_s(:, col_no) = phi_mat(:, I_s(col_no));
                    end
                    
                    % Estimating x_s_I_s from the current merged support
                    x_s_I_s = pinv((phi_mat_I_s')*phi_mat_I_s)*(phi_mat_I_s')*y;
                    
                    
                    % Filling up the values of x_s from x_s_I_s at the correct positions of
                    % the support
                    
                    x_s = zeros(size(x));
                    for col_no = 1:length(I_s)
                        x_s(I_s(col_no)) = x_s_I_s(col_no);
                    end
                    
                    % Updating the residual
                    residual_curr = y - phi_mat*x_s;
                    
                    % Stopping criteria, checking for threshold
                    if norm(residual_curr) < thresh
                        break
                    end
                    
                end
                
                rel_error_cell_stomp(m,n,iter_index) = rel_error_cell_stomp(m,n,iter_index)+norm(x_s - x)/norm(x);
                pes_error_cell_stomp(m,n,iter_index) = pes_error_cell_stomp(m,n,iter_index)+(max(length(x_pos), length(I_s)-length(intersect(x_pos, I_s))))/max(length(x_pos), length(I_s));
                
            end
            
        end
        
    end
    % Averaging over monte carlo sims
    rel_error_cell_omp(:,:,iter_index)=rel_error_cell_omp(:,:,iter_index)/1000;
    pes_error_cell_omp(:,:,iter_index)=pes_error_cell_omp(:,:,iter_index)/1000;
    rel_error_cell_stomp(:,:,iter_index)=rel_error_cell_stomp(:,:,iter_index)/1000;
    pes_error_cell_stomp(:,:,iter_index)=pes_error_cell_stomp(:,:,iter_index)/1000;
end

for s0 = 1:6
    figure(1);
    subplot(6,2,s0);
    imagesc(flipud(rel_error_cell_omp(:,:,s0)));
    title(['s0 = ' num2str(s0)]);
    colormap('gray');colorbar;
    caxis([0 1]);
    ax = gca;
    ax.YTickLabel = [50 40 30 20 10];
    
    subplot(6,2,s0+6);
    imagesc(flipud(rel_error_cell_stomp(:,:,s0)));
    title(['s0 = ' num2str(s0)]);
    colormap('gray');colorbar;
    caxis([0 1]);
    ax = gca;
    ax.YTickLabel = [50 40 30 20 10];
    
    figure(2);
    subplot(6,2,s0);
    imagesc(flipud(pes_error_cell_omp(:,:,s0)));
    title(['s0 = ' num2str(s0)]);
    colormap('gray');colorbar; 
    caxis([0 1]);
    ax = gca;
    ax.YTickLabel = [50 40 30 20 10];
    
    subplot(6,2,s0+6);
    imagesc(flipud(pes_error_cell_stomp(:,:,s0)));
    title(['s0 = ' num2str(s0)]);
    colormap('gray');colorbar; 
    caxis([0 1]);
    ax = gca;
    ax.YTickLabel = [50 40 30 20 10];
    
end


save('StOMP_data.mat');