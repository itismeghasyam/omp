clear all;
close all;
clc;

% Written in Matlab 2017a

% Meghasyam Tummalacherla
% ECE 285: Special Topics in Linear Algebra
% Part of the project Orthogonal Matching Pursuit (OMP) and variants

%% Orthogonal Matching pursuit: Setting up the Problem

% Define the size of the random matrix phi = mxn
m = 50; % no. of rows
n = 80; % no. of columns

% The entries of phi are i.i.d entries draw from a standard normal
% distribution
phi_mat = randn(m,n);

% Making phi_mat a sample from the Uniform Spherical Ensemble (USE)
for col_no = 1:size(phi_mat,2)
    phi_mat(:, col_no) = phi_mat(:, col_no)/norm(phi_mat(:, col_no));
end

% Sparsity of the sparse vector
s0 = 9; % No of non zero elements in x

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

% figure(1);
% subplot(2,1,1);
% stem(x);
% title('The sparse vector x');

% Generating the noise n_eta
sigma_eta = 1; % The variance of the noise signal
mean_eta = 0; % Mean of noise is zero (given)
n_eta = sigma_eta*randn([m 1]) + mean_eta;
n_eta = n_eta/norm(n_eta);% Making the noise unit norm
norm_n_eta = 0.5; % Setting the norm of the noise
n_eta = norm_n_eta*n_eta; 

%% The algorithm OMP
y = phi_mat*x + n_eta;
figure(1);
% subplot(2,1,2);
stem(y);
title(['The measured signal y, ' '||noise||_2 = ' num2str(norm_n_eta)]);
xlabel('index');
ylabel('magnitude');

% Setting up initial parameters for the OMP
gamma_curr = y; % gamma_0, initial residual
I_set = []; % set of columns that contribute to y
col_set = []; % index of the columns
max_iter = 100; % Maximum number of iterations
iter_no = 1; % initializing the iteration variable
error_iter = Inf; % initial error 

% Stopping threshold for error
if norm_n_eta > 0
    thresh = norm_n_eta;
else
thresh = 10^-10; 
end

norm_gamma = [norm(gamma_curr)]; % Initializing the norm of residual

while((iter_no < max_iter+1) && error_iter > thresh) % Checking for stopping criteria
%     iter_no
    max_ip = 0;
    
    % Finding the most contributing column
    for col_no = 1:size(phi_mat,2)
        curr_ip = abs(gamma_curr'*phi_mat(:,col_no));
        if curr_ip > max_ip
            max_ip = curr_ip;
            curr_col = col_no;
        end
    end
    col_set = [col_set curr_col]; % Merging basis indices
    I_set = [I_set phi_mat(:,curr_col)]; % Finding the merged basis

    % Finding the projection of the measured signal onto the merged basis
    z_curr = projection_meg(y, I_set);
    
    % Updating the residual
    gamma_curr = y - z_curr;
    
    % Checking for stopping criteria
    if norm(gamma_curr) > norm_gamma(end)
        break
    else
        norm_gamma = [norm_gamma norm(gamma_curr)];
    end
    
    % Update error and iteration index
    error_iter = norm(gamma_curr);
    iter_no = iter_no+1;
end

% Finding the final estimate of x
col_set = unique(col_set);
I_set = [];

% Finding x from I_set;
for col=1:length(col_set)
    I_set = [I_set phi_mat(:, col_set(col))];
end

% Solving the problem in reduced dimension of I_set
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

% Finding the reduced dimension version of x
x_reduced_dim = pinv(phi_reduced_dim)*y_reduced_dim;

% Filling up the sparse vector x with the elements from the reduced dimension
% of x
x_est = zeros(size(x));
for pos = 1:length(col_set)
    x_est(col_set(pos)) = x_reduced_dim(pos);
end

% Plots to compare the actual signal and the estimated signal x_est
figure(2);
subplot(2,1,1);
stem(x);
title(['Actual signal, s = ' num2str(s0)]);
xlabel('index');
ylabel('magnitude');

subplot(2,1,2);
stem(x_est);
title(['Estimated signal, ' '||noise||_2 = ' num2str(norm_n_eta)]);
xlabel('index');
ylabel('magnitude');


    

