clear all;
close all;
clc;

% Written in Matlab 2017a

% Meghasyam Tummalacherla
% ECE 285: Special Topics in Linear Algebra
% Part of the project Orthogonal Matching Pursuit (OMP) and variants

%% CoSaMP: Setting up the Problem

% Define the size of the random matrix phi = mxn
m = 50; % no. of rows
n = 90; % no. of columns

% The entries of phi are i.i.d entries draw from a standard normal
% distribution
phi_mat = randn(m,n);

% Making phi_mat a sample from the Uniform Spherical Ensemble (USE)
for col_no = 1:size(phi_mat,2)
    phi_mat(:, col_no) = phi_mat(:, col_no)/norm(phi_mat(:, col_no));
end

% Sparsity of the sparse vector
s0 = 15; % No of non zero elements in x
s0_actual = s0;

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


%% The algorithm CoSaMP
u = phi_mat*x + n_eta;
% figure(1);
% subplot(2,1,2);
% stem(u);
% title('The measured signal u');

% Setting up initial parameters for the CoSaMP:
gamma_curr = u; % gamma_0
x_curr = 0*x; % Initial estimate of x
iter_no = 1;% Initialize the iteration index
error_iter_main = Inf; % Initialize the error in the iteration to a large no (Infinity here), as it will be updated later

% Stopping threshold for error
if norm_n_eta > 0
    thresh = norm_n_eta;
else
thresh = 10^-10; 
end

T_col_set = []; % Initialize the set of columns that are a support for the signal proxy
sup_col_now = []; % Initialize the suport of the current estimate of x
max_iter = 50; % Maximum number of iterations
err_vec = []; % Intialize the vector of errors across iterations
sparsity_est = s0; % The sparsity of the required answer ( e.g.,  we can have a 12 sparse estimate of a 15 sparse vector) 
s0 = sparsity_est; % Setting the sparsity of the required answer to the one specified above

while(iter_no < max_iter) % Check for stopping criteria
    iter_no = iter_no + 1; % Update iteration
    
    y = phi_mat'*gamma_curr; % Form signal proxy
    
    [val col_now] = sort(abs(y), 'descend'); % Sort the column using absolute values
    col_now = col_now(1:2*s0); % Find the 2*s0 most prominent components of the proxy
    col_now = col_now(:); % Making the set a column vector for consistency in the calculations below
    
    T_col_set = unique([sup_col_now; col_now]); % Merging the supports (in indices of the columns)
    
    % Finding the basis set of the merged support
    T_set =[];
    for col=1:length(T_col_set)
        T_set = [T_set phi_mat(:, T_col_set(col))];
    end
    
    % Stopping criteria, to see if any new data is added, if not then exit
    % the loop
    if length(sup_col_now) == length(T_col_set)
        break
    end
    
    % Finding the estimate of x in reduced dimension using the basis of the
    % merged support
    x_reduced_dim = pinv(T_set)*u;
    
    % Finding the s0 most promininent components of the reduced
    % dimension estimate of x, which will lead to the s0 approximation of x
    [val pos] = sort(abs(x_reduced_dim), 'descend');
    x_curr = zeros(size(x));
    for pos_no = 1:s0
        x_curr(T_col_set(pos(pos_no))) = x_reduced_dim(pos(pos_no));
    end
    
    % Updating the column support (sparsity s0) of the current estimate of x
    [val sup_col_now] = sort(abs(x_curr), 'descend');
    sup_col_now = sup_col_now(1:s0);
    
    % Updating the residual
    gamma_curr = u-phi_mat*x_curr;
    
    % Updating the error vector
    err_vec = [ err_vec norm(x_curr - x)];
    
    % Checking for the stopping criteria
    if norm(gamma_curr) < thresh
        break
    end   
    
end

% Plots to compare the actual signal and the estimated signal x_curr
figure(2);
subplot(2,1,1);
stem(x);
title(['Actual signal, s = ' num2str(s0_actual)]);
xlabel('index');
ylabel('magnitude');

subplot(2,1,2);
stem(x_curr);
title(['Estimated signal, ' '||noise||_2 = ' num2str(norm_n_eta)]);
xlabel('index');
ylabel('magnitude');
