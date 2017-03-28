clear all;
close all;
clc;

% Written in Matlab 2017a

% Meghasyam Tummalacherla
% ECE 285: Special Topics in Linear Algebra
% Part of the project Orthogonal Matching Pursuit (OMP) and variants

%% Stage Wise Orthogonal Matching pursuit: Setting up the Problem

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
% 
% figure(1);
% subplot(2,1,1);
% stem(x);
% title('The sparse vector x');

% Generating the noise n_eta
sigma_eta = 1; % The variance of the noise signal
mean_eta = 0; % Mean of noise is zero (given)
n_eta = sigma_eta*randn([m 1]) + mean_eta;
n_eta = n_eta/norm(n_eta);% Making the noise unit norm
norm_n_eta = 0; % Setting the norm of the noise
n_eta = norm_n_eta*n_eta; 

%% The algorithm StOMP

y = phi_mat*x + n_eta; % The measured signal
n = size(phi_mat,1); % The size of the required signal x
max_stages = 10; % Max no.of stages for our StOMP algorithm

% Stopping threshold for error
if norm_n_eta > 0
    thresh = norm_n_eta;
else
thresh = 10^-10; 
end

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
    phi_mat_I_s = zeros(n, length(I_s)); 
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

% Plots to compare the actual signal and the estimated signal x_s
figure(2);
subplot(2,1,1);
stem(x);
title(['Actual signal, s = ' num2str(s0)]);
xlabel('index');
ylabel('magnitude');

figure(2);
subplot(2,1,2);
stem(x_s);
title(['Estimated signal, ' '||noise||_2 = ' num2str(norm_n_eta)]);
xlabel('index');
ylabel('magnitude');
