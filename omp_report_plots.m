clear all;
close all;
clc;

% Written in Matlab 2017a

% Meghasyam Tummalacherla
% ECE 285: Special Topics in Linear Algebra
% Orthogonal Matching Pursuit (OMP) and variants
% Code that plots the OMP performance plots from already generated data

%{

This code plots the already generated plots from the file omp_data.mat
Which are the plots of relative error (figure 3), and Probability of error
in support (PES) (figure 4), for a given sparsity s0, with the matrix phi
from the USE, the non zero elements of the sparse vector x are picked
uniformly from the range [-10 -1] U [1 10]. There is no noise added


%}

load('omp_data.mat');
for s0 = 1:9
    figure(1);
    subplot(3,3,s0);
    imagesc(flipud(rel_error_cell(:,:,s0)));
    title(['s0 = ' num2str(s0)]);
    colormap('gray');colorbar;
    ax = gca;
    ax.YTickLabel = [50 40 30 20 10];
    
    figure(2);
    subplot(3,3,s0);
    imagesc(flipud(pes_error_cell(:,:,s0)));
    title(['s0 = ' num2str(s0)]);
    colormap('gray');colorbar;
    ax = gca;
    ax.YTickLabel = [50 40 30 20 10];
    
end
