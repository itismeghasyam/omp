clear all;
close all;
clc;

% Written in Matlab 2017a

% Meghasyam Tummalacherla
% ECE 285: Special Topics in Linear Algebra
% Orthogonal Matching Pursuit (OMP) and variants
% Code that plots the StOMP vs OMP performance plots from already generated data

%{

This code plots the already generated plots from the file StOMP_data.mat
Which are the plots of relative error (figure 6) of OMP vs StOMP, and 
Probability of error in support (PES) (figure 7) of OMP vs StOMP, for a 
given sparsity s0, with the matrix phi from the USE, the non zero elements 
of the sparse vector x are picked uniformly from the range [-10 -1] U [1 10].
There is noise added of norm ||n||, labeled in the figure

%}


load('StOMP_data.mat');

for s0 = 1:6
    figure(1);
    subplot(6,2,2*s0-1);
    imagesc(flipud(rel_error_cell_omp(:,:,s0)));
    title(['OMP, ' 's0 = ' num2str(s0_vec(s0)) ', ||\eta||_2 = ' num2str(n_eta_norm_vec(s0))]);
    colormap('gray');colorbar;caxis([0 1]);
    ax = gca;
    ax.YTickLabel = [50 40 30 20 10];
    
    subplot(6,2,2*s0);
    imagesc(flipud(rel_error_cell_stomp(:,:,s0)));
    title(['StOMP, ' 's0 = ' num2str(s0_vec(s0)) ', ||\eta||_2 = ' num2str(n_eta_norm_vec(s0))]);
    colormap('gray');colorbar;caxis([0 1]);
    ax = gca;
    ax.YTickLabel = [50 40 30 20 10];
    
    figure(2);
    subplot(6,2,2*s0-1);
    imagesc(flipud(pes_error_cell_omp(:,:,s0)));
    title(['OMP, ' 's0 = ' num2str(s0_vec(s0)) ', ||\eta||_2 = ' num2str(n_eta_norm_vec(s0))]);
    colormap('gray');colorbar;caxis([0 1]);
    ax = gca;
    ax.YTickLabel = [50 40 30 20 10];
    
    subplot(6,2,2*s0);
    imagesc(flipud(pes_error_cell_stomp(:,:,s0)));
    title(['StOMP, ' 's0 = ' num2str(s0_vec(s0)) ', ||\eta||_2 = ' num2str(n_eta_norm_vec(s0))]);
    colormap('gray');colorbar;caxis([0 1]);
    ax = gca;
    ax.YTickLabel = [50 40 30 20 10];
    
end
