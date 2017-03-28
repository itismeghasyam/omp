% Written in Matlab 2017a

% Code to find the projetion of a vector y onto the columns of the basis I

function y_proj = projection_meg(y, I);

% Define the projection matrix
projection_mat = I*pinv((I'*I))*I';

% Projection of y onto I
y_proj = projection_mat*y;

end