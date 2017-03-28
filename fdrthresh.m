% Written in Matlab 2017a

% This function was inspired by and heavily based on the original function
% fdrthresh(z,q) from the Stanford SparseLab Version:100, created by Iddo Drori in
% 2006. Email: sparselab@stanford.edu

% This version is written by Meghasyam Tummalacherla, when he as a graduate student at 
% UC San Diego, Email: meghasyam@ucsd.edu

function thresh = fdrthresh(z, q)

az = abs(z); % Calculate the absolute values of z
[sz,iz] = sort(az,'descend'); % Sort the absolute valued array
pobs = erfc(sz./sqrt(2)); % Find the p-values using erfc

% Create the pnull vector to cross check with the p-values
N = 1:length(z); 
pnull =  N ./length(z);

% converting both the vectors into column vectors so that the dimensions
% match during the vector inequality
pnull = pnull(:);
pobs = pobs(:);

% Set of indices that satisfy the vector inequality, i.e the fraction of
% values that are false positives
pset = ((pobs) <= (q .* pnull));

% Finding the threshold

if any(pset)% If the pset is non-empty
    izmax  = iz(max(N(pset))); % Find the index in the original absolute value array from the sorted array indices
    thresh = az(izmax); % Set threshold with the found index
    
else % If the pset is empty
    thresh = max(az) + .01; % Set threshold as 0.1 greater than the maximum absolute value of the original array
end

% Limiting the threshold value between 2 and 3
if thresh<2
    thresh = 2;
end

if thresh>3
    thresh = 3;
end

end
