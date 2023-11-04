function [x ,gx] = n_discretized(mu, sigma, grid_length, support)

% function: n_discretized
% created June 2023
% discretized normal distribution

% in:
% mu>0, sigma>0: parameters 
% grid_length = 1,2,...: length of the grid
% support \in (0,1): size of support of the discretized distribution, in
% terms of share of total mass of a normal cdf with param mu and sigma (centered).

% out:
% x: grid vector of random values
% gx: grid vector of prob. mass function 

%%

% log normal match quality

% bounds for temp. grid
x1 = norminv((1-support)/2, mu, sigma);
x2 = norminv(1-(1-support)/2, mu, sigma);

% temp grid for rand values
x_tmp = linspace(x1, x2, grid_length+1)';

% temp grid for cdf
cdf_tmp = [normcdf(x_tmp, mu, sigma)];

% initialize vectors for final grid vectors 
x = zeros(grid_length,1);
gx = zeros(grid_length,1);

% fill in
for i = 1:grid_length
    x(i) = 1/2 * (x_tmp(i) + x_tmp(i+1));
    gx(i) = cdf_tmp(i+1) - cdf_tmp(i);
end

% rescale to sumup to one
% gx = gx/support;

tmp = 1/2*(1 - sum(gx));
gx(1) = gx(1) + tmp;
gx(end) = gx(end) + tmp;

% check
if sum(gx<0)>0 || abs(1-sum(gx))>1e-8
disp('discretized distribution is not well defined')
end

end