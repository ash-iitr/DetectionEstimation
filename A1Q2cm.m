clc;
clear;
close all;

%%Parameters

N=100;
noise_variance = 1;

sigma_square_est = zeros(1,N);
for l=1:N
    x = sqrt(noise_variance)*randn(N,1);
    sigma_square_est(l) = (1/N) * sum(x.^2);
end

%%Histogram Plot
histogram(sigma_square_est)
hold on
[vals,edges] = histcounts(sigma_square_est);
edges = edges(2:end) - (edges(2)-edges(1))/2;
plot(edges,vals)
hold off
title('Histogram for variance estimator For N=100')
xlabel('Variance estimator')

%%Mean and variance based on the simulation
mean_estimate = mean(sigma_square_est);
variance_estimate = var(sigma_square_est);