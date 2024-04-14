% Define the parameters
a = 0.9;
sigma_u_sq = 1;
n_max = 100;

% Initialize M[n|n]
M = zeros(1, n_max+1);
M(1) = 1; % M[-1|-1] = 1

for n = 1:n_max
    % Compute the Kalman gain
    sigma_n_sq = n + 1;
    K = (a^2 * M(n) + sigma_u_sq) / ((a^2 * M(n) + sigma_u_sq) + sigma_n_sq);

    % Compute M[n|n]
    M(n+1) = (1 - K) * (a^2 * M(n) + sigma_u_sq);
end

% Plot M[n|n]
figure;
plot(0:n_max, M);
title('M[n|n] for n');
xlabel('n');
ylabel('M[n|n]');