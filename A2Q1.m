%Q1
% Parameters
a = 0.98;
sigma_u_squared = 1;
mu_s = 5;
sigma_s_squared = 1;

% Initial state estimate
x_hat_0 = mu_s;
P_0 = sigma_s_squared;

% Number of time steps
n_steps = 101;

% Initialize arrays to store results
x_hat = zeros(1, n_steps);
P = zeros(1, n_steps);
K = zeros(1, n_steps);
x_hat_minus = zeros(1, n_steps);
P_minus = zeros(1, n_steps);

% Kalman filter equations
for n = 1:n_steps
    if n == 1
        x_hat(n) = x_hat_0;
        P(n) = P_0;
    else
        % Prediction step
        x_hat_minus(n) = a * x_hat(n-1);
        P_minus(n) = a^2 * P(n-1) + sigma_u_squared;
        
        % Update step
        K(n) = P_minus(n) / (P_minus(n) + sigma_s_squared);
        x_hat(n) = x_hat_minus(n) + K(n) * (mu_s - x_hat_minus(n));
        P(n) = (1 - K(n)) * P_minus(n);
    end
    
    % Print results
    fprintf('Step %d:\n', n);
    fprintf('Predicted state: %.4f\n', x_hat_minus(n));
    fprintf('Predicted state covariance: %.4f\n', P_minus(n));
    fprintf('Kalman gain: %.4f\n', K(n));
    fprintf('State estimate: %.4f\n', x_hat(n));
    fprintf('State estimate variance: %.4f\n\n', P(n));
end



% Define parameters
a = 0.98;
sigma_u_sq = 1;
mu_s = 5;
sigma_s_sq = 1;

% Iterate from n = 0 to 100
n_max = 100;

% Part (a): Generate and plot s[n] and x[n]


figure(1);
for sigma_n_sq = [0.9, 1, 1.2]
    % Generate s[n] and x[n]
    s = zeros(1, n_max+1);
    x = zeros(1, n_max+1);
    s(1) = mu_s;
    x(1) = s(1) + sqrt(sigma_n_sq) * randn;
    for n = 1:n_max
        s(n+1) = a * s(n) + sqrt(sigma_u_sq) * randn;
        x(n+1) = s(n+1) + sqrt(sigma_n_sq) * randn;
    end

    % Plot s[n] and x[n]
    subplot(3, 1, find(sigma_n_sq == [0.9, 1, 1.2]));
    plot(0:n_max, s, 'b-', 0:n_max, x, 'r--');
    legend('s[n]', 'x[n]');
    title(['\sigma_n^2 = ', num2str(sigma_n_sq)]);
    xlabel('n');
    ylabel('Magnitude');
end

% Part (b): Compute and plot Kalman filter variables
figure(2);
for sigma_n_sq = [0.9, 1, 1.2]
    % Initialize Kalman filter variables
    x_hat = zeros(1, n_max+1);
    P = zeros(1, n_max+1);
    K = zeros(1, n_max+1);
    x_hat(1) = mu_s;
    P(1) = sigma_s_sq;

    for n = 1:n_max
        % Prediction step
        x_hat_pred = a * x_hat(n);
        P_pred = a^2 * P(n) + sigma_u_sq;

        % Correction step
        K(n) = P_pred / (P_pred + sigma_n_sq);
        x_hat(n+1) = x_hat_pred + K(n) * (x(n) - x_hat_pred);
        P(n+1) = (1 - K(n)) * P_pred;
    end

    % Plot Kalman filter variables
    subplot(3, 2, 2*find(sigma_n_sq == [0.9, 1, 1.2])-1);
    plot(0:n_max, x, 'r--', 0:n_max, x_hat(1:end), 'b-');
    legend('x[n]', 'x\_hat[n]');
    title(['\sigma_n^2 = ', num2str(sigma_n_sq), ' - Prediction']);

    subplot(3, 2, 2*find(sigma_n_sq == [0.9, 1, 1.2]));
    plot(1:n_max+1, P, 'b-', 1:n_max+1, K, 'r--');
    legend('P[n]', 'K[n]');
    title(['\sigma_n^2 = ', num2str(sigma_n_sq), ' - MSE and Gain']);
end

% Part (c): Compute and plot the innovation sequence, its autocorrelation, and PSD
figure(3);
for sigma_n_sq = [0.9, 1, 1.2]
    % Initialize Kalman filter variables
    x_hat = zeros(1, n_max+1);
    P = zeros(1, n_max+1);
    K = zeros(1, n_max+1);
    x_hat(1) = mu_s;
    P(1) = sigma_s_sq;

    % Compute the innovation sequence
    innovation = zeros(1, n_max+1);
    for n = 1:n_max
        % Prediction step
        x_hat_pred = a * x_hat(n);
        P_pred = a^2 * P(n) + sigma_u_sq;

        % Correction step
        K(n) = P_pred / (P_pred + sigma_n_sq);
        x_hat(n+1) = x_hat_pred + K(n) * (x(n) - x_hat_pred);
        P(n+1) = (1 - K(n)) * P_pred;

        % Compute the innovation
        innovation(n) = x(n) - x_hat_pred;
    end

    % Compute the autocorrelation of the innovation sequence
    r_innovation = xcorr(innovation, 'normalized');

    % Compute the power spectral density of the innovation sequence
    [psd_innovation, f] = pwelch(innovation, [], [], [], 1);

    % Plot the results
    subplot(3, 3, 3*find(sigma_n_sq == [0.9, 1, 1.2])-2);
    plot(0:n_max, innovation);
    title(['\sigma_n^2 = ', num2str(sigma_n_sq), ' - Innovation']);
    xlabel('n');
    ylabel('Magnitude');

    subplot(3, 3, 3*find(sigma_n_sq == [0.9, 1, 1.2])-1);
    plot(0:n_max*2, r_innovation);
    title(['\sigma_n^2 = ', num2str(sigma_n_sq), ' - Autocorrelation']);
    xlabel('Lag');
    ylabel('Magnitude');

    subplot(3, 3, 3*find(sigma_n_sq == [0.9, 1, 1.2]));
    plot(f, psd_innovation);
    title(['\sigma_n^2 = ', num2str(sigma_n_sq), ' - PSD']);
    xlabel('Frequency');
    ylabel('Magnitude');
end
