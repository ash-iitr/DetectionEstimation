% Define the parameters
sigma_sq = 1; % Noise variance
P_FA = 1e-8; % Probability of false alarm

% Part (a): f_0 = 0.2, N = 50
f_0 = 0.2;
N = 50;
gamma = chi2inv(1 - P_FA, 2*N); % Decision threshold
A_range = linspace(0, 5, 100); % Range of signal amplitudes
P_D_a = 1 - ncx2cdf(gamma, 2*N, 2*N*A_range.^2/(2*sigma_sq), 'upper'); % Probability of detection

% Part (b): f_0 = 0.25, N = 25
f_0 = 0.25;
N = 25;
gamma = chi2inv(1 - P_FA, 2*N);
A_range = linspace(0, 5, 100);
P_D_b = 1 - ncx2cdf(gamma, 2*N, 2*N*A_range.^2/(2*sigma_sq), 'upper');

% Part (c): f_0 = 0.4, N = 10
f_0 = 0.4;
N = 10;
gamma = chi2inv(1 - P_FA, 2*N);
A_range = linspace(0, 5, 100);
P_D_c = 1 - ncx2cdf(gamma, 2*N, 2*N*A_range.^2/(2*sigma_sq), 'upper');

% Part (d): f_0 = 0.5, N = 30
f_0 = 0.5;
N = 30;
gamma = chi2inv(1 - P_FA, 2*N);
A_range = linspace(0, 5, 100);
P_D_d = 1 - ncx2cdf(gamma, 2*N, 2*N*A_range.^2/(2*sigma_sq), 'upper');

% Plot the results
figure;
subplot(2, 2, 1); plot(A_range, P_D_a); title('f_0 = 0.2, N = 50');
xlabel('Signal Amplitude, A'); ylabel('Probability of Detection, P_D');
subplot(2, 2, 2); plot(A_range, P_D_b); title('f_0 = 0.25, N = 25');
xlabel('Signal Amplitude, A'); ylabel('Probability of Detection, P_D');
subplot(2, 2, 3); plot(A_range, P_D_c); title('f_0 = 0.4, N = 10');
xlabel('Signal Amplitude, A'); ylabel('Probability of Detection, P_D');
subplot(2, 2, 4); plot(A_range, P_D_d); title('f_0 = 0.5, N = 30');
xlabel('Signal Amplitude, A'); ylabel('Probability of Detection, P_D');
