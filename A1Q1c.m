clc;
clear;
close all;

% Parameters
A = -25:1:24;
N = 1000; % Number of observations
sigma = 0.1; % Noise variance

% White Gaussian Noise samples
noise = sqrt(sigma) * randn(N, 1);

% Initialization
x = zeros(N, 1);
MSE = zeros(length(A), 4);

for A_id = 1:length(A)
    for k = 1:N
        x(k) = A(A_id) + noise(k);
    end

    % Estimators
    A1 = x(1);
    A2 = mean(x);
    A3 = 0.5 * mean(x);
    A4 = (A(A_id)^2 / (A(A_id)^2 + (sigma/N))) * mean(x);

    Est = [A1, A2, A3, A4];

    % Calculate Mean Squared Error (MSE) for each estimator
    for l = 1:length(Est)
        MSE(A_id, l) = mean(abs(Est(l) - A(A_id))^2);
    end
end

% Plotting
subplot(2,2,1)
plot(A, MSE(:,1))
xlabel('DC Level A')
ylabel('MSE -->')
title('Estimator - {A1}')

subplot(2,2,2)
plot(A, MSE(:,2))
xlabel('DC Level  A')
ylabel('MSE -->')
title('Estimator - {A2}')

subplot(2,2,3)
plot(A, MSE(:,3))
xlabel('DC Level A')
ylabel('MSE -->')
title('Estimator - {A3}')

subplot(2,2,4)
plot(A, MSE(:,4))
xlabel('DC Level A')
ylabel('MSE -->')
title('Estimator - {A4}')
