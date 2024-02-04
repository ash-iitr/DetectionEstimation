clc;
clear;
close all;

N=1000;
noise_variance = 0.1; % sigma squared
A = 1;
w = noise_variance * (randn(N,1));

x = A + w;

estA(1) = x(1); % matlab indexing starts from 1, not from 0
estA(2) = (1/N) * sum(x);
estA(3) = 0.5 * estA(2);
estA(4) = (A^2 / (A^2 + noise_variance/N )) * estA(2);

subplot(411)
p = histogram(x, 'Normalization', 'probability'); hold on;
% plot(1:length(x),estA(1));
bin_values = p.BinEdges(1:end-1) + p.BinWidth/2;
bin_count = p.Values;
plot(bin_values,bin_count,'LineWidth',1.5); hold on;
histogram(estA(1),N)
ylabel('pdf')
xlabel('estimates - Estimator 1')

subplot(412)
p = histogram(x, 'Normalization', 'probability'); hold on;
% plot(1:length(x),estA(1));
bin_values = p.BinEdges(1:end-1) + p.BinWidth/2;
bin_count = p.Values;
plot(bin_values,bin_count,'LineWidth',1.5); hold on;
histogram(estA(2),N)
ylabel('pdf')
xlabel('estimates - Estimator 2')

subplot(413)
p = histogram(x, 'Normalization', 'probability'); hold on;
% plot(1:length(x),estA(1));
bin_values = p.BinEdges(1:end-1) + p.BinWidth/2;
bin_count = p.Values;
plot(bin_values,bin_count,'LineWidth',1.5); hold on;
histogram(estA(3),N)
ylabel('pdf')
xlabel('estimates - Estimator 3')

subplot(414)
p = histogram(x, 'Normalization', 'probability'); hold on;
% plot(1:length(x),estA(1));
bin_values = p.BinEdges(1:end-1) + p.BinWidth/2;
bin_count = p.Values;
plot(bin_values,bin_count,'LineWidth',1.5); hold on;
histogram(estA(4),N)
ylabel('pdf')
xlabel('estimates - Estimator 4')


