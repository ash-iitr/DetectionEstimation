clc;
clear;
close all;
A=1;
variance = 1;
N =[10,100];
f = 0.01:0.005:0.49;
csum = zeros(length(f),length(N));
for k = 1:length(N)
    for l = 1:length(f)
        for n = 1:N(k)
            csum(l,k) = csum(l,k) + (4*pi^2*A^2*n^2*(sin(2*pi*l*n))^2);
            
        end
        CRLB(l,k) = variance/csum(l,k);
    end
end

figure(1);
plot(f,CRLB(:,1));
title('CRLB for N = 10');
 xlabel('frequency')
 ylabel('CRLB Value')
figure(2);
plot(f,CRLB(:,2));
title('CRLB for N = 100');
 xlabel('frequency')
 ylabel('CRLB Value')

