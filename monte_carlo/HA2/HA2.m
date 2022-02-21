%% Theoretical values
c_theory = [4,12,36,100,284,780,2172,5916,16268, 44100];

%% estimating c_n(2) for n = 1,2,3,4,5,6,7...
% Naive version

N = 1000;
x0 = [0,0];
n_lim = 10;
c2 = zeros(n_lim,1);
for current_n=1:n_lim
    sum = 0;
    for i=1:N
        sum = sum + naive(x0,current_n)*4^current_n;
    end
    c2(current_n) = sum/N;
end

c2
c_theory'
fprintf("\nEstimate of c_nu(2)")
figure
hold on
plot(c2)
title("Estimate of c_n(2) for N=" +N+ " (N_{SA}/N)")
ylabel("N_{SA}/N")
xlabel("Length of SAW [n]")
hold off

%% Free neighbours only, naive

N = 1000;
x0 = [0,0];
n_lim = 20;
c2 = zeros(n_lim,1);
c2s_SIS = zeros(n_lim,10);

for j=1:10
    for current_n=1:n_lim
        omegas = zeros(N,1);
        for i=1:N
            [cc, ww] = free_neighbours(x0,current_n);
            omegas(i) = cc*ww(end);
        end
        c2(current_n) = mean(omegas);
    end
    c2s_SIS(:,j) = c2;
end

c_mean_SIS = mean(c2s_SIS');

figure
hold on
plot(c2)
title("Estimate of c_n(2) for N=" +N+ " (N_{SA}/N)")
ylabel("N_{SA}/N")
xlabel("Length of SAW [n]")
hold off

%% SISR
clc
close all

N = 1000;
x0 = [0,0];
n_lim = 11;
c2 = zeros(n_lim,1);
c2s = zeros(n_lim-1,10);
for i=1:10
    %Initialization
    omegas = ones(N,2);
    c2(1) = 1;
    X = zeros(1,2,N);
    
    for current_n=1:n_lim-1
        [X,omegas] = sis(X, omegas);
        c2(current_n+1) = mean(omegas(:,2));
        indexes = randsample(N,N,true,omegas(:,2)/sum(omegas(:,2)));
        X = X(:,:,indexes);
        omegas(:,1) = omegas(:,2);
    end

    c2 = c2(2:end);
    c2s(:,i) = c2;
end

c_mean_SISR = mean(c2s');

figure
hold on
plot(c2)
title("Estimate of c_n(2) for N=" +N+ " (N_{SA}/N)")
ylabel("N_{SA}/N")
xlabel("Length of SAW [n]")
hold off


%SISR_diff = abs(c_theory - c_mean_SISR);
%SIS_diff = abs(c_theory - c_mean_SIS);
%nbr_n = [1,2,3,4,5,6,7,8,9];
%T = table(nbr_n, c_theory,c_mean_SIS,c_mean_SISR)
c_theory,c_mean_SIS,c_mean_SISR%,SIS_diff,SISR_diff

%% Question 6
% Estimate A_2, mu_2, gamma_2

% Relation
%           { A_d*mu_d^n*n^(gamma_d-1), d=1,2,3, d>=5
% c_n(d) ~  { 
A = zeros(length(c2s),1);
mu = zeros(length(c2s),1);
gamma = zeros(length(c2s),1);
for i=1:length(c2s)
    c_mean = log(c2s(:,i));
    x = (1:length(c_mean))';
    alpha = c_mean + log(x);
    X = [ones(length(c_mean),1), x, log(x)];
    beta = pinv(X)*alpha;
    A(i) = exp(beta(1));
    mu(i) = exp(beta(2));
    gamma(i) = beta(3);
end

mean_A = mean(A);
mean_mu = mean(mu);
mean_gamma = mean(gamma);
var_A = var(A);
var_mu = var(mu);
var_gamma = var(gamma);

fprintf("\nA_2 = %f\nmu_2 = %f\ngamma_2 = %f\n",mean_A,mean_mu,mean_gamma)
fprintf("\nVariance of A_2 = %f\nVariance of mu_2 = %f\nVariance of gamma_2 = %f\n",var_A,var_mu,var_gamma)


% b = x'\c_mean';
% figure
% hold on
% plot(c_mean)
% plot(b*x)
% hold off
% legend("log of c_2(2)","linear regression")


%% SISR multidimensional
clc
close all
tic
d = 3;
N = 1000;
n_lim = 10;
c2 = zeros(n_lim,1);
c2s = zeros(n_lim,10);
for i=1:10
    %Initialization
    omegas = ones(N,2);
    c2(1) = 1;
    X = zeros(1,d,N);
    
    for current_n=1:n_lim
        [X,omegas] = sis(X, omegas);
        
        indexes = randsample(N,N,true,omegas(:,2)/sum(omegas(:,2)));
        X = X(:,:,indexes);
        c2(current_n+1) = mean(omegas(:,2));
        omegas(:,1) = omegas(:,2);
        
    end

    c2 = c2(2:end);
    c2s(:,i) = c2;
end
toc

c_mean_SISR = mean(c2s');

figure
hold on
plot(c2)
title("Estimate of c_n(2) for N=" +N+ " (N_{SA}/N)")
ylabel("N_{SA}/N")
xlabel("Length of SAW [n]")
hold off

%% Estimating A, mu, gamma

A = zeros(length(c2s),1);
mu = zeros(length(c2s),1);
gamma = zeros(length(c2s),1);
for i=1:length(c2s)
    c_mean = log(c2s(:,i));
    x = (1:length(c_mean))';
    alpha = c_mean + log(x);
    X = [ones(length(c_mean),1), x, log(x)];
    beta = pinv(X)*alpha;
    A(i) = exp(beta(1));
    mu(i) = exp(beta(2));
    gamma(i) = beta(3);
end

mean_A = mean(A);
mean_mu = mean(mu);
mean_gamma = mean(gamma);
var_A = var(A);
var_mu = var(mu);
var_gamma = var(gamma);

fprintf("\nA_%d = %f\nmu_%d = %f\ngamma_%d = %f\n",d,mean_A,d,mean_mu,d,mean_gamma)
fprintf("\nVariance of A_%d = %f\nVariance of mu_%d = %f\nVariance of gamma_%d = %f\n",d,var_A,d,var_mu,d,var_gamma)

%% Estimating population measurement
% Estimate tau_k = E[ X_k | Y_{0:k} ], k = 1,2,3,...,50
close all
clc

load population.mat
% X is correct values for comparison
% Y is simulation with n=50 of the model, measurements 

N = 100000;
n = 50;

% Initialization
tau = zeros(n+1,1);
part = unifrnd(1/5, 3/5,N,1);

% Observation density Y_k|X_k
p = @(x,y) unifpdf(y, x/2, x);
w = p(part,Y(1));
tau(1) = sum(part.*w)/sum(w);
ind = randsample(N,N,true,w);
part = part(ind);

%intervals = zeros(n+1, 5);  % Steps are rows, tau, taulower, tauupper,
intervals = zeros(n+1, 4);  % conf.interval width, boolean if real X
                            % is inside conf. interval
[xx, I] = sort(part);
cw = cumsum(w(I)/sum(w));
Ilower=find(cw<= 0.025,1);
Iupper=find(cw>= 0.975,1);
taulower = xx(Ilower);
tauupper = xx(Iupper);
%intervals(1,:) = [tau(1), taulower, tauupper, tauupper-taulower, (X(1) >= taulower & X(1) < tauupper)];
intervals(1,:) = [tau(1), taulower, tauupper, (X(1) >= taulower & X(1) < tauupper)];
for i=1:n
    part = unifrnd(1/2,3,N,1).*part.*(1 - part);
    w = p(part,Y(i+1));
    tau(i+1) = sum(part.*w)/sum(w);
    
    [xx, I] = sort(part);
    cw = cumsum(w(I)/sum(w));
    Ilower=find(cw>= 0.025,1);
    Iupper=find(cw>= 0.975,1);
    taulower = xx(Ilower);
    tauupper = xx(Iupper);
    %intervals(i+1,:) = [tau(i+1), taulower, tauupper, tauupper-taulower...
    %    (X(i+1) >= taulower & X(i+1) < tauupper)];
    intervals(i+1,:) = [tau(i+1), taulower, tauupper,...
        (X(i+1) >= taulower & X(i+1) < tauupper)];
    ind = randsample(N,N,true,w);
    part = part(ind);
end
figure
hold on
plot(0:n, tau)
plot(0:n,X)
legend("\tau_n", "True values")
title("Estimate of relative population size for generation k, real values as reference")
xlabel("Generations [k]")
ylabel("Relative population size")
hold off

intervals
fprintf("\nPart of real X_n that are inside conf. interval of tau_n: %f\n",...
sum(intervals(:,4))/(n+1));





