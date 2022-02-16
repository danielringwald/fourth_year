%% Load
load powercurve_D236.mat
%% 2a Crude Monte Carlo
%close all
% Generate Weibull distributed random numbers: wblrn
% wblpdf
% wblcdf

lambdas = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7]; %Scale
ks = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];             %Shape
len = length(lambdas);
N = 10000;

winds = zeros(12, N);

for j=1:N
    winds(:, j) = wblrnd(lambdas, ks);
end

Ps = zeros(1, len);
powers = zeros(len, N);
for i=1:len
    P_month = zeros(1,N);
    for j=1:N
        power = P(winds(i,j));
        P_month(j) = power;
    end
    Ps(i) = mean(P_month);
    powers(i,:) = P_month;
end

sigmas = zeros(1,len);
for i=1:len
    sigmas(i) = std(powers(i,:));
end

intervals = zeros(3,len);
intervals(1,:) = Ps;
for i=1:len
    intervals(2,i) = -1.96*sigmas(i)/sqrt(N);
    intervals(3,i) = 1.96*sigmas(i)/sqrt(N);
end
intervals = intervals /1e6;
figure
hold on
plot(intervals(1,:))
plot(intervals(1,:)+intervals(2,:),'--r')
plot(intervals(1,:)+intervals(3,:),'--r')
hold off
xticks(1:12)
xticklabels({'Jan','Feb','March','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
ylabel("Power output [MW]")
legend("\tau","95% confidence interval")
title("STANDARD MONTE CARLO: Power output from wind turbine per month", "N = "+N+", 95% confidence interval")
sigmas
fprintf("\nUP: STD DOWN: VARIANCE\n")
sigmas.^2
I_upper = intervals(1,:) + intervals(3,:)
I_lower = intervals(1,:) + intervals(2,:)
I_upper - I_lower
%% 2a Truncated version
a = 3;
b = 30;
lambdas = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7]; %Scale
ks = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];             %Shape
month = length(lambdas);
N = 10000;

Xs = zeros(month,N);
Ps = zeros(1, month);
powers = zeros(month, N);
sigmas = zeros(1,month);

for i=1:month
    u = (wblcdf(b,lambdas(i),ks(i)) - wblcdf(a,lambdas(i),ks(i)))*rand(1,N) + wblcdf(a,lambdas(i),ks(i));
    Xs(i,:) = wblinv(u,lambdas(i),ks(i));
    powers(i,:) = P(Xs(i,:))*(wblcdf(b,lambdas(i),ks(i)) - wblcdf(a,lambdas(i),ks(i)));
    Ps(i) = mean(powers(i,:));
    sigmas(i) = std(powers(i,:));
end

intervals = zeros(3,month);
intervals(1,:) = Ps;
for i=1:month
    intervals(2,i) = -1.96*sigmas(i)/sqrt(N);
    intervals(3,i) = 1.96*sigmas(i)/sqrt(N);
end
intervals = intervals /1e6;

figure
hold on
plot(intervals(1,:))
plot(intervals(1,:)+intervals(2,:),'--r')
plot(intervals(1,:)+intervals(3,:),'--r')
hold off
xticks(1:12)
xticklabels({'Jan','Feb','March','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
ylabel("Power output [MW]")
legend("\tau","95% confidence interval")
title("TRUNCATED MONTE CARLO: Power output from wind turbine per month", "N = "+N+", 95% confidence interval")
sigmas
fprintf("\nUP: STD DOWN: VARIANCE\n")
sigmas.^2
I_upper = intervals(1,:) + intervals(3,:)
I_lower = intervals(1,:) + intervals(2,:)
I_upper - I_lower

%% 2b Control variate

lambdas = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7]; %Scale
ks = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];             %Shape

%Exact variance and expected value for the Weibull distribution
[M, V] = wblstat(lambdas, ks);

%Define length and amount of samples
len = length(lambdas);
N = 10000;

%Allocate memory for wind 
winds = zeros(12, N);

%Generate random Weibull distributed winds
for i=1:len
    winds(i, :) = wblrnd(lambdas(i), ks(i), [1,N]);
end

%Create matrix with all the power outputs per month, this is the objective
%function
powers = zeros(len, N);
for i=1:len
    P_month = zeros(1,N);
    for j=1:N
        power = P(winds(i,j));
        P_month(j) = power;
    end
    powers(i,:) = P_month;
end

%Add the fix with the control variate
Z = zeros(len,N);
for i=1:len
    beta = -cov(powers(i,:), winds(i,:))./ V(i);
    beta = beta(1,2);
    Z(i,:) = powers(i,:) + beta*(winds(i,:) - M(i));
end

%Find tau from the new Z
Ps = zeros(1, len);
for i=1:12
    Ps(i) = mean(Z(i,:));
end

%Find standard deviation of 
sigmasb = zeros(1,len);
for i=1:len
    sigmasb(i) = std(Z(i,:));
end

intervals = zeros(3,len);
intervals(1,:) = Ps;
for i=1:len
    intervals(2,i) = -1.96*sigmasb(i)/sqrt(N);
    intervals(3,i) = 1.96*sigmasb(i)/sqrt(N);
end

intervals = intervals /1e6;         %Make it into MW from W
figure
hold on
plot(intervals(1,:))
plot(intervals(1,:)+intervals(2,:),'--r')
plot(intervals(1,:)+intervals(3,:),'--r')
hold off
xticks(1:12)
xticklabels({'Jan','Feb','March','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
ylabel("Power output [MW]")
legend("\tau","95% confidence interval")
title("CONTROL VARIATE: Power output from wind turbine per month", "N = "+N+", 95% confidence interval")
sigmasb
fprintf("\nUP: STD DOWN: VARIANCE\n")
sigmasb.^2

% 2c????

% lambdas = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7]; %Scale
% ks = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];             %Shape
% len = length(lambdas);
% N = 5000;
% 
% winds = zeros(12, N);
% w = zeros(12,N);
% for j=1:N
%     winds(:, j) = wblrnd(lambdas, ks);
%     w(:, j) = wblrnd(2, ks);
% end
% 
% Ps = zeros(1, len);
% powers = zeros(len, N);
% for i=1:len
%     P_month = zeros(1,N);
%     for j=1:N
%         
%         power = P(winds(i,j)).*w(i,j);
%         P_month(j) = power;
%     end
%     Ps(i) = mean(P_month);
%     powers(i,:) = P_month;
% end
% 
% sigmas = zeros(10,len);
% for i=1:len
%     sigmas(i) = std(powers(i,:));
% end
% 
% 
% intervals = zeros(3,len);
% intervals(1,:) = Ps;
% for i=1:len
%     intervals(2,i) = -1.96*sigmas(i)/sqrt(N);
%     intervals(3,i) = 1.96*sigmas(i)/sqrt(N);
% end
% intervals = intervals /1e6;
% figure
% hold on
% plot(intervals(1,:))
% plot(intervals(1,:)+intervals(2,:),'--r')
% plot(intervals(1,:)+intervals(3,:),'--r')
% hold off
% xticks(1:12)
% xticklabels({'jan','feb','mars','apr','may','jun','jul','aug','sep','oct','nov','dec'})
% ylabel("Power output [MW]")
% dist1 = intervals(3,:)*2;
% mean4 = mean(dist1);
% mean4

%% 2c Importance sampling

lambdas = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7]; %Scale
ks = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];             %Shape
month = length(lambdas);
N = 10000;

f = @(x, month) wblpdf(x, lambdas(month), ks(month));
g = @(x, month) normpdf(x, 10, 7);

Frand = @(month) wblrnd(lambdas(month), ks(month), 1, N);
%Nrand = @(month) normrnd(10, 7, 1, N);

draw = zeros(month,N);
phiomega = zeros(month,N);
y = zeros(month,N);
tau = zeros(12,1);
std1 = zeros(12,1);
for i=1:month
   draw(i,:) = normrnd(13, 5, 1, N); 
   phiomega(i,:) = P(draw(i,:))'.*(f(draw(i,:),i)./g(draw(i,:),i));
   tau(i) = mean(phiomega(i,:));
   std1(i) = std(phiomega(i,:));
end
intervals = zeros(3,month);
intervals(1,:) = tau;
for i=1:month
    intervals(2,i) = -1.96*std1(i)/sqrt(N);
    intervals(3,i) = 1.96*std1(i)/sqrt(N);
end
intervals = intervals /1e6;
figure
hold on
plot(intervals(1,:))
plot(intervals(1,:)+intervals(2,:),'--r')
plot(intervals(1,:)+intervals(3,:),'--r')
hold off

xticks(1:12)
xticklabels({'Jan','Feb','March','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
ylabel("Power output [MW]")
legend("\tau","95% confidence interval")
title("IMPORTANCE SAMPLING: Power output from wind turbine per month", "N = "+N+", 95% confidence interval")
std1'
fprintf("\nUP: STD DOWN: VARIANCE\n")
std1.^2'


%% 2d antithetic sampling
a = 3;
b = 30;
lambdas = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7]; %Scale
ks = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];             %Shape
month = length(lambdas);
N = 10000;

%U =@(j) (wblcdf(b,lambdas(j),ks(j)) - wblcdf(a,lambdas(j),ks(j)))*rand(1,N/2) + wblcdf(a,lambdas(j),ks(j));
W = zeros(month,N/2);
tau = zeros(1,month);
std1 = zeros(1,12);
F = @(x,l,k) wblcdf(x, l,k);
Fx =@(x,l,k) (F(x,l,k) - F(a,l,k)) / (F(b,l,k) - F(a,l,k));
xhatu = 0;
xhatt = 0;

for j=1:month
%     for i=1:N/2
%         u = rand(1);
%         t = 1 - u;
%         
%         for l=3:30
%             if Fx(l,lambdas(j),ks(j)) >=u
%                 xhatu = l;
%                 break;
%             end
%         end
%         for l=3:30
%             if Fx(l,lambdas(j),ks(j)) >=t
%                 xhatt = l;
%                 break;
%             end
%         end
%         V(i) = xhatu;
%         Vtilde(i) = xhatt;
%     end
    U = rand(1,N/2);
    T = 1 - U;
    
    V = wblinv(U,lambdas(j), ks(i));
    Vtilde = wblinv(T,lambdas(j), ks(i));
    W = (P(V) + P(Vtilde))./2;
    tau(j) = mean(W);
    std1(j) = std(W);
end

intervals = zeros(3,month);
intervals(1,:) = tau;
for i=1:month
    intervals(2,i) = -1.96*std1(i)/sqrt(N/2);
    intervals(3,i) = 1.96*std1(i)/sqrt(N/2);
end
intervals = intervals/1e6;
figure
hold on
plot(intervals(1,:))
plot(intervals(1,:)+intervals(2,:),'--r')
plot(intervals(1,:)+intervals(3,:),'--r')
hold off
xticks(1:12)
xticklabels({'Jan','Feb','March','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'})
ylabel("Power output [MW]")
legend("\tau","95% confidence interval")
title("ANTHITHETIC SAMPLING: Power output from wind turbine per month", "N = "+N+", 95% confidence interval")
std1
fprintf("\nUP: STD DOWN: VARIANCE\n")
std1.^2
%% 2e
a = 3;
b = 30;
lambdas = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7]; %Scale
ks = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];             %Shape
month = length(lambdas);

F = @(x,l,k) wblcdf(x, l,k);
prob = zeros(1,12);
for i=1:month
    prob(i) = F(b,lambdas(i),ks(i)) - F(a,lambdas(i),ks(i));
end

bar(prob)
title("Analyical probabiliuty that the power plant generates power")
ylabel("Probability [%]")
xticks(1:12)
xticklabels({'jan','feb','mars','apr','may','jun','jul','aug','sep','oct','nov','dec'})
ylabel("Power output [MW]")

% V = @(l,k) wblrnd(l,k,N,1);
% prob = zeros(1,12);
% for i=1:month
%     v = V(lambdas(i),ks(i));
%     power = P(v);
%     k5 = length(P(v)==0);
%     prob(i) =k5/N;
%     %s3 = sum(v(v<3));
%     %s30 = sum(v(v>30));
%     %prob(i) = (s3+s30)/N;
% end
% 
% bar(1-prob);
figure
hold on
plot(0:0.1:35, wblpdf(0:0.1:35,lambdas(6),ks(6)))
plot(0:0.1:35,P(0:0.1:35)/(1e8))
hold off
legend("Probability density of wind speeds in June", "Power curve for Test D236")
xlabel("x")
ylabel("Probability (for PDF) and normalized power output (for power curve)")
title("Normalized power output and probability density or wind speeds in June")
%% 2f

lambdas = [11.7 10.7 10.1 8.8 8.6 8.9 8.6 8.9 10.0 10.9 11.7 11.7]; %Scale
ks = [2.0 2.0 2.0 1.9 1.9 1.9 1.9 1.9 2.0 1.9 2.0 2.0];             %Shape
month = length(lambdas);

d = 236;
rho = 1.225;
Ptot = @(v) ((rho*pi)/2).*((d^2)/4)*v.^3;
E = zeros(1,12);
Ptot1 = zeros(month,1);
M=3;

for i =1:month
    E(i) = gamma(1+M/(ks(i)))*lambdas(i)^M;
    Ptot1(i) = Ptot(E(i));
end

ratio = (tau./Ptot1');
std3 = zeros(month,1);

for i =1:month
    std3(i)=std1(i)./Ptot1(i);
end

intervals = zeros(3,month);
intervals(1,:) = ratio;

for i=1:month
    intervals(2,i) = -1.96*std3(i)/sqrt(N);
    intervals(3,i) = 1.96*std3(i)/sqrt(N);
end

figure
hold on
plot(1e6*intervals(1,:))
plot((intervals(1,:)+intervals(2,:))*1e6,'--r')
plot((intervals(1,:)+intervals(3,:))*1e6,'--r')
hold off

ratio
I_upper = intervals(1,:) + intervals(3,:)
I_lower = intervals(1,:) + intervals(2,:)
std3'

%% 2g

month = 12;
Prob = zeros(month,1);
for j = 1:month
    Prob(j)=wblcdf(30,lambdas(j),ks(j))-wblcdf(3,lambdas(j),ks(j));
end
bar(1:12,Prob,0.5)

capfac = zeros(month,1);
line1 = zeros(month,1);
line2 = zeros(month,1);
line3 = zeros(month,1);
for i = 1:month
    capfac(i)=tau(i)/(15*1e6);
    line1(i)=0.2;
    line2(i)= 0.4;
    line3(i)= 0.9;
end
figure
hold on
plot(capfac)
plot(line1,'--r')
plot(line2,'--r')
hold off
