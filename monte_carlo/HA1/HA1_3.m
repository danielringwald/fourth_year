%% Load
load powercurve_D236.mat

%% Define funcs
l = 10.05;
k = 1.95;

F = @(x) wblcdf(x,l,k);
f = @(x) wblpdf(x,l,k);

alpha = 0.638;
p = 3;
q = 1.5;

F12 = @(v1, v2) F(v1).*F(v2).*(1 + alpha.*(1 - F(v1).^p).^q .* (1-F(v2).^p).^q);
f12 = @(v1, v2) f(v1).*f(v2).*( 1 + alpha.*(1-F(v1).^p).^(q-1).*(1-F(v2).^p).^(q-1).*...
    (F(v1).^p.*(1+p*q)-1).*(F(v2).^p.*(1+p*q)-1) );


%% 3a Importance sampling to calc. E(P(V1)+P(V2))

x = linspace(0,30);
N = 5000;

phif = @(x)P(x)'.*f(x);
[val, idx] = max(phif(x));
mu = x(idx);
sig = 5;
draw = normrnd(mu, sig,1,N);

g = @(x) normpdf(x, mu, sig);

phiomega = phif(draw)./g(draw);
tau = 2*mean(phiomega);

fprintf("\nExpectation value of P(V_1) + P(V_2): %f MW\n", tau/1e6);

%% 3b
N = 10000;
x = linspace(0,30);
y = linspace(0,30);
f_grid = zeros(length(x));

[X,Y] = meshgrid(x,y);
for i=1:length(x)
    for j=1:length(y)
        f_grid(i,j) = f12(x(i),y(j));
    end
end

mus = [8,8];
sigmas = [23 0;0 23];
multi_randf = mvnpdf([X(:),Y(:)],mus,sigmas);

phiomega = zeros(1,N);
for i=1:N
    u = mvnrnd(mus,sigmas);
    omega = f12(u(1),u(2))./mvnpdf(u,mus,sigmas);

    phiomega(i) = P(u(1)).*P(u(2))*omega;
end

tau_mv = mean(phiomega);
cov = tau_mv - (tau/2)^2;
figure 
hold on
surf(X,Y,f_grid)
surf(X,Y,reshape(multi_randf,length(x),length(x))/norm(reshape(multi_randf,length(x),length(x))))
hold off

fprintf("\nCov(P(V_1,V_2)): %f\n", cov)

%% 3c

x = linspace(0,30);
N = 5000;

phif = @(x) P(x).^2'.*f(x);
[val, idx] = max(phif(x));
mu = x(idx);
sig = 5;
draw = normrnd(mu, sig,1,N);

g = @(x) normpdf(x, mu, sig);

phiomega = phif(draw)./g(draw);

tau2 = 2*mean(phiomega);
variability = 2*(tau^2)-2*tau2+2*cov;

fprintf("\nVariability, V(P(V_1)+P(V_2)): %f\n", variability)
fprintf("\nStandard devation, V(P(V_1)+P(V_2)): %f\n", sqrt(variability))

%% 3d

x = linspace(0,30);
N = 5000;

phif = @(x)P(x)'.*f(x);
[val, idx] = max(phif(x));
mu = x(idx);
sig = 5;
draw = normrnd(mu, sig,1,N);

g = @(x) normpdf(x, mu, sig);

phiomega = phif(draw)./g(draw);
phiomega2 = phiomega + phiomega;

Probg = sum(phiomega2> 15)/N;
Probl = sum(phiomega2 < 15)/N;

fprintf("\nProb. > 15M W: %f\nProb. < 15 MW: %f\n",Probg,Probl)


