

% Fix stream of random numbers
s1 = RandStream.create('mrg32k3a','Seed', 22);
s0 = RandStream.setGlobalStream(s1);
% 
% % Create data set
% n = 100; p = 6;
% correlation = 0.6;
% Sigma = correlation*ones(p) + (1 - correlation)*eye(p);
% mu = zeros(p,1);
% X = mvnrnd(mu, Sigma, n);
% % Model is lin.comb. of first three variables plus noise
% y = X(:,1) + 0.05.*X(:,2) + X(:,3) + 2*randn(n,1);
% 
% % Preprocess data
% X = normalize(X);
% y = center(y);

%t=3e-4;
% k=100
%if y(3e-4) is dep. 

%check y=x1 is dep on X2, X3 and X4
%construct y with K samples
%X2-X3,X4 at time t: K samples
%X2-X3,X4 at time t-dt: K samples
%X2-X3,X4 at time t-2dt: K samples
% 30 input variables.

[a,b,c]=size(T3D)
X=zeros(a*c,b);   %%%10*100*3
for i=1:a
    for j=1:c
    X(c*(i-1)+j,:)=T3D(i,j);
    end
end
X=X'
y=T2D(1,:)'
% Run LAR
[beta info] = lar(X, y, 0, true, true);

% Find best fitting model
[bestAIC bestIdx] = min(info.AIC);
best_s = info.s(bestIdx);

% Plot results
h1 = figure(1);
plot(info.s, beta, '.-');
xlabel('s'), ylabel('\beta', 'Rotation', 0)
line([best_s best_s], [-6 14], 'LineStyle', ':', 'Color', [1 0 0]);
legend('1','2','3','4','5','6',2);

% Restore random stream
RandStream.setGlobalStream(s0);
