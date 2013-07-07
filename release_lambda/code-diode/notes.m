
Original code:

mu = [1 2;-3 -5];
sigma = cat(3,[2 0;0 .5],[1 0;0 1]);
p = ones(1,2)/2;
obj = gmdistribution(mu,sigma,p);

ezsurf(@(x,y)pdf(obj,[x y]),[-10 10],[-10 10])

(0.063, 0.758) and (0.884, 0.21)

mu = [0.063 0.758; 0.884 0.21];
sigma = cat(3,[ 0.0881 -0.0013;-0.0013 0.1199],[0.0022 0.0122;0.0122 0.0759]);
p = ones(1,2)/2;
obj = gmdistribution(mu,sigma,p);

ezsurf(@(x,y)pdf(obj,[x y]),[-0.4 1.6],[-0.4 1.6])



    0.0881   -0.0013
   -0.0013    0.1199


ans(:,:,2) =

    0.0022    0.0122
    0.0122    0.0759




Sample generation
MU1 = [1 2];
SIGMA1 = [2 0; 0 .5];
MU2 = [-3 -5];
SIGMA2 = [1 0; 0 1];
X = [mvnrnd(MU1,SIGMA1,1000);
mvnrnd(MU2,SIGMA2,1000)];
scatter(X(:,1),X(:,2),10,'.')


fitting


options = statset('Display','final');
obj = gmdistribution.fit(X,2,'Options',options);
hold on
h = ezcontour(@(x,y)pdf(obj,[x y]),[-8 6],[-8 6]);
hold off


ComponentMeans = obj.mu
ComponentMeans =
    0.9391    2.0322
   -2.9823   -4.9737

ComponentCovariances = obj.Sigma
ComponentCovariances(:,:,1) =
    1.7786   -0.0528
   -0.0528    0.5312
ComponentCovariances(:,:,2) =
    1.0491   -0.0150
   -0.0150    0.9816

MixtureProportions = obj.PComponents
MixtureProportions =
    0.5000    0.5000



    he two-component model minimizes the Akaike information:

AIC = zeros(1,4);
obj = cell(1,4);
for k = 1:4
    obj{k} = gmdistribution.fit(X,k);
    AIC(k)= obj{k}.AIC;
end

[minAIC,numComponents] = min(AIC);
numComponents
numComponents =
     2

model = obj{2}
model = 
Gaussian mixture distribution
with 2 components in 2 dimensions
Component 1:
Mixing proportion: 0.500000
Mean:     0.9391    2.0322
Component 2:
Mixing proportion: 0.500000
Mean:    -2.9823   -4.9737
Both the Akaike and Bayes information are negative log-likelihoods for the data with penalty terms for the number of estimated parameters. You can use them to determine an appropriate number of components for a model when the number of components is unspecified.

Simulating Gaussian Mixtures


ans(:,:,1) =

    0.0881   -0.0013
   -0.0013    0.1199


ans(:,:,2) =

    0.0022    0.0122
    0.0122    0.0759



     I'm looking for an efficient algorithm (i.e. Matlab-efficient) to
> sample N samples from a multivariate Gaussian mixture model that has
> M component densities.


Hi Rudolph -

When you say, "multivariate gaussian mixture", I'll assume that you mean that 
each component is a MVN.

There's lots of ways to do this, including a loop over 1:N. The following code 
is just to prove that it can be vectorized (the loop is only over 1:M, which 
ought to be short). There's no error checking in this code, by the way. There 
are several ways to do the "stacked matrix multiply", and in fact a loop over 
1:N may be the fastest, depending on N, M, and d. The first way I did it uses 
more memory to avoid the 1:N loop. The second uses FIND and a 1:M loop. 
There's also a sparse solution.

By the way, in the next release, MVNRND will accepts multiple mean and cov 
inputs, so this becomes a two-liner.

Hope this helps.

- Peter Perkins
   The MathWorks, Inc.


% some test parameters
M = 3;
mu = [-5 -5; 0 0; 5 5];
sigma = cat(3, [2 0; 0 1], [2 -.2; -.2 2], [1 .9; .9 1]);
y = mvgmmrnd(mu,sigma,[1 1 1],100000);

function y = mvgmmrnd(mu,sigma,p,n)
%MVGMMRND Random vectors from a mixture of multivariate normals.
% MU is an M-by-D matrix of means for the M component normals
% SIGMA is a D-by-D-by-M array of covariance matrices for the
% M component normals.
% P is an M-by-1 vector of component mixing probabilities.
% N is the desired number of random vectors.

[M,d] = size(mu);

% randomly pick from the components
[dum,compon] = histc(rand(n,1), [0; cumsum(p(:))./sum(p)]);

% generate random vectors from the selected components with a
% "stacked" matrix multiply
for i = 1:M
     Rt(i,:,:) = chol(sigma(:,:,i)); % holds the transposed cholesky factors
end
Z = repmat(randn(n,d), [1,1,d]);
y = squeeze(sum(Z.*Rt(compon,:,:),2)) + mu(compon,:);

% another way to generate the random vectors
% y = zeros(n,d);
% for i = 1:M
% mbrs = find(compon == i);
% ni = length(mbrs);
% y(mbrs,:) = randn(ni,d) * chol(sigma(:,:,i)) + repmat(mu(i,:),ni,1);
% end