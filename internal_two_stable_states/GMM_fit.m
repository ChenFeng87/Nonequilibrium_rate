%%
%input = simdata(:, [2 7]);
randn('seed',2021)
%Constant diffusion coefficient ？？？？Zimo不理解 为啥是常数扩散系数
fprintf('########## SIMULATION ##########\n');  
a = 0; b = 1.0; k0 = 1; S0 = 0.5; l = 4; x0 = 0.5; y0 = 0.5; tmin = 0; 
tmax = 1e4; h = 0.1; V = 35; r = 4;  % tmax 原先是8e7，Zimo改为8e5
input = euler_simD_r(a, b, k0, S0, l, x0, y0, tmin, tmax, h, V, r);

r0 = h * r;
t = tmin:r0:tmax;
t0 = t;
%%
fit GMM for 2D
fprintf('########## GMM MODEL ##########\n');
k = 2;
options = statset('MaxIter', 2000);

%set start

S.mu = [0, 1; 1, 0];
%S.Sigma = ones(1, 2, k);
S.Sigma = zeros(2, 2, k);
for i = 1:k
    S.Sigma(:, :, i) = 0.1 * diag(ones(1,2));
end
S.ComponentProportion = ones(1, k) / k;

gmm_2D = fitgmdist(input, k, 'Options', options, 'Start', S);

gmm = gmm_2D;

plot surf
fsurf(@(x, y)-log(reshape(pdf(gmm_2D, [x(:) y(:)]), size(x))), [-0.5 2])

