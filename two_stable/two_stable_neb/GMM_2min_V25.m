%%
randn('seed',1998) % the random seed
rand('seed',1998) % the random seed
a = 0; b = 1; k0 = 1; S0 = 0.5; l = 4; x0 = 0.5; y0 = 0.5; tmin = 0;
tmax = 4e5; h = 0.001; V = 25; r = 1;  % V=15 tmax = 4e7
r0 = h * r;
t = tmin:r0:tmax;
t0 = t;

MFPT_static=cell(20,2)

for os=1:20
    fprintf('number=%d',os); 
    fprintf('########## SIMULATION ##########\n');
    input = euler_simD_r(a, b, k0, S0, l, x0, y0, tmin, tmax, h, V, r);

    fprintf('########## DIFFUSION COEFFICIENT ##########\n'); % ZiMo: input 大小相应改变
    [xd, yd, Dkx, Dky, Dkxy, DDx, DDy, DDxy, Dxm, Dym, Dxym, Dkxm, Dkym, Dkxym, DDxm, DDym, DDxym] = evaluateD_2D(input(1:2e7, :), 1/25, r0);
    D = (DDxm + DDym) / 2;
    D30 = D;

    fprintf('########## GMM MODEL ##########\n');
    k = 2;
    options = statset('MaxIter', 2000);

    S.mu = [0, 1; 1, 0];
    S.Sigma = zeros(2, 2, k);
    for i = 1:k
        S.Sigma(:, :, i) =  diag(ones(1, 2));
    end
    S.ComponentProportion = ones(1, k) / k;
    gmm_2D = fitgmdist(input, k, 'Options', options, 'Start', S);
    gmm = gmm_2D;
    gmm_pdf = @(x, y)reshape(pdf(gmm, [x(:) y(:)]), size(x));
    energy = @(x, y) -log(gmm_pdf(x, y));

    minima_array=gmm.mu;

    s0 = (minima_array(1, 1:2) + minima_array(2, 1:2)) / 2;
    s0 = s0';
    v0 = minima_array(1, 1:2) - minima_array(2, 1:2);
    v0 = v0 / sqrt(sum(v0.^2));
    v0 = v0';

    fprintf('########## SADDLE1 ##########\n');
    [saddle, v, judge_saddle] = saddle_GMM_2D(gmm, s0, v0);
     % 判断鞍点有没有找到
    if judge_saddle ~= 0
        continue;
    end

    saddle_pd = gmm_pdf(saddle(1,1),saddle(2,1));
    saddle_e = -log(saddle_pd);

    minimaEnergy=zeros(size(minima_array,1),1)
    for i=1:size(minima_array,1)
        minimaEnergy(i,1)=-log(gmm_pdf(minima_array(i,1),minima_array(i,2)));
    end

    % 计算hession矩阵

    minimaHessian = calculateHessian(gmm,minima_array);
    barrierHessian = calculateHessian(gmm,saddle');

    barrierMinima = [1, 2];
    MFPT = MFPTmatrixWithHessian(minimaEnergy,minimaHessian,saddle_e,barrierHessian,[barrierMinima(:,1) [1:size(barrierMinima(:,1),1)]' barrierMinima(:,2)],1,1/D30);
    MFPT_22 = MFPTmatrixWithReciprocal(minimaEnergy,minimaHessian,saddle_e,barrierHessian,[barrierMinima(:,1) [1:size(barrierMinima(:,1),1)]' barrierMinima(:,2)],1,1/D30);
    fprintf('########## MFPT ##########\n')
    MFPT
    fprintf('########## MFPT4---1/k ##########\n')
    MFPT_22
    % record
    MFPT_static{os,1}=MFPT;
    MFPT_static{os,2}=MFPT_22;
end
save('two_V25_neb.mat','MFPT_static')






% [x, y] = meshgrid(0:0.01:2.5, 0:0.01:2.5);
% z = zeros(size(x));
% for i = 1:size(z, 1)
%     for j=1:size(z,2)
%         z(i, j) = gmm_pdf(x(i, j), y(i, j));
%     end
% end
% ez = -log(z);
% ez1 = ez;
% for i = 1:size(ez, 1)
%     for j = 1:size(ez, 2)
%         if ez(i, j) >= 7
%             ez1(i, j) = 7;
%         end
%     end
% end
% figure(1);
% s = surf(x, y,ez1);
% set(s, 'LineStyle','none');
% hold on;
% pr10 = plot3(saddle1(1), saddle1(2), saddle1_e, 'go', 'markerfacecolor', 'g');
% hold on;
% pr11 = plot3(saddle2(1), saddle2(2), saddle2_e, 'go', 'markerfacecolor', 'g');
% hold on;