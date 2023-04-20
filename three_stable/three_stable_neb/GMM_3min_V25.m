%%
randn('seed',1998) % the random seed
rand('seed',1998) % the random seed
a = 1; b = 1; k0 = 1; S0 = 0.5; l = 4; x0 = 0.5; y0 = 0.5; tmin = 0;
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
    k = 3;
    options = statset('MaxIter', 2000);

    S.mu = [0, 2; 1, 1; 2, 0];
    S.Sigma = zeros(2, 2, k);
    for i = 1:k
        S.Sigma(:, :, i) = 0.1 * diag(ones(1, 2));
    end
    S.ComponentProportion = ones(1, k) / k;
    gmm_2D = fitgmdist(input, k, 'Options', options, 'Start', S);
    gmm = gmm_2D;
    gmm_pdf = @(x, y)reshape(pdf(gmm, [x(:) y(:)]), size(x));
    energy = @(x, y) -log(gmm_pdf(x, y));

    minima_array=gmm.mu;

    s1 = (minima_array(1, 1:2) + minima_array(2, 1:2)) / 2;
    s1 = s1';
    v1 = minima_array(1, 1:2) - minima_array(2, 1:2);
    v1 = v1 / sqrt(sum(v1.^2));
    v1 = v1';

    s2 = (minima_array(2, 1:2) + minima_array(3, 1:2)) / 2;
    s2 = s2';
    v2 = minima_array(2, 1:2) - minima_array(3, 1:2);
    v2 = v2 / sqrt(sum(v2.^2));
    v2 = v2';

    fprintf('########## SADDLE1 ##########\n');
    [saddle1, V1, judge_saddle1] = saddle_GMM_2D(gmm, s1, v1);

    fprintf('########## SADDLE2 ##########\n');
    [saddle2, V2, judge_saddle2] = saddle_GMM_2D(gmm, s2, v2);

    saddle1_pd = gmm_pdf(saddle1(1,1),saddle1(2,1));
    saddle1_e = -log(saddle1_pd);
    saddle2_pd = gmm_pdf(saddle2(1,1),saddle2(2,1));
    saddle2_e = -log(saddle2_pd);

    minimaEnergy=zeros(size(minima_array,1),1)
    for i=1:size(minima_array,1)
        minimaEnergy(i,1)=-log(gmm_pdf(minima_array(i,1),minima_array(i,2)));
    end

    % 计算hession矩阵

    minimaHessian = calculateHessian(gmm,minima_array);
    barrierHessian = calculateHessian(gmm,[saddle1'; saddle2']);

    barrierMinima = [1, 2; 2, 3];
    MFPT = MFPTmatrixWithHessian(minimaEnergy,minimaHessian,[saddle1_e; saddle2_e],barrierHessian,[barrierMinima(:,1) [1:size(barrierMinima(:,1),1)]' barrierMinima(:,2)],1,1/D30);
    MFPT_22 = MFPTmatrixWithReciprocal(minimaEnergy,minimaHessian,[saddle1_e; saddle2_e],barrierHessian,[barrierMinima(:,1) [1:size(barrierMinima(:,1),1)]' barrierMinima(:,2)],1,1/D30);
    fprintf('########## MFPT ##########\n')
    MFPT
    fprintf('########## MFPT4---1/k ##########\n')
    MFPT_22
    % record
    MFPT_static{os,1}=MFPT;
    MFPT_static{os,2}=MFPT_22;
end
save('three_V25_neb.mat','MFPT_static')






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