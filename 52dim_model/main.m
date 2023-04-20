% clear all; clc;
% addpath('/share/inspurStorage/home1/chenfeng/matla_package/code_figure')
% addpath('/home/biansr/chenfeng/code_figure')
addpath('D:\projects\baiyubo_matlab\code_figure')
randn('seed',1998) % the random seed
rand('seed',1998) % the random seed
MFPT_static=cell(20,4)
saddle_static=cell(20,2)
%%parameters of model 
IM=load ('matrixCF.txt');
N = size(IM,1); %The dimension of the system.
nill=zeros(N,N); %Hill functiom
S=zeros(N,N); % threshold
AB=zeros(N,N); %the scale factors for the activation and inhibition
for i=1:N
    for j=1:N
        if IM(i,j)==1
            AB(i,j)=0.37;
            S(i,j)=0.5;
            nill(i,j)=3;
        elseif IM(i,j)==-1
            AB(i,j)=0.5;
            S(i,j)=0.5;
            nill(i,j)=3;
        end
    end
end
K=ones(N,1);
tmin = 0; 
tmax = 2e5; % 2e5
h = 0.005; V_dd = 0.0275; %%%%%%%%%%%%% parameters for step and D
t0 = tmin:h:tmax;
all_length = tmax/h;
save_step = 100000; % 100000
save_number = all_length/save_step;
for os=1:1
%     x0=ones(N,1)*0.5;
%     fprintf('number=%d',os);
%     fprintf('########## SIMULATION ##########\n');  
%     input=zeros(all_length+1,2);
%     for ot=1:save_number
%         input_orial = euler_simD_r(IM,AB,S,nill,K, x0, save_step, h, V_dd); % generate a long trajectory
%         input(save_step*ot-save_step+1:save_step*ot+1,1)=input_orial(:,16); % Gata6
%         input(save_step*ot-save_step+1:save_step*ot+1,2)=input_orial(:,3); % Nonog
%         x0=input_orial(save_step+1,:)';
%         clear input_orial
%     end

%%
    fprintf('########## KDE MODEL ##########\n');
    p = kde(input(1:10:40000001, :)', 'rot' ); % 40000001
    s_points = [2, 0; 1, 1];
    minima_array = find_minima_kde(p, s_points, 2);
    
    % 判断稳定点个数
    if size(minima_array,1) ~= 2
        continue;
    end

    s0 = (minima_array(1, 1:2) + minima_array(2, 1:2)) / 2;
    s0 = s0';
    v0 = minima_array(1, 1:2) - minima_array(2, 1:2);
    v0 = v0 / sqrt(sum(v0.^2));
    v0 = v0';
    %%
    fprintf('########## SADDLE ##########\n');
    [saddle, v, judge_saddle] = saddle_kde_2D(p, s0, v0);  %ZiMo:从稳定点连线中点找鞍点，不理解原理，好像来自一篇文献，包括前面从数据中得到flux也是来自某文献
   
    % 判断鞍点有没有找到
    if judge_saddle ~= 0
        continue;
    end
    %%
    fprintf('########## PATH & SIMULATION MFPT ##########\n');
    num = zeros(1, size(input, 1));

    num((input(:, 1) - minima_array(1, 1)).^2 + (input(:, 2) - minima_array(1, 2)).^2 <= 1e-4) = 1;
    num((input(:, 1) - minima_array(2, 1)).^2 + (input(:, 2) - minima_array(2, 2)).^2 <= 1e-4) = 2;

    index = zeros(2, size(input, 1));
    index(1, :) = 1:size(input, 1);
    index(2, :) = num;

    id = true(1, size(input, 1));
    for i = 1:size(input, 1)
        if num(i) == 0
            id(i) = false;
        end
    end

    index = index(:, id);
    tindex = index;
    ind = true(1, size(index, 2));
    for i = 2:(size(index, 2) - 1)
        if index(2, i) == index(2, i - 1) && index(2, i) == index(2, i + 1)
            ind(i) = false;
        end
    end
    if index(2, 1) == index(2, 2)
        ind(1) = false;
    end
    if index(2, size(index, 2)) == index(2, size(index, 2) - 1)
        ind(size(index, 2)) = false;
    end

    index = index(:, ind);


    %Create a way matrix, the 1st column is the number of the starting point,
    %the 2nd column is the number of ending point, the 3rd column is state
    %number of the starting point, the 4th column is state number of the ending
    %point.
    way = zeros(size(index, 2) - 1, 4);

    for i = 1:(size(index, 2) - 1)
        if index(2, i) ~= index(2, i + 1)
            way(i, 1) = index(1, i);
            way(i, 2) = index(1, i + 1);
            way(i, 3) = index(2, i);
            way(i, 4) = index(2, i + 1);
        end
    end

    way = way(way(:, 1) ~= 0, :);
    way1 = way(way(:, 3) == 1, :);
    way2 = way(way(:, 3) == 2, :);
    way1(:, 3) = [];
    way2(:, 3) = [];


    %%
    %index for calculate mfpt
    tind = true(1, size(tindex, 2));
    for i = 2:size(tindex, 2)
        if tindex(2, i) == tindex(2, i - 1)
            tind(i) = false;
        end
    end

    tindex = tindex(:, tind);

    %%
    %Create 2 matrices of paths among each 
    way_1to2 = way1(way1(:, 3) == 2, 1:2);
    way_2to1 = way2(way2(:, 3) == 1, 1:2);


    %%
    %Record FPT  and calculate Measured MFPT
    FPT = zeros(size(tindex, 2) - 1, 4);

    for i = 1:(size(tindex, 2) - 1)
        if tindex(2, i) ~= tindex(2, i + 1)
            FPT(i, 1) = tindex(1, i);
            FPT(i, 2) = tindex(1, i + 1);
            FPT(i, 3) = tindex(2, i);
            FPT(i, 4) = tindex(2, i + 1);
        end
    end

    FPT = FPT(FPT(:, 1) ~= 0, :);

    FPT1 = FPT(FPT(:, 3) == 1, :);
    FPT2 = FPT(FPT(:, 3) == 2, :);

    FPT1(:, 3) = [];
    FPT2(:, 3) = [];

    nFPT1_2 = FPT1(FPT1(:, 3) == 2, 1:2);
    nFPT2_1 = FPT2(FPT2(:, 3) == 1, 1:2);

    FPT1_2 = zeros(size(nFPT1_2, 1), 1);
    FPT2_1 = zeros(size(nFPT2_1, 1), 1);

    for i = 1:size(FPT1_2, 1)
        FPT1_2(i) = t0(nFPT1_2(i, 2)) - t0(nFPT1_2(i, 1));
    end

    for i = 1:size(FPT2_1, 1)
        FPT2_1(i) = t0(nFPT2_1(i, 2)) - t0(nFPT2_1(i, 1));
    end

    mMFPT1_2 = mean(FPT1_2);
    mMFPT2_1 = mean(FPT2_1);


    %%
    %Divide each path into almost equidistantly into (n+1) points

    n = max(2 * max(way(:, 2) - way(:, 1)), 500);
    WAY_1to2 = zeros(n + 1, 2, size(way_1to2, 1));
    WAY_2to1 = zeros(n + 1, 2, size(way_2to1, 1));

    %Calculate WAY_1to2
    for i = 1:size(way_1to2, 1)
        tmp_num = way_1to2(i, 1):way_1to2(i, 2);
        len = 1:(way_1to2(i, 2) - way_1to2(i, 1));

        for j = 1:size(len, 2)
            len(j) = sqrt((input(tmp_num(j), 1) - ...
                input(tmp_num(j + 1), 1))^2 + (input(tmp_num(j), 2) - ...
                input(tmp_num(j + 1), 2))^2);
        end

        len_per = floor(n * len / sum(len));
        rem = n - sum(len_per);

        if rem > 0
            for j = 1:rem
                len_per(j) = len_per(j) + 1;
            end
        end

        len_num = zeros(1, size(len_per, 2) + 1);
        for j = 1:size(len_per, 2)
            len_num(j + 1) = sum(len_per(1:j));
        end

        for j = 1:size(len_per, 2)
            for k = len_num(j):(len_num(j+1) - 1)
                lambda = (k - len_num(j)) / (len_num(j + 1) - len_num(j));
                WAY_1to2(k + 1, 1, i) = (1 - lambda) * input(tmp_num(j), 1) ...
                    + lambda * input(tmp_num(j + 1), 1);
                WAY_1to2(k + 1, 2, i) = (1 - lambda) * input(tmp_num(j), 2) ...
                    + lambda * input(tmp_num(j + 1), 2);
            end
        end
        WAY_1to2(n + 1, 1, i) = input(tmp_num(size(tmp_num, 2)), 1);
        WAY_1to2(n + 1, 2, i) = input(tmp_num(size(tmp_num, 2)), 2);

    end

    %Calculate WAY_2to1
    for i = 1:size(way_2to1, 1)
        tmp_num = way_2to1(i, 1):way_2to1(i, 2);
        len = 1:(way_2to1(i, 2) - way_2to1(i, 1));

        for j = 1:size(len, 2)
            len(j) = sqrt((input(tmp_num(j), 1) - ...
                input(tmp_num(j + 1), 1))^2 + (input(tmp_num(j), 2) - ...
                input(tmp_num(j + 1), 2))^2);
        end

        len_per = floor(n * len / sum(len));
        rem = n - sum(len_per);

        if rem > 0
            for j = 1:rem
                len_per(j) = len_per(j) + 1;
            end
        end

        len_num = zeros(1, size(len_per, 2) + 1);
        for j = 1:size(len_per, 2)
            len_num(j + 1) = sum(len_per(1:j));
        end

        for j = 1:size(len_per, 2)
            for k = len_num(j):(len_num(j+1) - 1)
                lambda = (k - len_num(j)) / (len_num(j + 1) - len_num(j));
                WAY_2to1(k + 1, 1, i) = (1 - lambda) * input(tmp_num(j), 1) ...
                    + lambda * input(tmp_num(j + 1), 1);
                WAY_2to1(k + 1, 2, i) = (1 - lambda) * input(tmp_num(j), 2) ...
                    + lambda * input(tmp_num(j + 1), 2);
            end
        end
        WAY_2to1(n + 1, 1, i) = input(tmp_num(size(tmp_num, 2)), 1);
        WAY_2to1(n + 1, 2, i) = input(tmp_num(size(tmp_num, 2)), 2);

    end


    %%
    %Calculate the mean path
    mean_1to2 = mean(WAY_1to2, 3);
    mean_2to1 = mean(WAY_2to1, 3);

    %%
    %Calculate the saddle point for each mean path

    fprintf('########## PATH ENERGY & SADDLE ##########\n');
    %probablity density of paths

    pd_1to2 = zeros(1, size(mean_1to2, 1));
    pd_2to1 = zeros(1, size(mean_2to1, 1));

    %energy of paths

    e_1to2 = zeros(1, size(mean_1to2, 1));
    e_2to1 = zeros(1, size(mean_2to1, 1));

    for i = 1:size(mean_1to2, 1)
        pd_1to2(i) = evaluate(p, [mean_1to2(i, 1); mean_1to2(i, 2)]);
        e_1to2(i) = -log(pd_1to2(i));
    end

    for i = 1:size(mean_2to1, 1)
        pd_2to1(i) = evaluate(p, [mean_2to1(i, 1); mean_2to1(i, 2)]);
        e_2to1(i) = -log(pd_2to1(i));
    end

    [saddle_pd_1to2, saddle_num_1to2] = min(pd_1to2);
    [saddle_pd_2to1, saddle_num_2to1] = min(pd_2to1);

    saddle_pd = evaluate(p, saddle);
    saddle_e = -log(saddle_pd);
    %%
    %MFPTnew
    fprintf('########## MFPT ##########\n');
    saddleCoords = zeros(2, 2);
    saddleCoords(1, :) =  mean_1to2(saddle_num_1to2, :);
    saddleCoords(2, :) =  mean_2to1(saddle_num_2to1, :); 

    %saddleEnergy
    saddleEnergy = zeros(1, 2);
    for i = 1:2
        saddleEnergy(i) = -log(evaluate(p, [saddleCoords(i, 1); saddleCoords(i, 2)]));
    end

    %saddleMinima
    saddleMinima = [1, 2; 2, 1];

    minimaHessian_kde = calculateHessian_kde_2D(p, minima_array(:, 1:2));
    saddleHessian_kde = calculateHessian_kde_2D(p, saddleCoords);
    [MFPT2,judge_1] = mMFPTmatrixWithHessian(minima_array(:, 3),minimaHessian_kde,saddleEnergy',saddleHessian_kde,[saddleMinima(:,1) [1:size(saddleMinima(:,1),1)]' saddleMinima(:,2)],1,1/V_dd);
    
    % 判断特征值符不符合条件
    if judge_1 ~= 0
        continue;
    end
    %%
    %MFPTNEB
    barrierMinima = [1, 2];
    barrierHessian_kde = calculateHessian_kde_2D(p, saddle');
    [MFPT_NEB,judge_2] = MFPTmatrixWithHessian(minima_array(:, 3),minimaHessian_kde,saddle_e,barrierHessian_kde,[barrierMinima(:,1) [1:size(barrierMinima(:,1),1)]' barrierMinima(:,2)],1,1/V_dd);
    if judge_2 ~= 0
        continue;
    end

%% picture
%flux based on data
fprintf('########## FLUX ##########\n');
len = 1 / 15;
%len = d;
[xs, ys, pr, flux_x, flux_y, F_data, flux_x1, flux_y1] = flux_2D(input, t0, len);
xf = xs + len / 2;
yf = ys + len / 2;

fprintf('########## LANDSCAPE ##########\n');
[x, y] = meshgrid(0:0.01:3.0, 0:0.01:2.0);
z = zeros(size(x));
for i = 1:size(z, 1)
    z(i, :) = evaluate(p, [x(i, :); y(i, :)]);
end
ez = -log(z);
%%
%landscape new_2D
ez1 = ez;
for i = 1:size(ez, 1)
    for j = 1:size(ez, 2)
        if ez(i, j) >= 7
            ez1(i, j) = 7;
        end
    end
end
figure(5);
se = pcolor(x, y, ez1);
set(se, 'LineStyle','none');
colorbar;
hold on;
cs1= 1:2:size(xf,1);
quiver(xf(cs1,cs1), yf(cs1,cs1), flux_x(cs1,cs1), flux_y(cs1,cs1), 'w', 'linewidth', 0.8);
hold on;

pe2 = plot(mean_1to2(:, 1), mean_1to2(:, 2), 'r.-');
hold on;
pe3 = plot(mean_2to1(:, 1), mean_2to1(:, 2), 'y.-');
hold on;

pe6 = plot(mean_1to2(saddle_num_1to2, 1), mean_1to2(saddle_num_1to2, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
pe7 = plot(mean_2to1(saddle_num_2to1, 1), mean_2to1(saddle_num_2to1, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;

pe10 = plot(saddle(1), saddle(2), 'go', 'markerfacecolor', 'g');
hold on;

axis equal;
name = {'XLim', 'YLim', 'Xtick', 'Ytick'};
value = {[0 3.0], [0 2.0], (0:0.5:3.0), (0:0.5:2.0)};
set(gca, name, value);
xlabel('Protein a');
ylabel('Protein b');

figure(4)
sl = surf(x, y, ez1);
set(sl, 'LineStyle','none');
hold on;

pl1 = plot3(mean_1to2(:, 1), mean_1to2(:, 2), e_1to2, 'r.-');
hold on;
pl2 = plot3(mean_2to1(:, 1), mean_2to1(:, 2), e_2to1, 'y.-');
hold on;
pl3 = plot3(mean_1to2(saddle_num_1to2, 1), mean_1to2(saddle_num_1to2, 2), ...
    e_1to2(saddle_num_1to2)+0.05, 'ro', 'markerfacecolor', 'r');
hold on;
pl4 = plot3(mean_2to1(saddle_num_2to1, 1), mean_2to1(saddle_num_2to1, 2), ...
    e_2to1(saddle_num_2to1)+0.05, 'yo', 'markerfacecolor', 'y');
hold on;
pl5 = plot3(saddle(1), saddle(2), saddle_e+0.05, 'go', 'markerfacecolor', 'g');
hold on;

% axis equal;
name = {'XLim', 'YLim', 'Xtick', 'Ytick'};
value = {[0 3.0], [0 2.0], (0:0.5:3.0), (0:0.5:2.0)};
set(gca, name, value);
xlabel('GATA6');
ylabel('NANOG');

    
%%
    %record
    MFPT30_kde = {MFPT2, mMFPT1_2, mMFPT2_1, MFPT_NEB};
    fprintf('########## minima_location ##########\n')
    minima_array
    fprintf('########## MFPT2 ##########\n')
    MFPT2
    fprintf('########## mMFPT1_2 ##########\n')
    mMFPT1_2
    fprintf('########## mMFPT2_1 ##########\n')
    mMFPT2_1
    fprintf('########## MFPT_NEB ##########\n')
    MFPT_NEB
    
    saddle_static{os,1}=saddleCoords;
    saddle_static{os,2}=saddle;
    
    MFPT_static{os,1}=MFPT2;
    MFPT_static{os,2}=mMFPT1_2;
    MFPT_static{os,3}=mMFPT2_1;
    MFPT_static{os,4}=MFPT_NEB;
end
save('two_D0275_wai.mat','MFPT_static','saddle_static')
