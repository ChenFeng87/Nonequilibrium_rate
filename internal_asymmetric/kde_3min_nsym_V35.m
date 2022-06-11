%%
fprintf('########## SIMULATION ##########\n');
a = 1; b1 = 0.7; b2 = 1; k0 = 1; S0 = 0.5; l = 4; x0 = 0.5; y0 = 0.5; tmin = 0; 
tmax = 4e7; h = 0.1; V = 35; r =4;
input = euler_simD_r_nsym(a, b1, b2, k0, S0, l, x0, y0, tmin, tmax, h, V, r);

r0 = h * r;
t = tmin:r0:tmax;
t0 = t;

fprintf('########## DIFFUSION COEFFICIENT ##########\n');
[xd, yd, Dkx, Dky, Dkxy, DDx, DDy, DDxy, Dxm, Dym, Dxym, Dkxm, Dkym, Dkxym, DDxm, DDym, DDxym] = evaluateD_2D(input(1:1e7, :), 1/25, r0);
D = (DDxm + DDym) / 2;
D45 = D;

%%
fprintf('########## KDE MODEL ##########\n');
p = kde(input(1:100:100000001, :)', 'rot' );
s_points = [0, 2; 1, 1; 1.5, 0];
minima_array = find_minima_kde(p, s_points, 3);


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
[saddle1, V1] = saddle_kde_2D(p, s1, v1);

fprintf('########## SADDLE2 ##########\n');
[saddle2, V2] = saddle_kde_2D(p, s2, v2);

%%
fprintf('########## PATH & SIMULATION MFPT ##########\n');
num = zeros(1, size(input, 1));

num((input(:, 1) - minima_array(1, 1)).^2 + (input(:, 2) - minima_array(1, 2)).^2 <= 1e-4) = 1;
num((input(:, 1) - minima_array(2, 1)).^2 + (input(:, 2) - minima_array(2, 2)).^2 <= 1e-4) = 2;
num((input(:, 1) - minima_array(3, 1)).^2 + (input(:, 2) - minima_array(3, 2)).^2 <= 1e-4) = 3;

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
way3 = way(way(:, 3) == 3, :);

way1(:, 3) = [];
way2(:, 3) = [];
way3(:, 3) = [];


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
%Create 4 matrices of paths among each 

way1_2 = way1(way1(:, 3) == 2, 1:2);
way2_1 = way2(way2(:, 3) == 1, 1:2);
way2_3 = way2(way2(:, 3) == 3, 1:2);
way3_2 = way3(way3(:, 3) == 2, 1:2);


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
FPT3 = FPT(FPT(:, 3) == 3, :);

FPT1(:, 3) = [];
FPT2(:, 3) = [];
FPT3(:, 3) = [];

nFPT1_2 = FPT1(FPT1(:, 3) == 2, 1:2);
nFPT2_1 = FPT2(FPT2(:, 3) == 1, 1:2);
nFPT2_3 = FPT2(FPT2(:, 3) == 3, 1:2);
nFPT3_2 = FPT3(FPT3(:, 3) == 2, 1:2);

FPT1_2 = zeros(size(nFPT1_2, 1), 1);
FPT2_1 = zeros(size(nFPT2_1, 1), 1);
FPT2_3 = zeros(size(nFPT2_3, 1), 1);
FPT3_2 = zeros(size(nFPT3_2, 1), 1);


for i = 1:size(FPT1_2, 1)
    FPT1_2(i) = t0(nFPT1_2(i, 2)) - t0(nFPT1_2(i, 1));
end

for i = 1:size(FPT2_1, 1)
    FPT2_1(i) = t0(nFPT2_1(i, 2)) - t0(nFPT2_1(i, 1));
end

for i = 1:size(FPT2_3, 1)
    FPT2_3(i) = t0(nFPT2_3(i, 2)) - t0(nFPT2_3(i, 1));
end

for i = 1:size(FPT3_2, 1)
    FPT3_2(i) = t0(nFPT3_2(i, 2)) - t0(nFPT3_2(i, 1));
end


mMFPT1_2 = mean(FPT1_2);
mMFPT2_1 = mean(FPT2_1);
mMFPT2_3 = mean(FPT2_3);
mMFPT3_2 = mean(FPT3_2);


%%
%Divide each path into almost equidistantly into (n+1) points

n = max(2 * max(way(:, 2) - way(:, 1)), 1000);

WAY1_2 = zeros(n + 1, 2, size(way1_2, 1));
WAY2_1 = zeros(n + 1, 2, size(way2_1, 1));
WAY2_3 = zeros(n + 1, 2, size(way2_3, 1));
WAY3_2 = zeros(n + 1, 2, size(way3_2, 1));

%Calculate WAY1_2

for i = 1:size(way1_2, 1)
    tmp_num = way1_2(i, 1):way1_2(i, 2);
    len = 1:(way1_2(i, 2) - way1_2(i, 1));
    
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
            WAY1_2(k + 1, 1, i) = (1 - lambda) * input(tmp_num(j), 1) ...
                + lambda * input(tmp_num(j + 1), 1);
            WAY1_2(k + 1, 2, i) = (1 - lambda) * input(tmp_num(j), 2) ...
                + lambda * input(tmp_num(j + 1), 2);
        end
    end
    WAY1_2(n + 1, 1, i) = input(tmp_num(size(tmp_num, 2)), 1);
    WAY1_2(n + 1, 2, i) = input(tmp_num(size(tmp_num, 2)), 2);
    
end


%Calculate WAY2_1

for i = 1:size(way2_1, 1)
    tmp_num = way2_1(i, 1):way2_1(i, 2);
    len = 1:(way2_1(i, 2) - way2_1(i, 1));
    
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
            WAY2_1(k + 1, 1, i) = (1 - lambda) * input(tmp_num(j), 1) ...
                + lambda * input(tmp_num(j + 1), 1);
            WAY2_1(k + 1, 2, i) = (1 - lambda) * input(tmp_num(j), 2) ...
                + lambda * input(tmp_num(j + 1), 2);
        end
    end
    WAY2_1(n + 1, 1, i) = input(tmp_num(size(tmp_num, 2)), 1);
    WAY2_1(n + 1, 2, i) = input(tmp_num(size(tmp_num, 2)), 2);
    
end


%Calculate WAY2_3

for i = 1:size(way2_3, 1)
    tmp_num = way2_3(i, 1):way2_3(i, 2);
    len = 1:(way2_3(i, 2) - way2_3(i, 1));
    
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
            WAY2_3(k + 1, 1, i) = (1 - lambda) * input(tmp_num(j), 1) ...
                + lambda * input(tmp_num(j + 1), 1);
            WAY2_3(k + 1, 2, i) = (1 - lambda) * input(tmp_num(j), 2) ...
                + lambda * input(tmp_num(j + 1), 2);
        end
    end
    WAY2_3(n + 1, 1, i) = input(tmp_num(size(tmp_num, 2)), 1);
    WAY2_3(n + 1, 2, i) = input(tmp_num(size(tmp_num, 2)), 2);
    
end


%Calculate WAY3_2

for i = 1:size(way3_2, 1)
    tmp_num = way3_2(i, 1):way3_2(i, 2);
    len = 1:(way3_2(i, 2) - way3_2(i, 1));
    
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
            WAY3_2(k + 1, 1, i) = (1 - lambda) * input(tmp_num(j), 1) ...
                + lambda * input(tmp_num(j + 1), 1);
            WAY3_2(k + 1, 2, i) = (1 - lambda) * input(tmp_num(j), 2) ...
                + lambda * input(tmp_num(j + 1), 2);
        end
    end
    WAY3_2(n + 1, 1, i) = input(tmp_num(size(tmp_num, 2)), 1);
    WAY3_2(n + 1, 2, i) = input(tmp_num(size(tmp_num, 2)), 2);
    
end


%%
%Calculate the mean path

mean1_2 = mean(WAY1_2, 3);
mean2_1 = mean(WAY2_1, 3);
mean2_3 = mean(WAY2_3, 3);
mean3_2 = mean(WAY3_2, 3);


%%
%Calculate the saddle point for each mean path

fprintf('########## PATH ENERGY & SADDLE ##########\n');
%probablity density of paths

pd1_2 = zeros(1, size(mean1_2, 1));
pd2_1 = zeros(1, size(mean2_1, 1));
pd2_3 = zeros(1, size(mean2_3, 1));
pd3_2 = zeros(1, size(mean3_2, 1));

%energy of paths

e1_2 = zeros(1, size(mean1_2, 1));
e2_1 = zeros(1, size(mean2_1, 1));
e2_3 = zeros(1, size(mean2_3, 1));
e3_2 = zeros(1, size(mean3_2, 1));


for i = 1:size(mean1_2, 1)
    pd1_2(i) = evaluate(p, [mean1_2(i, 1); mean1_2(i, 2)]);
    e1_2(i) = -log(pd1_2(i));
end

for i = 1:size(mean2_1, 1)
    pd2_1(i) = evaluate(p, [mean2_1(i, 1); mean2_1(i, 2)]);
    e2_1(i) = -log(pd2_1(i));
end

for i = 1:size(mean2_3, 1)
    pd2_3(i) = evaluate(p, [mean2_3(i, 1); mean2_3(i, 2)]);
    e2_3(i) = -log(pd2_3(i));
end

for i = 1:size(mean3_2, 1)
    pd3_2(i) = evaluate(p, [mean3_2(i, 1); mean3_2(i, 2)]);
    e3_2(i) = -log(pd3_2(i));
end

[saddle_pd1_2, saddle_num1_2] = min(pd1_2);
[saddle_pd2_1, saddle_num2_1] = min(pd2_1);
[saddle_pd2_3, saddle_num2_3] = min(pd2_3);
[saddle_pd3_2, saddle_num3_2] = min(pd3_2);

saddle1_pd = evaluate(p, saddle1);
saddle1_e = -log(saddle1_pd);
saddle2_pd = evaluate(p, saddle2);
saddle2_e = -log(saddle2_pd);

%%
%flux based on data
fprintf('########## FLUX ##########\n');
len = 1 / 15;
[xs, ys, pr, flux_x, flux_y, F_data, flux_x1, flux_y1] = flux_2D(input, t0, len);
xf = xs + len / 2;
yf = ys + len / 2;


%%
fprintf('########## LANDSCAPE ##########\n');
[x, y] = meshgrid(0:0.01:2.5, 0:0.01:2.5);
z = zeros(size(x));
for i = 1:size(z, 1)
    z(i, :) = evaluate(p, [x(i, :); y(i, :)]);
end
ez = -log(z);


%%
%probability new
figure(1);
s = surf(x, y, z);
set(s, 'LineStyle','none');
colorbar;
hold on;

pr1 = plot3(mean1_2(:, 1), mean1_2(:, 2), pd1_2, 'r.-');
hold on;
pr2 = plot3(mean2_1(:, 1), mean2_1(:, 2), pd2_1, 'y.-');
hold on;
pr3 = plot3(mean2_3(:, 1), mean2_3(:, 2), pd2_3, 'r.-');
hold on;
pr4 = plot3(mean3_2(:, 1), mean3_2(:, 2), pd3_2, 'y.-');
hold on;

pr5 = plot3(mean1_2(saddle_num1_2, 1), mean1_2(saddle_num1_2, 2), ...
    pd1_2(saddle_num1_2), 'ro', 'markerfacecolor', 'r');
hold on;
pr6 = plot3(mean2_1(saddle_num2_1, 1), mean2_1(saddle_num2_1, 2), ...
    pd2_1(saddle_num2_1), 'yo', 'markerfacecolor', 'y');
hold on;
pr7 = plot3(mean2_3(saddle_num2_3, 1), mean2_3(saddle_num2_3, 2), ...
    pd2_3(saddle_num2_3), 'ro', 'markerfacecolor', 'r');
hold on;
pr8 = plot3(mean3_2(saddle_num3_2, 1), mean3_2(saddle_num3_2, 2), ...
    pd3_2(saddle_num3_2), 'yo', 'markerfacecolor', 'y');
hold on;

pr9 = plot3(saddle1(1), saddle1(2), saddle1_e, 'go', 'markerfacecolor', 'g');
hold on;
pr10 = plot3(saddle2(1), saddle2(2), saddle2_e, 'go', 'markerfacecolor', 'g');
hold on;

name = {'XLim', 'YLim', 'Xtick', 'Ytick'};
value = {[0 2.5], [0 2.5], (0:0.5:2.5), (0:0.5:2.5)};
set(gca, name, value);
xlabel('Protein a');
ylabel('Protein b');
legend([pr1 pr5 pr2 pr6 pr9], {'A-M-B', 'saddle path', 'B-M-A', ...
    'saddle path', 'saddle'}, 'Location', 'north', 'NumColumns', 3);
legend('boxoff')


%%
%probability new_2D
figure(2);
sp = pcolor(x, y, z);
set(sp, 'LineStyle','none');
colorbar;
hold on;

quiver(xf, yf, flux_x, flux_y, 'w');
hold on;

p1 = plot(mean1_2(:, 1), mean1_2(:, 2), 'r.-');
hold on;
p2 = plot(mean2_1(:, 1), mean2_1(:, 2), 'y.-');
hold on;
p3 = plot(mean2_3(:, 1), mean2_3(:, 2), 'r.-');
hold on;
p4 = plot(mean3_2(:, 1), mean3_2(:, 2), 'y.-');
hold on;

p5 = plot(mean1_2(saddle_num1_2, 1), mean1_2(saddle_num1_2, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
p6 = plot(mean2_1(saddle_num2_1, 1), mean2_1(saddle_num2_1, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;
p7 = plot(mean2_3(saddle_num2_3, 1), mean2_3(saddle_num2_3, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
p8 = plot(mean3_2(saddle_num3_2, 1), mean3_2(saddle_num3_2, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;

p9 = plot(saddle1(1), saddle1(2), 'go', 'markerfacecolor', 'g');
hold on;
p10 = plot(saddle2(1), saddle2(2), 'go', 'markerfacecolor', 'g');
hold on;

name = {'XLim', 'YLim', 'Xtick', 'Ytick'};
value = {[0 2.5], [0 2.5], (0:0.5:2.5), (0:0.5:2.5)};
set(gca, name, value);
xlabel('Protein a');
ylabel('Protein b');
legend([p1 p5 p2 p6 p9], {'A-M-B', 'saddle path', 'B-M-A', ...
    'saddle path', 'saddle'}, 'Location', 'north', 'NumColumns', 3);
legend('boxoff')

%%
%landscape new
ez1 = ez;
for i = 1:size(ez, 1)
    for j = 1:size(ez, 2)
        if ez(i, j) >= 7
            ez1(i, j) = 7;
        end
    end
end
figure(4)
sl = surf(x, y, ez1);
set(sl, 'LineStyle','none');
hold on;

o1 = plot3(mean1_2(:, 1), mean1_2(:, 2), e1_2, 'r.-');
hold on;
o2 = plot3(mean2_1(:, 1), mean2_1(:, 2), e2_1, 'y.-');
hold on;
o3 = plot3(mean2_3(:, 1), mean2_3(:, 2), e2_3, 'r.-');
hold on;
o4 = plot3(mean3_2(:, 1), mean3_2(:, 2), e3_2, 'y.-');
hold on;

o5 = plot3(mean1_2(saddle_num1_2, 1), mean1_2(saddle_num1_2, 2), ...
    e1_2(saddle_num1_2), 'ro', 'markerfacecolor', 'r');
hold on;
o6 = plot3(mean2_1(saddle_num2_1, 1), mean2_1(saddle_num2_1, 2), ...
    e2_1(saddle_num2_1), 'yo', 'markerfacecolor', 'y');
hold on;
o7 = plot3(mean2_3(saddle_num2_3, 1), mean2_3(saddle_num2_3, 2), ...
    e2_3(saddle_num2_3), 'ro', 'markerfacecolor', 'r');
hold on;
o8 = plot3(mean3_2(saddle_num3_2, 1), mean3_2(saddle_num3_2, 2), ...
    e3_2(saddle_num3_2), 'yo', 'markerfacecolor', 'y');
hold on;

o9 = plot3(saddle1(1), saddle1(2), saddle1_e, 'go', 'markerfacecolor', 'g');
hold on;
o10 = plot3(saddle2(1), saddle2(2), saddle2_e, 'go', 'markerfacecolor', 'g');
hold on;

name = {'XLim', 'YLim', 'Xtick', 'Ytick'};
value = {[0 2.5], [0 2.5], (0:0.5:2.5), (0:0.5:2.5)};
set(gca, name, value);
xlabel('Protein a');
ylabel('Protein b');
legend([o1 o5 o2 o6 o9], {'A-M-B', 'saddle path', 'B-M-A', ...
    'saddle path', 'saddle'}, 'Location', 'north', 'NumColumns', 3);
legend('boxoff')


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

quiver(xf, yf, flux_x, flux_y, 'w');
hold on;

pe2 = plot(mean1_2(:, 1), mean1_2(:, 2), 'r.-');
hold on;
pe3 = plot(mean2_1(:, 1), mean2_1(:, 2), 'y.-');
hold on;
pe4 = plot(mean2_3(:, 1), mean2_3(:, 2), 'r.-');
hold on;
pe5 = plot(mean3_2(:, 1), mean3_2(:, 2), 'y.-');
hold on;

pe6 = plot(mean1_2(saddle_num1_2, 1), mean1_2(saddle_num1_2, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
pe7 = plot(mean2_1(saddle_num2_1, 1), mean2_1(saddle_num2_1, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;
pe8 = plot(mean2_3(saddle_num2_3, 1), mean2_3(saddle_num2_3, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
pe9 = plot(mean3_2(saddle_num3_2, 1), mean3_2(saddle_num3_2, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;

pe10 = plot(saddle1(1), saddle1(2), 'go', 'markerfacecolor', 'g');
hold on;
pe11 = plot(saddle2(1), saddle2(2), 'go', 'markerfacecolor', 'g');
hold on;

%quiver(0.52, 1.75, 0.5 * cosd(-45), 0.5 * sind(-45), 2, 'r', 'LineWidth', 1);

axis equal;
name = {'XLim', 'YLim', 'Xtick', 'Ytick'};
value = {[0 2.5], [0 2.5], (0:0.5:2.5), (0:0.5:2.5)};
set(gca, name, value);
xlabel('Protein a');
ylabel('Protein b');
legend([pe2 pe6 pe3 pe7 pe10], {'A-M-B', 'saddle path', 'B-M-A', ...
    'saddle path', 'saddle'}, 'Location', 'north', 'NumColumns', 3);
legend('boxoff')
text(0.3, 2.1, 'A', 'FontSize', 13);
text(1.2, 1.2, 'M', 'FontSize', 13);
text(1.7, 0.15, 'B', 'FontSize', 13);

%%
%MFPTnew
fprintf('########## MFPT ##########\n');
saddleCoords = zeros(4, 2);
saddleCoords(1, :) =  mean1_2(saddle_num1_2, :);
saddleCoords(2, :) =  mean2_1(saddle_num2_1, :); 
saddleCoords(3, :) =  mean2_3(saddle_num2_3, :); 
saddleCoords(4, :) =  mean3_2(saddle_num3_2, :); 

%saddleEnergy
saddleEnergy = zeros(1, 4);
for i = 1:4
    saddleEnergy(i) = -log(evaluate(p, [saddleCoords(i, 1); saddleCoords(i, 2)]));
end

%saddleMinima
saddleMinima = [1, 2; 2, 1; 2, 3; 3, 2];


%Calculate mean first passage times with gamma = 1
minimaHessian_kde = calculateHessian_kde_2D(p, minima_array(:, 1:2));
saddleHessian_kde = calculateHessian_kde_2D(p, saddleCoords);
MFPT2 = mMFPTmatrixWithHessian(minima_array(:, 3),minimaHessian_kde,saddleEnergy',saddleHessian_kde,[saddleMinima(:,1) [1:size(saddleMinima(:,1),1)]' saddleMinima(:,2)],1,1/D);
MFPT4 = mMFPTmatrixWithReciprocal(minima_array(:, 3),minimaHessian_kde,saddleEnergy',saddleHessian_kde,[saddleMinima(:,1) [1:size(saddleMinima(:,1),1)]' saddleMinima(:,2)],1,1/D);


%%
%MFPTNEB
barrierMinima = [1, 2; 2, 3];
barrierHessian_kde = calculateHessian_kde_2D(p, [saddle1'; saddle2']);
MFPT1 = MFPTmatrixWithHessian(minima_array(:, 3),minimaHessian_kde,[saddle1_e; saddle2_e],barrierHessian_kde,[barrierMinima(:,1) [1:size(barrierMinima(:,1),1)]' barrierMinima(:,2)],1,1/D);
MFPT3 = MFPTmatrixWithReciprocal(minima_array(:, 3),minimaHessian_kde,[saddle1_e; saddle2_e],barrierHessian_kde,[barrierMinima(:,1) [1:size(barrierMinima(:,1),1)]' barrierMinima(:,2)],1,1/D);


%%
MFPT45_kde = {MFPT2,MFPT4,mMFPT1_2,mMFPT2_1,mMFPT2_3,mMFPT3_2,MFPT1,MFPT3};
1/D
sound(sin(2*pi*25*(1:4000)/100));