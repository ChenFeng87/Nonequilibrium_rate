%%
randn('seed',2022)
fprintf('########## SIMULATION ##########\n');
a = 0.95; b = 0.05; k0 = 1; S0 = 0.5; l = 4; x0 = 0.5; y0 = 0.5; tmin = 0; 
tmax = 4e7; h = 0.1; V = 35; r = 4;
input = euler_simD_r(a, b, k0, S0, l, x0, y0, tmin, tmax, h, V, r);

r0 = h * r;
t = tmin:r0:tmax;
t0 = t;

fprintf('########## DIFFUSION COEFFICIENT ##########\n');
[xd, yd, Dkx, Dky, Dkxy, DDx, DDy, DDxy, Dxm, Dym, Dxym, Dkxm, Dkym, Dkxym, DDxm, DDym, DDxym] = evaluateD_2D(input(1:1e7, :), 1/25, r0);
D = (DDxm + DDym) / 2;


%%
%flux based on data
fprintf('########## FLUX ##########\n');
len = 1 / 25;
%len = d;
[xs, ys, pr, flux_x, flux_y, F_data, flux_x1, flux_y1] = flux_2D(input, t0, len);
xf = xs + len / 2;
yf = ys + len / 2;

%%
fprintf('########## KDE MODEL ##########\n');
p = kde(input(1:100:100000001, :)', 'rot' );
s_points = [0.1, 0.1; 0, 1; 1, 0; 1, 1];
minima_array = find_minima_kde(p, s_points, 4, 1e-5);

s1 = (minima_array(1, 1:2) + minima_array(2, 1:2)) / 2;
s1 = s1';
v1 = minima_array(1, 1:2) - minima_array(2, 1:2);
v1 = v1 / sqrt(sum(v1.^2));
v1 = v1';

s2 = (minima_array(1, 1:2) + minima_array(3, 1:2)) / 2;
s2 = s2';
v2 = minima_array(1, 1:2) - minima_array(3, 1:2);
v2 = v2 / sqrt(sum(v2.^2));
v2 = v2';

s3 = (minima_array(2, 1:2) + minima_array(4, 1:2)) / 2;
s3 = s3';
v3 = minima_array(2, 1:2) - minima_array(4, 1:2);
v3 = v3 / sqrt(sum(v3.^2));
v3 = v3';

s4 = (minima_array(3, 1:2) + minima_array(4, 1:2)) / 2;
s4 = s4';
v4 = minima_array(3, 1:2) - minima_array(4, 1:2);
v4 = v4 / sqrt(sum(v4.^2));
v4 = v4';

fprintf('########## SADDLE1 ##########\n');
[saddle1, V1] = saddle_kde_2D(p, s1, v1);

fprintf('########## SADDLE2 ##########\n');
[saddle2, V2] = saddle_kde_2D(p, s2, v2);

fprintf('########## SADDLE3 ##########\n');
[saddle3, V3] = saddle_kde_2D(p, s3, v3);

fprintf('########## SADDLE4 ##########\n');
[saddle4, V4] = saddle_kde_2D(p, s4, v4);


%%
fprintf('########## PATH & SIMULATION MFPT ##########\n');
num = zeros(1, size(input, 1));

num((input(:, 1) - minima_array(1, 1)).^2 + (input(:, 2) - minima_array(1, 2)).^2 <= 1e-4) = 1;
num((input(:, 1) - minima_array(2, 1)).^2 + (input(:, 2) - minima_array(2, 2)).^2 <= 1e-4) = 2;
num((input(:, 1) - minima_array(3, 1)).^2 + (input(:, 2) - minima_array(3, 2)).^2 <= 1e-4) = 3;
num((input(:, 1) - minima_array(4, 1)).^2 + (input(:, 2) - minima_array(4, 2)).^2 <= 1e-4) = 4;

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
way4 = way(way(:, 3) == 4, :);

way1(:, 3) = [];
way2(:, 3) = [];
way3(:, 3) = [];
way4(:, 3) = [];


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
%Create 8 matrices of paths among each 

way1_2 = way1(way1(:, 3) == 2, 1:2);
way1_3 = way1(way1(:, 3) == 3, 1:2);
way2_1 = way2(way2(:, 3) == 1, 1:2);
way2_4 = way2(way2(:, 3) == 4, 1:2);
way3_1 = way3(way3(:, 3) == 1, 1:2);
way3_4 = way3(way3(:, 3) == 4, 1:2);
way4_2 = way4(way4(:, 3) == 2, 1:2);
way4_3 = way4(way4(:, 3) == 3, 1:2);


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
FPT4 = FPT(FPT(:, 3) == 4, :);

FPT1(:, 3) = [];
FPT2(:, 3) = [];
FPT3(:, 3) = [];
FPT4(:, 3) = [];

nFPT1_2 = FPT1(FPT1(:, 3) == 2, 1:2);
nFPT1_3 = FPT1(FPT1(:, 3) == 3, 1:2);
nFPT2_1 = FPT2(FPT2(:, 3) == 1, 1:2);
nFPT2_4 = FPT2(FPT2(:, 3) == 4, 1:2);
nFPT3_1 = FPT3(FPT3(:, 3) == 1, 1:2);
nFPT3_4 = FPT3(FPT3(:, 3) == 4, 1:2);
nFPT4_2 = FPT4(FPT4(:, 3) == 2, 1:2);
nFPT4_3 = FPT4(FPT4(:, 3) == 3, 1:2);

FPT1_2 = zeros(size(nFPT1_2, 1), 1);
FPT1_3 = zeros(size(nFPT1_3, 1), 1);
FPT2_1 = zeros(size(nFPT2_1, 1), 1);
FPT2_4 = zeros(size(nFPT2_4, 1), 1);
FPT3_1 = zeros(size(nFPT3_1, 1), 1);
FPT3_4 = zeros(size(nFPT3_4, 1), 1);
FPT4_2 = zeros(size(nFPT4_2, 1), 1);
FPT4_3 = zeros(size(nFPT4_3, 1), 1);

for i = 1:size(FPT1_2, 1)
    FPT1_2(i) = t0(nFPT1_2(i, 2)) - t0(nFPT1_2(i, 1));
end

for i = 1:size(FPT1_3, 1)
    FPT1_3(i) = t0(nFPT1_3(i, 2)) - t0(nFPT1_3(i, 1));
end

for i = 1:size(FPT2_1, 1)
    FPT2_1(i) = t0(nFPT2_1(i, 2)) - t0(nFPT2_1(i, 1));
end

for i = 1:size(FPT2_4, 1)
    FPT2_4(i) = t0(nFPT2_4(i, 2)) - t0(nFPT2_4(i, 1));
end

for i = 1:size(FPT3_1, 1)
    FPT3_1(i) = t0(nFPT3_1(i, 2)) - t0(nFPT3_1(i, 1));
end

for i = 1:size(FPT3_4, 1)
    FPT3_4(i) = t0(nFPT3_4(i, 2)) - t0(nFPT3_4(i, 1));
end

for i = 1:size(FPT4_2, 1)
    FPT4_2(i) = t0(nFPT4_2(i, 2)) - t0(nFPT4_2(i, 1));
end

for i = 1:size(FPT4_3, 1)
    FPT4_3(i) = t0(nFPT4_3(i, 2)) - t0(nFPT4_3(i, 1));
end


mMFPT1_2 = mean(FPT1_2);
mMFPT1_3 = mean(FPT1_3);
mMFPT2_1 = mean(FPT2_1);
mMFPT2_4 = mean(FPT2_4);
mMFPT3_1 = mean(FPT3_1);
mMFPT3_4 = mean(FPT3_4);
mMFPT4_2 = mean(FPT4_2);
mMFPT4_3 = mean(FPT4_3);


%%
%Divide each path into almost equidistantly into (n+1) points

n = max(2 * max(way(:, 2) - way(:, 1)), 100);

WAY1_2 = zeros(n + 1, 2, size(way1_2, 1));
WAY1_3 = zeros(n + 1, 2, size(way1_3, 1));
WAY2_1 = zeros(n + 1, 2, size(way2_1, 1));
WAY2_4 = zeros(n + 1, 2, size(way2_4, 1));
WAY3_1 = zeros(n + 1, 2, size(way3_1, 1));
WAY3_4 = zeros(n + 1, 2, size(way3_4, 1));
WAY4_2 = zeros(n + 1, 2, size(way4_2, 1));
WAY4_3 = zeros(n + 1, 2, size(way4_3, 1));

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


%Calculate WAY1_3

for i = 1:size(way1_3, 1)
    tmp_num = way1_3(i, 1):way1_3(i, 2);
    len = 1:(way1_3(i, 2) - way1_3(i, 1));
    
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
            WAY1_3(k + 1, 1, i) = (1 - lambda) * input(tmp_num(j), 1) ...
                + lambda * input(tmp_num(j + 1), 1);
            WAY1_3(k + 1, 2, i) = (1 - lambda) * input(tmp_num(j), 2) ...
                + lambda * input(tmp_num(j + 1), 2);
        end
    end
    WAY1_3(n + 1, 1, i) = input(tmp_num(size(tmp_num, 2)), 1);
    WAY1_3(n + 1, 2, i) = input(tmp_num(size(tmp_num, 2)), 2);
    
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


%Calculate WAY2_4

for i = 1:size(way2_4, 1)
    tmp_num = way2_4(i, 1):way2_4(i, 2);
    len = 1:(way2_4(i, 2) - way2_4(i, 1));
    
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
            WAY2_4(k + 1, 1, i) = (1 - lambda) * input(tmp_num(j), 1) ...
                + lambda * input(tmp_num(j + 1), 1);
            WAY2_4(k + 1, 2, i) = (1 - lambda) * input(tmp_num(j), 2) ...
                + lambda * input(tmp_num(j + 1), 2);
        end
    end
    WAY2_4(n + 1, 1, i) = input(tmp_num(size(tmp_num, 2)), 1);
    WAY2_4(n + 1, 2, i) = input(tmp_num(size(tmp_num, 2)), 2);
    
end


%Calculate WAY3_1

for i = 1:size(way3_1, 1)
    tmp_num = way3_1(i, 1):way3_1(i, 2);
    len = 1:(way3_1(i, 2) - way3_1(i, 1));
    
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
            WAY3_1(k + 1, 1, i) = (1 - lambda) * input(tmp_num(j), 1) ...
                + lambda * input(tmp_num(j + 1), 1);
            WAY3_1(k + 1, 2, i) = (1 - lambda) * input(tmp_num(j), 2) ...
                + lambda * input(tmp_num(j + 1), 2);
        end
    end
    WAY3_1(n + 1, 1, i) = input(tmp_num(size(tmp_num, 2)), 1);
    WAY3_1(n + 1, 2, i) = input(tmp_num(size(tmp_num, 2)), 2);
    
end


%Calculate WAY3_4

for i = 1:size(way3_4, 1)
    tmp_num = way3_4(i, 1):way3_4(i, 2);
    len = 1:(way3_4(i, 2) - way3_4(i, 1));
    
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
            WAY3_4(k + 1, 1, i) = (1 - lambda) * input(tmp_num(j), 1) ...
                + lambda * input(tmp_num(j + 1), 1);
            WAY3_4(k + 1, 2, i) = (1 - lambda) * input(tmp_num(j), 2) ...
                + lambda * input(tmp_num(j + 1), 2);
        end
    end
    WAY3_4(n + 1, 1, i) = input(tmp_num(size(tmp_num, 2)), 1);
    WAY3_4(n + 1, 2, i) = input(tmp_num(size(tmp_num, 2)), 2);
    
end


%Calculate WAY4_2

for i = 1:size(way4_2, 1)
    tmp_num = way4_2(i, 1):way4_2(i, 2);
    len = 1:(way4_2(i, 2) - way4_2(i, 1));
    
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
            WAY4_2(k + 1, 1, i) = (1 - lambda) * input(tmp_num(j), 1) ...
                + lambda * input(tmp_num(j + 1), 1);
            WAY4_2(k + 1, 2, i) = (1 - lambda) * input(tmp_num(j), 2) ...
                + lambda * input(tmp_num(j + 1), 2);
        end
    end
    WAY4_2(n + 1, 1, i) = input(tmp_num(size(tmp_num, 2)), 1);
    WAY4_2(n + 1, 2, i) = input(tmp_num(size(tmp_num, 2)), 2);
    
end


%Calculate WAY4_3

for i = 1:size(way4_3, 1)
    tmp_num = way4_3(i, 1):way4_3(i, 2);
    len = 1:(way4_3(i, 2) - way4_3(i, 1));
    
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
            WAY4_3(k + 1, 1, i) = (1 - lambda) * input(tmp_num(j), 1) ...
                + lambda * input(tmp_num(j + 1), 1);
            WAY4_3(k + 1, 2, i) = (1 - lambda) * input(tmp_num(j), 2) ...
                + lambda * input(tmp_num(j + 1), 2);
        end
    end
    WAY4_3(n + 1, 1, i) = input(tmp_num(size(tmp_num, 2)), 1);
    WAY4_3(n + 1, 2, i) = input(tmp_num(size(tmp_num, 2)), 2);
    
end


%%
%Calculate the mean path

mean1_2 = mean(WAY1_2, 3);
mean1_3 = mean(WAY1_3, 3);
mean2_1 = mean(WAY2_1, 3);
mean2_4 = mean(WAY2_4, 3);
mean3_1 = mean(WAY3_1, 3);
mean3_4 = mean(WAY3_4, 3);
mean4_2 = mean(WAY4_2, 3);
mean4_3 = mean(WAY4_3, 3);


%%
%Calculate the saddle point for each mean path
fprintf('########## PATH ENERGY & SADDLE ##########\n');
% gmm_pdf = @(x, y)reshape(pdf(gmm_2D, [x(:) y(:)]), size(x));
% energy = @(x, y) -log(gmm_pdf(x, y));

%probablity density of paths

pd1_2 = zeros(1, size(mean1_2, 1));
pd1_3 = zeros(1, size(mean1_3, 1));
pd2_1 = zeros(1, size(mean2_1, 1));
pd2_4 = zeros(1, size(mean2_4, 1));
pd3_1 = zeros(1, size(mean3_1, 1));
pd3_4 = zeros(1, size(mean3_4, 1));
pd4_2 = zeros(1, size(mean4_2, 1));
pd4_3 = zeros(1, size(mean4_3, 1));

%energy of paths

e1_2 = zeros(1, size(mean1_2, 1));
e1_3 = zeros(1, size(mean1_3, 1));
e2_1 = zeros(1, size(mean2_1, 1));
e2_4 = zeros(1, size(mean2_4, 1));
e3_1 = zeros(1, size(mean3_1, 1));
e3_4 = zeros(1, size(mean3_4, 1));
e4_2 = zeros(1, size(mean4_2, 1));
e4_3 = zeros(1, size(mean4_3, 1));


for i = 1:size(mean1_2, 1)
    pd1_2(i) = evaluate(p, [mean1_2(i, 1); mean1_2(i, 2)]);
    e1_2(i) = -log(pd1_2(i));
end

for i = 1:size(mean1_3, 1)
    pd1_3(i) = evaluate(p, [mean1_3(i, 1); mean1_3(i, 2)]);
    e1_3(i) = -log(pd1_3(i));
end

for i = 1:size(mean2_1, 1)
    pd2_1(i) = evaluate(p, [mean2_1(i, 1); mean2_1(i, 2)]);
    e2_1(i) = -log(pd2_1(i));
end

for i = 1:size(mean2_4, 1)
    pd2_4(i) = evaluate(p, [mean2_4(i, 1); mean2_4(i, 2)]);
    e2_4(i) = -log(pd2_4(i));
end

for i = 1:size(mean3_1, 1)
    pd3_1(i) = evaluate(p, [mean3_1(i, 1); mean3_1(i, 2)]);
    e3_1(i) = -log(pd3_1(i));
end

for i = 1:size(mean3_4, 1)
    pd3_4(i) = evaluate(p, [mean3_4(i, 1); mean3_4(i, 2)]);
    e3_4(i) = -log(pd3_4(i));
end

for i = 1:size(mean4_2, 1)
    pd4_2(i) = evaluate(p, [mean4_2(i, 1); mean4_2(i, 2)]);
    e4_2(i) = -log(pd4_2(i));
end

for i = 1:size(mean4_3, 1)
    pd4_3(i) = evaluate(p, [mean4_3(i, 1); mean4_3(i, 2)]);
    e4_3(i) = -log(pd4_3(i));
end

[saddle_pd1_2, saddle_num1_2] = min(pd1_2);
[saddle_pd1_3, saddle_num1_3] = min(pd1_3);
[saddle_pd2_1, saddle_num2_1] = min(pd2_1);
[saddle_pd2_4, saddle_num2_4] = min(pd2_4);
[saddle_pd3_1, saddle_num3_1] = min(pd3_1);
[saddle_pd3_4, saddle_num3_4] = min(pd3_4);
[saddle_pd4_2, saddle_num4_2] = min(pd4_2);
[saddle_pd4_3, saddle_num4_3] = min(pd4_3);

saddle1_pd = evaluate(p, saddle1);
saddle1_e = -log(saddle1_pd);
saddle2_pd = evaluate(p, saddle2);
saddle2_e = -log(saddle2_pd);
saddle3_pd = evaluate(p, saddle3);
saddle3_e = -log(saddle3_pd);
saddle4_pd = evaluate(p, saddle4);
saddle4_e = -log(saddle4_pd);


%%
fprintf('########## LANDSCAPE ##########\n');
[x, y] = meshgrid(0:0.01:1.5, 0:0.01:1.5);
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
pr2 = plot3(mean1_3(:, 1), mean1_3(:, 2), pd1_3, 'r.-');
hold on;
pr4 = plot3(mean2_1(:, 1), mean2_1(:, 2), pd2_1, 'y.-');
hold on;
pr6 = plot3(mean2_4(:, 1), mean2_4(:, 2), pd2_4, 'r.-');
hold on;
pr7 = plot3(mean3_1(:, 1), mean3_1(:, 2), pd3_1, 'y.-');
hold on;
pr9 = plot3(mean3_4(:, 1), mean3_4(:, 2), pd3_4, 'r.-');
hold on;
pr11 = plot3(mean4_2(:, 1), mean4_2(:, 2), pd4_2, 'y.-');
hold on;
pr12 = plot3(mean4_3(:, 1), mean4_3(:, 2), pd4_3, 'y.-');
hold on;

pr13 = plot3(mean1_2(saddle_num1_2, 1), mean1_2(saddle_num1_2, 2), ...
    pd1_2(saddle_num1_2), 'ro', 'markerfacecolor', 'r');
hold on;
pr14 = plot3(mean1_3(saddle_num1_3, 1), mean1_3(saddle_num1_3, 2), ...
    pd1_3(saddle_num1_3), 'ro', 'markerfacecolor', 'r');
hold on;
pr16 = plot3(mean2_1(saddle_num2_1, 1), mean2_1(saddle_num2_1, 2), ...
    pd2_1(saddle_num2_1), 'yo', 'markerfacecolor', 'y');
hold on;
pr18 = plot3(mean2_4(saddle_num2_4, 1), mean2_4(saddle_num2_4, 2), ...
    pd2_4(saddle_num2_4), 'ro', 'markerfacecolor', 'r');
hold on;
pr19 = plot3(mean3_1(saddle_num3_1, 1), mean3_1(saddle_num3_1, 2), ...
    pd3_1(saddle_num3_1), 'yo', 'markerfacecolor', 'y');
hold on;
pr21 = plot3(mean3_4(saddle_num3_4, 1), mean3_4(saddle_num3_4, 2), ...
    pd3_4(saddle_num3_4), 'ro', 'markerfacecolor', 'r');
hold on;
pr23 = plot3(mean4_2(saddle_num4_2, 1), mean4_2(saddle_num4_2, 2), ...
    pd4_2(saddle_num4_2), 'yo', 'markerfacecolor', 'y');
hold on;
pr24 = plot3(mean4_3(saddle_num4_3, 1), mean4_3(saddle_num4_3, 2), ...
    pd4_3(saddle_num4_3), 'yo', 'markerfacecolor', 'y');
hold on;

pr25 = plot3(saddle1(1), saddle1(2), saddle1_pd, 'go', 'markerfacecolor', 'g');
hold on;
pr26 = plot3(saddle2(1), saddle2(2), saddle2_pd, 'go', 'markerfacecolor', 'g');
hold on;
pr27 = plot3(saddle3(1), saddle3(2), saddle3_pd, 'go', 'markerfacecolor', 'g');
hold on;
pr28 = plot3(saddle4(1), saddle4(2), saddle4_pd, 'go', 'markerfacecolor', 'g');
hold on;

name = {'XLim', 'YLim', 'Xtick', 'Ytick'};
value = {[0 1.5], [0 1.5], (0:0.5:1.5), (0:0.5:1.5)};
set(gca, name, value);
xlabel('Protein a');
ylabel('Protein b');
legend([pr1 pr13 pr4 pr16 pr25], {'small to large', 'saddle path', ...
    'large to small', 'saddle path', 'saddle'}, 'Location', 'north', ...
    'NumColumns', 3);
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
p2 = plot(mean1_3(:, 1), mean1_3(:, 2), 'r.-');
hold on;
p4 = plot(mean2_1(:, 1), mean2_1(:, 2), 'y.-');
hold on;
p6 = plot(mean2_4(:, 1), mean2_4(:, 2), 'r.-');
hold on;
p7 = plot(mean3_1(:, 1), mean3_1(:, 2), 'y.-');
hold on;
p9 = plot(mean3_4(:, 1), mean3_4(:, 2), 'r.-');
hold on;
p11 = plot(mean4_2(:, 1), mean4_2(:, 2), 'y.-');
hold on;
p12 = plot(mean4_3(:, 1), mean4_3(:, 2), 'y.-');
hold on;

p13 = plot(mean1_2(saddle_num1_2, 1), mean1_2(saddle_num1_2, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
p14 = plot(mean1_3(saddle_num1_3, 1), mean1_3(saddle_num1_3, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
p16 = plot(mean2_1(saddle_num2_1, 1), mean2_1(saddle_num2_1, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;
p18 = plot(mean2_4(saddle_num2_4, 1), mean2_4(saddle_num2_4, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
p19 = plot(mean3_1(saddle_num3_1, 1), mean3_1(saddle_num3_1, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;
p21 = plot(mean3_4(saddle_num3_4, 1), mean3_4(saddle_num3_4, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
p23 = plot(mean4_2(saddle_num4_2, 1), mean4_2(saddle_num4_2, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;
p24 = plot(mean4_3(saddle_num4_3, 1), mean4_3(saddle_num4_3, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;

p25 = plot(saddle1(1), saddle1(2), 'go', 'markerfacecolor', 'g');
hold on;
p26 = plot(saddle2(1), saddle2(2), 'go', 'markerfacecolor', 'g');
hold on;
p27 = plot(saddle3(1), saddle3(2), 'go', 'markerfacecolor', 'g');
hold on;
p28 = plot(saddle4(1), saddle4(2), 'go', 'markerfacecolor', 'g');
hold on;

name = {'XLim', 'YLim', 'Xtick', 'Ytick'};
value = {[0 1.5], [0 1.5], (0:0.5:1.5), (0:0.5:1.5)};
set(gca, name, value);
xlabel('Protein a');
ylabel('Protein b');
legend([p1 p13 p4 p16 p25], {'small to large', 'saddle path', ...
    'large to small', 'saddle path', 'saddle'}, 'Location', 'north', ...
    'NumColumns', 3);
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


o2 = plot3(mean1_2(:, 1), mean1_2(:, 2), e1_2, 'r.-');
hold on;
o3 = plot3(mean1_3(:, 1), mean1_3(:, 2), e1_3, 'r.-');
hold on;
o5 = plot3(mean2_1(:, 1), mean2_1(:, 2), e2_1, 'y.-');
hold on;
o7 = plot3(mean2_4(:, 1), mean2_4(:, 2), e2_4, 'r.-');
hold on;
o8 = plot3(mean3_1(:, 1), mean3_1(:, 2), e3_1, 'y.-');
hold on;
o10 = plot3(mean3_4(:, 1), mean3_4(:, 2), e3_4, 'r.-');
hold on;
o12 = plot3(mean4_2(:, 1), mean4_2(:, 2), e4_2, 'y.-');
hold on;
o13 = plot3(mean4_3(:, 1), mean4_3(:, 2), e4_3, 'y.-');
hold on;

o14 = plot3(mean1_2(saddle_num1_2, 1), mean1_2(saddle_num1_2, 2), ...
    e1_2(saddle_num1_2), 'ro', 'markerfacecolor', 'r');
hold on;
o15 = plot3(mean1_3(saddle_num1_3, 1), mean1_3(saddle_num1_3, 2), ...
    e1_3(saddle_num1_3), 'ro', 'markerfacecolor', 'r');
hold on;
o17 = plot3(mean2_1(saddle_num2_1, 1), mean2_1(saddle_num2_1, 2), ...
    e2_1(saddle_num2_1), 'yo', 'markerfacecolor', 'y');
hold on;
o19 = plot3(mean2_4(saddle_num2_4, 1), mean2_4(saddle_num2_4, 2), ...
    e2_4(saddle_num2_4), 'ro', 'markerfacecolor', 'r');
hold on;
o20 = plot3(mean3_1(saddle_num3_1, 1), mean3_1(saddle_num3_1, 2), ...
    e3_1(saddle_num3_1), 'yo', 'markerfacecolor', 'y');
hold on;
o22 = plot3(mean3_4(saddle_num3_4, 1), mean3_4(saddle_num3_4, 2), ...
    e3_4(saddle_num3_4), 'ro', 'markerfacecolor', 'r');
hold on;
o24 = plot3(mean4_2(saddle_num4_2, 1), mean4_2(saddle_num4_2, 2), ...
    e4_2(saddle_num4_2), 'yo', 'markerfacecolor', 'y');
hold on;
o25 = plot3(mean4_3(saddle_num4_3, 1), mean4_3(saddle_num4_3, 2), ...
    e4_3(saddle_num4_3), 'yo', 'markerfacecolor', 'y');
hold on;

o26 = plot3(saddle1(1), saddle1(2), saddle1_e, 'go', 'markerfacecolor', 'g');
hold on;
o27 = plot3(saddle2(1), saddle2(2), saddle2_e, 'go', 'markerfacecolor', 'g');
hold on;
o28 = plot3(saddle3(1), saddle3(2), saddle3_e, 'go', 'markerfacecolor', 'g');
hold on;
o29 = plot3(saddle4(1), saddle4(2), saddle4_e, 'go', 'markerfacecolor', 'g');
hold on;

name = {'XLim', 'YLim', 'Xtick', 'Ytick'};
value = {[0 1.5], [0 1.5], (0:0.5:1.5), (0:0.5:1.5)};
set(gca, name, value);
xlabel('Protein a');
ylabel('Protein b');
legend([o2 o14 o5 o17 o26], {'small to large', 'saddle path', ...
    'large to small', 'saddle path', 'saddle'}, 'Location', 'north', ...
    'NumColumns', 3);
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

l2 = plot(mean1_2(:, 1), mean1_2(:, 2), 'r.-');
hold on;
l3 = plot(mean1_3(:, 1), mean1_3(:, 2), 'r.-');
hold on;
l5 = plot(mean2_1(:, 1), mean2_1(:, 2), 'y.-');
hold on;
l7 = plot(mean2_4(:, 1), mean2_4(:, 2), 'r.-');
hold on;
l8 = plot(mean3_1(:, 1), mean3_1(:, 2), 'y.-');
hold on;
l10 = plot(mean3_4(:, 1), mean3_4(:, 2), 'r.-');
hold on;
l12 = plot(mean4_2(:, 1), mean4_2(:, 2), 'y.-');
hold on;
l13 = plot(mean4_3(:, 1), mean4_3(:, 2), 'y.-');
hold on;

l14 = plot(mean1_2(saddle_num1_2, 1), mean1_2(saddle_num1_2, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
l15 = plot(mean1_3(saddle_num1_3, 1), mean1_3(saddle_num1_3, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
l17 = plot(mean2_1(saddle_num2_1, 1), mean2_1(saddle_num2_1, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;
l19 = plot(mean2_4(saddle_num2_4, 1), mean2_4(saddle_num2_4, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
l20 = plot(mean3_1(saddle_num3_1, 1), mean3_1(saddle_num3_1, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;
l22 = plot(mean3_4(saddle_num3_4, 1), mean3_4(saddle_num3_4, 2), ...
    'ro', 'markerfacecolor', 'r');
hold on;
l24 = plot(mean4_2(saddle_num4_2, 1), mean4_2(saddle_num4_2, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;
l25 = plot(mean4_3(saddle_num4_3, 1), mean4_3(saddle_num4_3, 2), ...
    'yo', 'markerfacecolor', 'y');
hold on;

l26 = plot(saddle1(1), saddle1(2), 'go', 'markerfacecolor', 'g');
hold on;
l27 = plot(saddle2(1), saddle2(2), 'go', 'markerfacecolor', 'g');
hold on;
l28 = plot(saddle3(1), saddle3(2), 'go', 'markerfacecolor', 'g');
hold on;
l29 = plot(saddle4(1), saddle4(2), 'go', 'markerfacecolor', 'g');
hold on;

axis equal;
name = {'XLim', 'YLim', 'Xtick', 'Ytick'};
value = {[0 1.5], [0 1.5], (0:0.5:1.5), (0:0.5:1.5)};
set(gca, name, value);
xlabel('Protein a');
ylabel('Protein b');
legend([l2 l14 l5 l17 l26], {'small to large', 'saddle path', ...
    'large to small', 'saddle path', 'saddle'}, 'Location', 'north', ...
    'NumColumns', 3);
legend('boxoff')


%%
%MFPTnew
saddleCoords = zeros(8, 2);
saddleCoords(1, :) =  mean1_2(saddle_num1_2, :);
saddleCoords(2, :) =  mean1_3(saddle_num1_3, :); 
saddleCoords(3, :) =  mean2_1(saddle_num2_1, :); 
saddleCoords(4, :) =  mean2_4(saddle_num2_4, :);
saddleCoords(5, :) =  mean3_1(saddle_num3_1, :);
saddleCoords(6, :) =  mean3_4(saddle_num3_4, :); 
saddleCoords(7, :) =  mean4_2(saddle_num4_2, :); 
saddleCoords(8, :) =  mean4_3(saddle_num4_3, :);

%saddleEnergy
saddleEnergy = zeros(1, 8);
for i = 1:8
    saddleEnergy(i) = -log(evaluate(p, [saddleCoords(i, 1); saddleCoords(i, 2)]));
end

%saddleMinima
saddleMinima = [1, 2; 1, 3; 2, 1; 2, 4; 3, 1; 3, 4; 4, 2; 4, 3];

minimaHessian_kde = calculateHessian_kde_2D(p, minima_array(:, 1:2));
saddleHessian_kde = calculateHessian_kde_2D(p, saddleCoords);
MFPT2 = mMFPTmatrixWithHessian(minima_array(:, 3),minimaHessian_kde,saddleEnergy',saddleHessian_kde,[saddleMinima(:,1) [1:size(saddleMinima(:,1),1)]' saddleMinima(:,2)],1,1/D);
MFPT4 = mMFPTmatrixWithReciprocal(minima_array(:, 3),minimaHessian_kde,saddleEnergy',saddleHessian_kde,[saddleMinima(:,1) [1:size(saddleMinima(:,1),1)]' saddleMinima(:,2)],1,1/D);


%%
%MFPTNEB
barrierMinima = [1, 2; 1, 3; 2, 4; 3, 4];
barrierHessian_kde = calculateHessian_kde_2D(p, [saddle1'; saddle2'; saddle3'; saddle4']);
MFPT1 = MFPTmatrixWithHessian(minima_array(:, 3),minimaHessian_kde,[saddle1_e; saddle2_e; saddle3_e; saddle4_e],barrierHessian_kde,[barrierMinima(:,1) [1:size(barrierMinima(:,1),1)]' barrierMinima(:,2)],1,1/D);
MFPT3 = MFPTmatrixWithReciprocal(minima_array(:, 3),minimaHessian_kde,[saddle1_e; saddle2_e; saddle3_e; saddle4_e],barrierHessian_kde,[barrierMinima(:,1) [1:size(barrierMinima(:,1),1)]' barrierMinima(:,2)],1,1/D);


%%
MFPT15_kde = {MFPT2,MFPT4,mMFPT1_2,mMFPT1_3,mMFPT2_1,mMFPT2_4,mMFPT3_1,mMFPT3_4,mMFPT4_2,mMFPT4_3,MFPT1,MFPT3};
D15 = D;
1/D15
sound(sin(2*pi*25*(1:4000)/100));