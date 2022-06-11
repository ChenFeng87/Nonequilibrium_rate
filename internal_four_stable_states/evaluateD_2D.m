function [xd, yd, Dkx, Dky, Dkxy, DDx, DDy, DDxy, Dxm, Dym, Dxym, Dkxm, Dkym, Dkxym, DDxm, DDym, DDxym] = evaluateD_2D(input, a, dt)
%input is a n*2 matrix, every row is a point in 2D plane
%t is time corresponding to the point
%a is the side length of small square

n = size(input, 1);
mi = floor(min(min(input)));
ma = ceil(max(max(input)));
d = floor((ma - mi) / a) + 1;
pr = zeros(d);
xv2 = cell(d);
yv2 = cell(d);
xyv = cell(d);

Dx = cell(d);
Dy = cell(d);
Dxy = cell(d);

Dkx = zeros(d);
Dky = zeros(d);
Dkxy = zeros(d);

xv2m = zeros(1, 10);
yv2m = zeros(1, 10);
xyvm = zeros(1, 10);

for i = 1:d
    for j = 1:d
        Dx{i, j} = zeros(1, 10);
        Dy{i, j} = zeros(1, 10);
        Dxy{i, j} = zeros(1, 10);
        xv2{i, j} = zeros(1, 10);
        yv2{i, j} = zeros(1, 10);
        xyv{i, j} = zeros(1, 10);
    end
end

[xd, yd] = meshgrid(mi:a:ma, mi:a:ma);

for i = 1:(n - 10)
    kx = floor((input(i, 1) - mi) / a) + 1;
    ky = floor((input(i, 2) - mi) / a) + 1;
    pr(ky, kx) = pr(ky, kx) + 1;  % ZiMo: 统计每个点在哪个方格里面
    for j = 1:10
        xv2{ky, kx}(j) = xv2{ky, kx}(j) + (input(i + j, 1) - input(i, 1))^2;  % ZiMo: 计算后面十个位置与该位置的距离和
        yv2{ky, kx}(j) = yv2{ky, kx}(j) + (input(i + j, 2) - input(i, 2))^2;
        xyv{ky, kx}(j) = xyv{ky, kx}(j) + (input(i + j, 1) - input(i, 1)) * (input(i + j, 2) - input(i, 2));
        
        xv2m(j) = xv2m(j) + (input(i + j, 1) - input(i, 1))^2;  % ZiMo: 计算后面十个位置与该位置的距离和
        yv2m(j) = yv2m(j) + (input(i + j, 2) - input(i, 2))^2;
        xyvm(j) = xyvm(j) + (input(i + j, 1) - input(i, 1)) * (input(i + j, 2) - input(i, 2));
    end
end

for u = 1:d
    if sum(pr(u, :) == 0) ~= d
        break
    end
end
for v = 1:d
    if sum(pr(:, v) == 0) ~= d
        break
    end
end
for i = 1:d
    if sum(pr(d - i + 1, :) == 0) ~= d
        break
    end
end
for j = 1:d
    if sum(pr(:, d - j + 1) == 0) ~= d
        break
    end
end

w = min(u, v);
k = min(i, j);
xd = xd(w:(d - k + 1), w:(d - k + 1));
yd = yd(w:(d - k + 1), w:(d - k + 1));
pr = pr(w:(d - k + 1), w:(d - k + 1));
xv2 = xv2(w:(d - k + 1), w:(d - k + 1));
yv2 = yv2(w:(d - k + 1), w:(d - k + 1));
xyv = xyv(w:(d - k + 1), w:(d - k + 1));
Dx = Dx(w:(d - k + 1), w:(d - k + 1));
Dy = Dy(w:(d - k + 1), w:(d - k + 1));
Dxy = Dxy(w:(d - k + 1), w:(d - k + 1));
Dkx = Dkx(w:(d - k + 1), w:(d - k + 1));
Dky = Dky(w:(d - k + 1), w:(d - k + 1));
Dkxy = Dkxy(w:(d - k + 1), w:(d - k + 1));

DDx = zeros(size(pr));
DDy = zeros(size(pr));
DDxy = zeros(size(pr));

for i = 1:size(pr, 1)
    for j = 1:size(pr, 2)
        if pr(i, j) > 0
            Dx{i, j} = xv2{i, j} / pr(i, j);
            Dy{i, j} = yv2{i, j} / pr(i, j);
            Dxy{i, j} = xyv{i, j} / pr(i, j);
            
            Dkx(i, j) = (2 * (1:10) * Dx{i,j}' - 11 * sum(Dx{i,j})) / (330 * dt);
            Dky(i, j) = (2 * (1:10) * Dy{i,j}' - 11 * sum(Dy{i,j})) / (330 * dt);
            Dkxy(i, j) = (2 * (1:10) * Dxy{i,j}' - 11 * sum(Dxy{i,j})) / (330 * dt);
            
            DDx(i, j) = Dx{i, j}(1) / (2 * dt);
            DDy(i, j) = Dy{i, j}(1) / (2 * dt);
            DDxy(i, j) = Dxy{i, j}(1) / (2 * dt);
            
        end
    end
end

Dxm = xv2m / (n - 10);
Dym = yv2m / (n - 10);
Dxym = xyvm / (n - 10);

Dkxm = (2 * (1:10) * Dxm' - 11 * sum(Dxm)) / (330 * dt);
Dkym = (2 * (1:10) * Dym' - 11 * sum(Dym)) / (330 * dt);
Dkxym = (2 * (1:10) * Dxym' - 11 * sum(Dxym)) / (330 * dt);

DDxm = Dxm(1) / (2 * dt);
DDym = Dym(1) / (2 * dt);
DDxym = Dxym(1) / (2 * dt);

end