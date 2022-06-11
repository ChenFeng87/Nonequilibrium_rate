function [xd, yd, pr, flux_x, flux_y, F_data, flux_x1, flux_y1] = flux_2D(input, t, a)
%input is a n*2 matrix, every row is a point in 2D plane
%t is time corresponding to the point
%a is the side length of small square

n = size(t, 2);
mi = floor(min(min(input)));
ma = ceil(max(max(input)));
d = floor((ma - mi) / a) + 1;
pr = zeros(d);
flux_x = zeros(d);
flux_y = zeros(d);
[xd, yd] = meshgrid(mi:a:ma, mi:a:ma);

kx = floor((input(1, 1) - mi) / a) + 1;
ky = floor((input(1, 2) - mi) / a) + 1;
pr(ky, kx) = pr(ky, kx) + 1;
for i = 2:n
    kx0 = kx;
    ky0 = ky;
    kx = floor((input(i, 1) - mi) / a) + 1;
    ky = floor((input(i, 2) - mi) / a) + 1;
    pr(ky, kx) = pr(ky, kx) + 1;
    
    qx = kx - kx0;
    qy = ky - ky0;
    q = abs(qx) + abs(qy);
    if q == 1
        flux_x(ky0, kx0) = flux_x(ky0, kx0) + 0.5 * sign(qx);
        flux_y(ky0, kx0) = flux_y(ky0, kx0) + 0.5 * sign(qy);
        flux_x(ky, kx) = flux_x(ky, kx) + 0.5 * sign(qx);
        flux_y(ky, kx) = flux_y(ky, kx) + 0.5 * sign(qy);
    elseif q > 1
        if qx == 0
            flux_y(ky0, kx0) = flux_y(ky0, kx0) + 0.5 * sign(qy);
            flux_y(ky, kx) = flux_y(ky, kx) + 0.5 * sign(qy);
            for j = 1:(abs(qy) - 1)
                flux_y(ky0 + j * sign(qy), kx0) = flux_y(ky0 + j * sign(qy), kx0) + sign(qy);
            end
        elseif qy == 0
            flux_x(ky0, kx0) = flux_x(ky0, kx0) + 0.5 * sign(qx);
            flux_x(ky, kx) = flux_x(ky, kx) + 0.5 * sign(qx);
            for j = 1:(abs(qx) - 1)
                flux_x(ky0, kx0 + j * sign(qx)) = flux_x(ky0, kx0 + j * sign(qx)) + sign(qx);
            end
        else
            lambda_x = zeros(abs(qx), 2);
            lambda_y = ones(abs(qy), 2);
            kx_min = min(kx, kx0);
            ky_min = min(ky, ky0);
            for j = 1:abs(qx)
                lambda_x(j, 1) = (mi + a * (kx_min + j - 1) - input(i - 1, 1)) / (input(i, 1) - input(i - 1, 1));
            end
            for j = 1:abs(qy)
                lambda_y(j, 1) = (mi + a * (ky_min + j - 1) - input(i - 1, 2)) / (input(i, 2) - input(i - 1, 2));
            end
            lambda = [lambda_x; lambda_y];
            lambda = sortrows(lambda, 1);
            kx1 = kx0;
            ky1 = ky0;
            for j = 1:size(lambda, 1)
                if lambda(j, 2) == 0
                    flux_x(kx1, ky1) = flux_x(kx1, ky1) + 0.5 * sign(qx);
                    flux_x(kx1 + sign(qx), ky1) = flux_x(kx1 + sign(qx), ky1) + 0.5 * sign(qx);
                    kx1 = kx1 + sign(qx);
                else
                    flux_y(kx1, ky1) = flux_y(kx1, ky1) + 0.5 * sign(qy);
                    flux_y(kx1, ky1 + sign(qy)) = flux_y(kx1, ky1 + sign(qy)) + 0.5 * sign(qy);
                    ky1 = ky1 + sign(qy);
                end
            end
        end
    end
end

pr = pr / n;
flux_x = flux_x / t(n);
flux_y = flux_y / t(n);

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
flux_x = flux_x(w:(d - k + 1), w:(d - k + 1));
flux_y = flux_y(w:(d - k + 1), w:(d - k + 1));

F_data = zeros(size(flux_x));
for i = 1:size(F_data, 1)
    for j = 1:size(F_data, 2)
        F_data(i, j) = log(sqrt(flux_x(i, j)^2 + flux_y(i, j)^2) + 1e-16);
    end
end

flux_x1 = zeros(size(flux_x));
flux_y1 = zeros(size(flux_y));
for i = 1:size(flux_x1, 1)
    for j = 1:size(flux_x1, 2)
        if abs(flux_x(i, j)) + abs(flux_y(i, j)) ~=0
            flux_x1(i, j) = flux_x(i, j) / (sqrt(flux_x(i, j)^2 + flux_y(i, j)^2));
            flux_y1(i, j) = flux_y(i, j) / (sqrt(flux_x(i, j)^2 + flux_y(i, j)^2));
        end
    end
end

end