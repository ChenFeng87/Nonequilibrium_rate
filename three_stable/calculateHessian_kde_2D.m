function Hessian = calculateHessian_kde_2D(p, points)
D = getDim(p);

if D ~= 2
    error('Dimension is not equal to 2');
end

num_points = size(points, 1);
Hessian = zeros(num_points, D, D);
h = 1e-2;

for i =1:num_points
    x = points(i, :);
    u00 = -log(evaluate(p, x'));
    x1 = [h, 0];
    x2 = [0, h];
    u10 = -log(evaluate(p, (x + x1)'));
    u_10 = -log(evaluate(p, (x - x1)'));
    u01 = -log(evaluate(p, (x + x2)'));
    u0_1 = -log(evaluate(p, (x - x2)'));
    u11 = -log(evaluate(p, (x + x1 + x2)'));
    u1_1 = -log(evaluate(p, (x + x1 - x2)'));
    u_11 = -log(evaluate(p, (x - x1 + x2)'));
    u_1_1 = -log(evaluate(p, (x - x1 - x2)'));
    
    uxx = (u10 + u_10 - 2 * u00) / h^2;
    uyy = (u01 + u0_1 - 2 * u00) / h^2;
    uxy = (u11 + u_1_1 - u1_1 - u_11) / (4 * h^2);
    Hessian(i,:,:) = [uxx, uxy; uxy, uyy];
end
end

