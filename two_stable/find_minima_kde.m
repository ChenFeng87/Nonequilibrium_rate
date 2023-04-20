function minima_array = find_minima_kde(p, s_points, max_jumps, max_rms, similarity, iter_max, dx)
%p is kde, s_points is initial guess(M*N, M is number of points)

minima_array = [];
if nargin <= 1
    s_points = [];
end
if nargin <= 2
    max_jumps = 3;
end
if nargin <= 3
    max_rms = 1e-6; 
end
if nargin <= 4
    similarity = 1e-4;
end
if nargin <= 5
    iter_max = 1000;
end
if nargin <= 6
    dx = 5e-3;  % ZiMo:升高了5倍
end

N = getDim(p);
data = getPoints(p);
boundary = zeros(N, 2);
for i = 1:N
    boundary(i, 1) = min(data(i, :));
    boundary(i, 2) = max(data(i, :));
end

fprintf('########## MINIMA ##########\n');

for jump = 1:max_jumps
    iter = 1;
    converged = 0;
    x = zeros(1, N);
    for i = 1:N
        x(i) = (boundary(i, 2) - boundary(i, 1)) * rand(1) + boundary(i, 1);
    end
    
    if jump <= size(s_points, 1)
        x = s_points(jump, :);
    end
    
    while ((iter < iter_max) && (converged == 0))
        pr = evaluate(p, x');
        g = zeros(1, N);
        for i = 1:N
            x1 = zeros(1, N);
            x1(i) = 1e-6;
            g(i) = (evaluate(p, (x1 + x)') - pr) / 1e-6;
        end
        g = -g / pr;
        x = x - dx * g;
        rms = sum(g.^2);
        
        if rms < max_rms
            converged = 1;
            e = -log(evaluate(p, x'));
        end
        
        if rms ~= rms
            iter = iter_max;
            converged = 0;
        end
        
        for i = 1:N
            if ((x(i) > boundary(i, 2)) || (x(i) < boundary(i, 1)))
                iter = iter_max;
                converged = 0;
            end
            if x(i) ~= x(i)
                iter = iter_max;
                converged = 0;
            end
        end
        iter = iter + 1;
    end
    
    if converged == 1
        if size(minima_array, 1) == 0
            minima_array = [x, e];
            fprintf('Located minima 1\n');
        else
            accept = 1;
            for i = 1:size(minima_array, 1)
                if sum((x - minima_array(i, 1:N)).^2) < similarity
                    accept = 0;
                end
            end
            
            if accept == 1
                minima_array = [minima_array; [x, e]];
                w = size(minima_array, 1);
                fprintf('Located minima %d\n', w);
            end
        end
    end
end
fprintf('tebie: Located minima: \n');
iter

end