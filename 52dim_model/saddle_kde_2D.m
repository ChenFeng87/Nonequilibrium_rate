function [x, v, reccur_n] = saddle_kde_2D(p, x0, v0)
D = getDim(p);

if D ~= 2
    error('Dimension is not equal to 2');
end

I = diag(ones(1, D));

rms = 1e-5;  % 收敛条件加大10倍
dt = 0.01;
h = 1e-4;
n = 5e4; % 1e5
reccur_n=0;

beta = 1e-2;
gamma = 1e-2;

x = x0;
v = v0;

for i = 1:n
    x_ahead = x;
    v_ahead = v;
    
    g = zeros(D, 1);
    g1 = zeros(D, 1);
    g2 = zeros(D, 1);
    for j = 1:D
        x1 = zeros(D, 1);
        x1(j) = h;
        g(j)=(-log(evaluate(p,x_ahead+x1))+log(evaluate(p,x_ahead-x1)))/(2*h);
        g1(j)=(-log(evaluate(p,x_ahead+h*v_ahead+x1))+log(evaluate(p,x_ahead+h*v_ahead-x1)))/(2*h);
        g2(j)=(-log(evaluate(p,x_ahead-h*v_ahead+x1))+log(evaluate(p,x_ahead-h*v_ahead-x1)))/(2*h);
    end
    hxv = (g1 - g2) / (2 * h);
    
    xg = (I - 2 * (v_ahead * v_ahead')) * g;
    vg = (I - (v_ahead * v_ahead')) * hxv;
    
    x = x_ahead - dt * beta * xg;
    v = v_ahead - dt * gamma * vg;
    v = v / sqrt(sum(v.^2));
    
    x_rms = sum(xg.^2);
    v_rms = sum(vg.^2);
    m_rms = max(x_rms, v_rms);
    if x_rms <= rms
        fprintf('######### CONVERGED #########\n');
        fprintf('iter = %d\n', i);
        break
    end
end

if i == n
    fprintf('######### MAXIMUM ITERATION REACHED #########\n');
    fprintf('rms = %f\n', m_rms);
    reccur_n=1;
end


end