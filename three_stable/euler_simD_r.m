function x = euler_simD_r(a, b, k, S, l, x0, y0, tmin, tmax, h, V, r)


r0 = h * r;
t=tmin:r0:tmax;         %%the time points need to be calculated 
n=length(t);
x=zeros(n,2);   %%the process

x(1,1)=x0;
x(1,2)=y0;%%initial points

for i=2:n
    tmp1=x(i - 1, 1);
    tmp2=x(i - 1, 2);
    for j = 1:r
        n1 = randn(1); n2 = randn(1);
        d1 = (sqrt(0.1 + a * tmp1^l / (S^l + tmp1^l) + b * S^l / (S^l + tmp2^l)+k * tmp1)) * n1/ sqrt(V);
        d2 = (sqrt(0.1 + a * tmp2^l / (S^l + tmp2^l) + b * S^l / (S^l + tmp1^l)+k * tmp2)) * n2/ sqrt(V);
        tmp1 = tmp1 + h * (0.1 + a * tmp1^l / (S^l + tmp1^l) + b * S^l / (S^l + tmp2^l) - k * tmp1)  + d1 * sqrt(h);
        tmp2 = tmp2 + h * (0.1 + a * tmp2^l / (S^l + tmp2^l) + b * S^l / (S^l + tmp1^l) - k * tmp2)  + d2 * sqrt(h);
        tmp1 = max(tmp1, 0);
        tmp2 = max(tmp2, 0);
    end
    x(i, 1) = tmp1;
    x(i, 2) = tmp2;
end

end