randn('seed',2029)
fprintf('########## SIMULATION ##########\n');
a = 1; b = 1; k0 = 1; S0 = 0.5; l = 4; x0 = 0.5; y0 = 0.5; tmin = 0; 
%tmax = 1e6; h = 0.01; V = 30; r = 1;
%tmax = 0.05e7; h = 0.005; V = 80; r = 1;
%tmax = 1e6; h = 0.01; V = 40; r = 1;
%tmax = 1e7; h = 0.1; V = 40; r = 1;
tmax = 4e5; h = 0.1; V = 20; r = 4;
input = euler_simD_r(a, b, k0, S0, l, x0, y0, tmin, tmax, h, V, r);
r0 = h * r;
t = tmin:r0:tmax;
t0 = t;
plot(t, input(:,1));
