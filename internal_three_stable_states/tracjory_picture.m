
randn('seed',2020)

fprintf('########## SIMULATION ##########\n');
a = 1; b = 1; k0 = 1; S0 = 0.5; l = 4; x0 = 0.5; y0 = 0.5; tmin = 0; 
tmax = 1e5; h = 0.1; V = 35; r = 4;
input = euler_simD_r(a, b, k0, S0, l, x0, y0, tmin, tmax, h, V, r);

r0 = h * r;
t = tmin:r0:tmax;
t0 = t;
figure(1);
se1 = plot(t, input(:, 1), 'k.-');
% set(se1, 'LineStyle','none');
xlabel('Time');
ylabel('Expression');

figure(2);
se2 = plot(t, input(:, 2), 'k.-');
% set(se2, 'LineStyle','none');
xlabel('Time');
ylabel('Expression');

figure(3);
se3 = plot(t(1:1000), input(1:1000, 1), 'k.-');
% set(se3, 'LineStyle','none');
xlabel('Time');
ylabel('Expression');

figure(4);
se4 = plot(t(1:1000), input(1:1000, 2), 'k.-');
% set(se4, 'LineStyle','none');
xlabel('Time');
ylabel('Expression');

% 
% 
% %%
% %landscape new_2D
% ez1 = ez;
% for i = 1:size(ez, 1)
%     for j = 1:size(ez, 2)
%         if ez(i, j) >= 7
%             ez1(i, j) = 7;
%         end
%     end
% end
% figure(5);
% se = pcolor(x, y, ez1);
% set(se, 'LineStyle','none');
% colorbar;
% hold on;
% cs1= 1:2:size(xf,1);
% quiver(xf(cs1,cs1), yf(cs1,cs1), flux_x(cs1,cs1), flux_y(cs1,cs1), 'w', 'linewidth', 0.8);
% hold on;
% 
% pe2 = plot(mean1_2(:, 1), mean1_2(:, 2), 'r.-');
% hold on;
% pe3 = plot(mean2_1(:, 1), mean2_1(:, 2), 'y.-');
% hold on;
% pe4 = plot(mean2_3(:, 1), mean2_3(:, 2), 'r.-');
% hold on;
% pe5 = plot(mean3_2(:, 1), mean3_2(:, 2), 'y.-');
% hold on;
% 
% pe6 = plot(mean1_2(saddle_num1_2, 1), mean1_2(saddle_num1_2, 2), ...
%     'ro', 'markerfacecolor', 'r');
% hold on;
% pe7 = plot(mean2_1(saddle_num2_1, 1), mean2_1(saddle_num2_1, 2), ...
%     'yo', 'markerfacecolor', 'y');
% hold on;
% pe8 = plot(mean2_3(saddle_num2_3, 1), mean2_3(saddle_num2_3, 2), ...
%     'ro', 'markerfacecolor', 'r');
% hold on;
% pe9 = plot(mean3_2(saddle_num3_2, 1), mean3_2(saddle_num3_2, 2), ...
%     'yo', 'markerfacecolor', 'y');
% hold on;
% 
% pe10 = plot(saddle1(1), saddle1(2), 'go', 'markerfacecolor', 'g');
% hold on;
% pe11 = plot(saddle2(1), saddle2(2), 'go', 'markerfacecolor', 'g');
% hold on;
% 
% axis equal;
% name = {'XLim', 'YLim', 'Xtick', 'Ytick'};
% value = {[0 2.5], [0 2.5], (0:0.5:2.5), (0:0.5:2.5)};
% set(gca, name, value);
% xlabel('Protein a');
% ylabel('Protein b');
