clear all
clc
% addpath('D:\projects\baiyubo_matlab\code_figure')
load('data_52.mat');
rand('seed',1998) % the random seed
fprintf('########## GMM MODEL ##########\n');
k = 2;
options = statset('MaxIter', 2000);
S.mu = [2, 0; 1, 1];
S.Sigma = zeros(2, 2, k);
for i = 1:k
    S.Sigma(:, :, i) =  diag(ones(1, 2));
end
S.ComponentProportion = ones(1, k) / k;
gmm_new = fitgmdist(input(1:10:40000001, :), k, 'Options', options, 'Start', S);
gmm_pdf = @(x, y)reshape(pdf(gmm_new, [x(:) y(:)]), size(x));
sum1=0;
for i=1:20
    t=clock;
    for j=1:1000
        gmm_pdf(1.5+rand(1), 0.5*rand(1));
    end
    sum1=sum1+etime(clock,t);
    sum1/i
end







% fprintf('########## KDE MODEL ##########\n');
% p = kde(input(1:10:40000001, :)', 'rot' ); % 40000001
% sum1=0;
% for i=1:20
%     t=clock;
%     for j=1:1000
%         evaluate(p, [1.5+rand(1); 0.5*rand(1)]);
%     end
%     sum1=sum1+etime(clock,t);
%     sum1/i
% end




% sum1=0;
% for i=1:20
%     t=clock;
%     p = kde(input(1:10:40000001, :)', 'rot' ); % 40000001
%     sum1=sum1+etime(clock,t);
%     sum1/i
%     clear p
% end
% 
% sum2 = 0;
% for j =1:20
%     fprintf('########## GMM MODEL ##########\n');
%     t=clock;
%     k = 2;
%     options = statset('MaxIter', 2000);
%     S.mu = [2, 0; 1, 1];
%     S.Sigma = zeros(2, 2, k);
%     for i = 1:k
%         S.Sigma(:, :, i) =  diag(ones(1, 2));
%     end
%     S.ComponentProportion = ones(1, k) / k;
%     gmm_new = fitgmdist(input(1:10:40000001, :), k, 'Options', options, 'Start', S);
%     sum2=sum2+etime(clock,t);
%     sum2/j
%     clear gmm_new
% end


