function Y = f_MSE(X1 ,X2 ,dim,r,k,sigma)
% X1 ,X2 分别是 V1 和 V2 的特征
% dim 是降维的维数
% para_candi = [1,10,10^2,10^3,10^4];
% para_candi1 = [2:1:9];
% for i = 1:length(para_candi)
%     for j=1:length(para_candi1)
% sigma= para_candi(i);
%    r = para_candi1(i);
% k=12;
L1 = laplacian_eigen_L(X1, k, sigma);
L2 = laplacian_eigen_L(X2,  k, sigma);
[Y, a_rocord] = AlterOpti(L1, L2, dim, r);

%     end
% end
