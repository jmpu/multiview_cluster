clear all;
clc;
addpath(genpath('utils'))
addpath(genpath('Wiki_textimage'))
load('Wiki_textimage.mat');
load('Wikilabel.mat');



count=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% V1,V2,V3分别是三个视图
% V1=mfeatfou;
% V2=mfeatpix;
V1=NormalizeFea(image,1);
V2=NormalizeFea(text,1);

nFea = 2;
Data = cell(nFea,1);
Data{1,1} = V1;
Data{2,1} = V2;
% Data{3,1} = V3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Lapmatrix = cell(nFea,1);
options = [];
options.NeighborMode = 'KNN';
% options.gnd =  [zeros(200,1);ones(200,1);ones(200,1)*2;ones(200,1)*3;ones(200,1)*4;...
%       ones(200,1)*5;ones(200,1)*6;ones(200,1)*7;ones(200,1)*8;ones(200,1)*9];
options.WeightMode = 'HeatKernel';
options.t = 1;
options.k = 5;

for ifeat =1:1:nFea
  W = constructW(Data{ifeat},options); 
  W = full(W);
  D = diag(sum(W,2));
  L = D - W;
  L = D^(-1/2)*L*D^(-1/2);
  Lapmatrix{ifeat} = L;  
end
 para_candi = [10^(-3),10^(-2),0.1,1,10,10^2,10^4];     
 para_candi1= [10^(-6),10^(-3),10^(-2),1,10,100,10^4];
 para_candi2=[5,10,1000];
 para_candi3=[5];
  for i = 1:length(para_candi)
    for   j = 1:length(para_candi1)
        for   k = 1:length(para_candi2)
            for   m = 1:length(para_candi3)
para =[];
para.iter = 100;
para.alpha = para_candi(i);
para.mu = para_candi1(j);
para.R = para_candi2(k);
para.k = para_candi3(m);
para.L = Lapmatrix;
para.c = 10;


[U,obj] = MultiviewClustering(Data,para);
[~, label] = max(U,[],2);
gnd = Z;
   
res = bestMap(gnd,label);
AC = length(find(gnd == res))/length(gnd);
MIhat = MutualInfo(gnd,res);
  Result(i,j,k,m,1) = para.alpha;
      Result(i,j,k,m,2) = para.mu;
      Result(i,j,k,m,3) = para.R;
     Result(i,j,k,m,4) = para.k;
     Result(i,j,k,m,5) =AC;
     Result(i,j,k,m,6) =MIhat;
     count=count+1;
    fprintf('******index= %d/ %d, AC = %f*****\n',count, length(para_candi)*length(para_candi1)*length(para_candi2)*length(para_candi3),Result(i,j,k,m,5));
 end
 end
    end
  end
