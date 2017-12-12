function [U,obj]= MultiviewClustering(Data,para)

%  input  Data  cell nview*1      n*dim  
%  input  Lapmatrix cell nview*1  n*n
iter = para.iter;
R = para.R;
alpha = para.alpha; 
Lapmatrix = para.L;  
c = para.c;
mu = para.mu;
nView = length(Data);
%  initialize  theta
beta = zeros(nView,1);
beta(:) = 1/(nView);
% initialize L
L = GlobalLapmatrix(Lapmatrix,beta,R);

%%% initialize U

X = [];    
for i=1:1:nView
     X = [X,Data{i}];      
end
X = X';

[idx_U, A, evs] = CAN(X, c);
for i=1:length(idx_U)
    U(i,idx_U(i))=1;
end
U=U+0.2; 

V = cell(nView,1);
for i=1:1:nView
    
    V{i} = zeros(size(Data{i}',1),c);
end
%  initialize lambda

lambda = cell(nView,1);

for i=1:1:nView
    
    lambda{i} = zeros(size(Data{i}'));
end
    lambda_1 = zeros(size(U));
    
pho = 10;
mu_max = 10^10;  
obj = zeros(iter,1);

E = cell(nView,1);

for i=1:1:iter
    
    % update E
    for j=1:1:nView
        
        temp_e = Data{j}'- V{j}*U'+1/mu*lambda{j};
        [E{j}]=L21_solver(temp_e,1/mu);
        E{j} = real(E{j});
    end
    
    % update V
       for j=1:1:nView
           
           V{j} = (Data{j}'-E{j} + 1/mu*lambda{j})*U;
           V{j} = real(V{j});
       end
    
      % updata Z_1
      
       K = (U - 1/mu*lambda_1 - alpha/mu*L*U);
       [Z_1]=nonneg_L2(K);
       Z_1 = real(Z_1);
      % updata U
      
       H = 1/mu*lambda_1 + Z_1-(alpha/mu)*L*Z_1; 
       for j=1:1:nView
          H = H + (Data{j}'-E{j}+1/mu*lambda{j})'*V{j};
       end
       [N_u,D_u,Q_u] = svds(H,c);
       U = N_u*Q_u';
       U = real(U);
      % update beta
   for j=1:1:nView
         pv = trace(Z_1'*Lapmatrix{j}*U);
         beta(j) = (R*pv)^(1/(R-1));  
   end
     beta = beta/sum(beta,1);
     L = GlobalLapmatrix(Lapmatrix,beta,R); 
     
    % update parameter
   for j=1:1:nView
     lambda{j} = lambda{j} + mu*(Data{j}'-V{j}*U'-E{j});
   end
   lambda_1 = lambda_1 + mu*(Z_1 -U);
   
    mu = max(pho*mu,mu_max);
   % mu = pho*mu;
    temp_obj  = 0;
    for j=1:1:nView
        temp_obj  = temp_obj  + norm(Data{j}' - V{j}*U','fro')^2;
       
    end
         temp_obj  = temp_obj  + alpha*(trace(U'*L*U));
     if i==1
        obj(i,1) = temp_obj;
        obj_last = temp_obj;
    else
        obj(i,1)= abs(obj_last - temp_obj);
        obj_last = temp_obj;
     end 
%     fprintf('******iter= %d,obj = %f******,\n',i,obj(i,1));    
end

