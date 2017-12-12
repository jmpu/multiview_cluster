function [Y, a_rocord] = AlterOpti(L1, L2, dim, r)

%迭代优化权系数a with UX
a_prior=[0,1];
a=[0.5,0.5];
a_rocord=a;

i=0;
while ((a_prior-a)*(a_prior-a)'>0.00001&&i<5)
    i=i+1;
    a_prior=a;
    
    L=a(1)*L1+a(2)*L2;
    L=(L+L')/2;
    [Y, eig_value] = eigs(L,dim,'sm');
    
    % % normalization
    % Y = Y*diag(1./((sum(Y.^2).^0.5)));
    
    tr1=trace(Y'*L1*Y)+10^(-9);
    tr2=trace(Y'*L2*Y)+10^(-9);
    
    a(1)=(1/tr1)^(1/(r-1));
    a(2)=(1/tr2)^(1/(r-1));
    a_sum=sum(a);
    a(1)=a(1)/a_sum;
    a(2)=a(2)/a_sum;
    a;
    a_rocord=[a_rocord;a];  
end
