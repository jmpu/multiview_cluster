function L = laplacian_eigen_L(X, k, sigma)
%LAPLACIAN_EIGEN Performs non-linear dimensionality reduction using Laplacian Eigenmaps
%
%   [mappedX, mapping] = laplacian_eigen(X, no_dims, k, sigma, eig_impl)
%
% Performs non-linear dimensionality reduction using Laplacian Eigenmaps.
% The data is in matrix X, in which the rows are the observations and the
% columns the dimensions. The variable dim indicates the preferred amount
% of dimensions to retain (default = 2). The variable k is the number of 
% neighbours in the graph (default = 12).
% The reduced data is returned in the matrix mappedX.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology


    if ~exist('k', 'var')
        k = 12;
    end
	if ~exist('sigma', 'var')
		sigma = 1;
    end
    
    % Construct neighborhood graph
    % disp('Constructing neighborhood graph...');
    
    if size(X, 1) < 4000
        G = L2_distance(X', X');
        
        % Compute neighbourhood graph
        [~, ind] = sort(G, 'ascend'); 
        for i=1:size(G, 1)
            G(i, ind((2 + k):end, i)) = 0; 
        end
        G = sparse(double(G));
        %G = max(G, G');             % Make sure distance matrix is symmetric
        %G = (G+ G')/2;
    else
        G = find_nn(X, k);
        % G = max(G, G');  
        %G = (G+ G')/2;
    end
    
%     Gl = L2_distance(X', X');
%     % Compute neighbourhood graph
%     [~, ind] = sort(Gl, 'ascend'); 
%     for i=1:size(Gl, 1)
%         Gl(i, ind((2 + k):end, i)) = 0; 
%     end
%     Gl = sparse(double(Gl));
%     %Gl = max(Gl, Gl');             % Make sure distance matrix is symmetric
%     Gf = find_nn(X, k);

    G = G .^ 2;
    max_dist = max(max(G));
    G = G ./ max_dist;
    
    % Compute weights (W = G)
    % disp('Computing weight matrices...');
    
    % Compute Gaussian kernel (heat kernel-based weights)
    G(G ~= 0) = exp(-G(G ~= 0) / (2 * sigma ^ 2));     % 有待修改 因两点相等也会出现0
    G = (G+G')/2;             % Make sure distance matrix is symmetric
        
    % Construct diagonal weight matrix
    D = diag(sum(G, 2));
    
    % Compute Laplacian
    L = D - G;
    L(isnan(L)) = 0; D(isnan(D)) = 0;
	L(isinf(L)) = 0; D(isinf(D)) = 0;
    
%     % normalization on L
%     N=size(X,1);
%     for i=1:N
%         D(i,i)=D(i,i)^(-1/2);
%     end
%     L = D*L*D;