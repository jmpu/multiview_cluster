function  L = GlobalLapmatrix(Lapmatrix, a,R)

 %  by Qian Zhang
 
  nfeaSize = length(Lapmatrix);
  
  L = zeros(size(Lapmatrix{1}));
  
  for i =1:1:nfeaSize
      L = L + a(i)^R*Lapmatrix{i};
  end
  
  L = (L + L')/2;