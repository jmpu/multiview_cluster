
% The original view
OR = csvread('languages\EN-EN.feature');
% The translated view
TR = csvread('languages\EN-FR.feature');
% The translated view
SR = csvread('languages\EN-GR.feature');
% topic labels
label_raw = csvread('languages\EN.label');

% number of topics
K = 6;
% number of documents to sample
N = 1200;
tmp = randperm(length(label_raw));
selected = tmp(1:N);
label = label_raw(selected);
[label,ind] = sort(label,'ascend');
selected = selected(ind);

for i=1:K
    fprintf('%d documents sampled from Topic %d\n',nnz(label==i),i);
end

%construct the graph
S1 = sparse(OR(:,1),OR(:,2),OR(:,3));
S1 = S1(selected,:);
S1 = full(S1);
for i=1:size(S1,1)
    S1(i,:)=S1(i,:)/norm(S1(i,:));
end
% A1 = S*S';

S2 = sparse(TR(:,1),TR(:,2),TR(:,3));
S2 = S2(selected,:);
S2 = full(S2);
for i=1:size(S2,1)
    S2(i,:)=S2(i,:)/norm(S2(i,:));
end
% A2 = S*S';


S3 = sparse(SR(:,1),SR(:,2),SR(:,3));
S3 = S3(selected,:);
S3 = full(S3);
for i=1:size(S3,1)
    S3(i,:)=S3(i,:)/norm(S3(i,:));
end
% A3 = S*S';

