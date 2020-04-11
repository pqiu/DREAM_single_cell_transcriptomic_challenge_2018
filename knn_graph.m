function [adj,adj2] = knn_graph(X,K,working_mode,exclude_adj)
% [adj,adj2] = knn_graph(X,K,working_mode,exclude_adj)
%     KNN graph based on Euclidian distances
%     use (X)(n x p) to nearest neighbor graph, with number of neighbor being K
%     each row is one node, (n nodes in total)
%     p dimensional variable-space, 
%     working_mode : 'euclidean' (default)
%                    'corr'
%                    'abs_corr'
% 
%     out:Xmst (objects-1 x 2) link set between 'objects' indexed as rows in X
%         adj adjacency matrix



if ~exist('working_mode')
    working_mode = 'euclidean';
end
if isempty(intersect({'euclidean','corr','abs_corr'}, working_mode))
    working_mode = 'euclidean';
end
if ~isempty(intersect({'corr','abs_corr'}, working_mode))
    X = per_gene_normalization(X); % this makes computing correlation easier
    X = X./norm(X(1,:));
end

if ~exist('exclude_adj') || isempty(exclude_adj) 
    exclude_adj = sparse(size(X,1),size(X,1));
end

[nX,mX] = size(X);

adj = sparse(size(X,1),size(X,1)); 
adj2 = sparse(size(X,1),size(X,1)); 
count = 0; %fprintf('constructing a total of %d MST edges ... %6d', size(X,1)-1,count);
for i=1:nX
    if isequal(working_mode, 'euclidean')
        dist = comp_dist_euclidean(X,i,1:nX); 
    elseif isequal(working_mode, 'corr')
        dist = comp_dist_corr(X,i,1:nX);
    elseif isequal(working_mode, 'abs_corr')
        dist = comp_dist_abs_corr(X,i,1:nX); 
    end
    dist(i) = max(dist)+1; 
    dist = dist + exclude_adj(i,:).*(max(dist)+1);    
    dist = full(dist);
    
    [Y,I] = sort(dist,'ascend');
    adj(i,I(1:K))=1;
    adj(i,i)=0;
    adj(i,exclude_adj(i,:)==1)=0;
    
    adj(:,i) = (adj(:,i) | (adj(i,:)'));
    
    adj2(i,adj(i,:)==1) = dist( adj(i,:)==1);
    adj2(adj(i,:)==1,i) = dist( adj(i,:)==1)';
end
 
return




function dist = comp_dist_euclidean(X,ind1,ind2)
dist = zeros(length(ind1),length(ind2));
for i=1:length(ind1)
    dist(i,:) = sqrt(sum((repmat(X(ind1(i),:),length(ind2),1) - X(ind2,:)).^2,2)); 
end
return


function dist = comp_dist_corr(X,ind1,ind2)
dist = zeros(length(ind1),length(ind2));
corr = X(ind1,:)*X(ind2,:)';
dist = 1-corr; 
return


function dist = comp_dist_abs_corr(X,ind1,ind2)
dist = zeros(length(ind1),length(ind2));
corr = X(ind1,:)*X(ind2,:)';
dist = 1-abs(corr); 
return

