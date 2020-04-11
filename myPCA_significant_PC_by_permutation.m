function [num_significant_PCs, max_random_variance] = myPCA_significant_PC_by_permutation(data_used, num_PC, perm_iter, isdisplay)
    if ~exist('num_PC')
        num_PC = min(size(data_used)); 
    end
    if ~exist('perm_iter')
        perm_iter = 30; 
    end
    
    if ~exist('isdisplay')
        isdisplay = 0; 
    end

    
    highest_random_variance = zeros(1,perm_iter);
    fprintf('Random permutation %d iter for computing signficant PC ... %5d ', perm_iter, 0);
    tic
    for i=1:perm_iter
        perm_data = zeros(size(data_used));
        [~,I] = sort(rand(size(data_used)),2);
        perm_data = data_used((I-1)*size(data_used,1)+repmat((1:size(data_used,1))',1,size(data_used,2)));
        [U,S,V] = svds(perm_data,num_PC);
        highest_random_variance(i) = max(diag(S));
        fprintf('\b\b\b\b\b\b%5d ', i);
    end
    toc
    max_random_variance = max(highest_random_variance);

    
    [U,S,V] = svds(data_used, num_PC);
    PCA_cell_embeddings = V(:,:)*S;
    PCA_gene_loadings = U(:,:);           
    PCA_variances = diag(S);
    
    num_significant_PCs = sum(PCA_variances>=max_random_variance);

end

