function [PCA_cell_embeddings,PCA_gene_loadings,PCA_variances] = myPCA(data_used, num_PC)
    [U,S,V] = svds(data_used, num_PC);
    PCA_cell_embeddings = V(:,:)*S;
    PCA_gene_loadings = U(:,:);           
    PCA_variances = diag(S);
end
