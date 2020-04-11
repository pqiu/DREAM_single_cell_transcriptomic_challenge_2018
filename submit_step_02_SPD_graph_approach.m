%% load needed data
load raw_data scRNAseq_counts scRNAseq_normalized scRNAseq_gene_names bdtnp_gene_names scRNAseq_binarized_distMap

%% look at scRNAseq data, see which genes can be predicted by the landmarks
rng(100)

predictors_ind = find(ismember(scRNAseq_gene_names,bdtnp_gene_names));

data = scRNAseq_normalized;
X = [ones(1,size(data,2)); data(predictors_ind,:)];
Y = data;
Y_pred = (X'*(pinv(X')*Y'))';
pred_err = sum((Y_pred-Y).^2,2);

pred_err_null = [];
for iter = 1:10
    Y_perm = zeros(size(Y));
    [~,I] = sort(rand(size(Y)),2);
    Y_perm = Y((I-1)*size(Y,1)+repmat((1:size(Y,1))',1,size(Y,2)));
    Y_pred_null = (X'*(pinv(X')*Y_perm'))';
    pred_err_null = [pred_err_null,sum((Y_pred_null-Y_perm).^2,2)];
end

err_diff = mean(pred_err_null,2)-pred_err;
mu = mean(err_diff);
sigma = std(err_diff);
useful_genes_ind = union(find(err_diff>mu+sigma),predictors_ind);



%% PCA
data = per_gene_normalization(scRNAseq_normalized(useful_genes_ind,:));
gene_names = scRNAseq_gene_names(useful_genes_ind,:);
[num_significant_PCs, max_random_variance] = myPCA_significant_PC_by_permutation(data);
[PCA_cell_embeddings,PCA_gene_loadings,PCA_variances] = myPCA(data, num_significant_PCs);



%% kNN
distance = 'euclidean';
[IDX,D] = knnsearch(PCA_cell_embeddings,PCA_cell_embeddings,'k',5+1,'Distance',distance);
IDX(:,1) = []; D(:,1) = [];
G = knn2jaccard(IDX);
adj = G + G';
% [adj,adj2] = knn_graph(PCA_cell_embeddings,5);
% [regulation_matrix, kept_ind, components_size,components] = extract_connected_component(G+G'~=0);



%% iterative feature selection / elimination
[score, score2] = scoring_SPD_progression_similarity(adj, scRNAseq_binarized_distMap);
current_scores = sum(score);
selected_feature_ind = [];

for iter = 1:84
    scores_tmp = zeros(1,84)+Inf;
    for i=1:84
        [iter,i]
        if ismember(i, selected_feature_ind)
            continue;
        end
        scores_tmp(i) = sum(scoring_SPD_progression_similarity(adj, scRNAseq_binarized_distMap(setdiff(1:84,[i,selected_feature_ind]),:)));
    end
    [best_score_tmp,ind_to_remove] = min(scores_tmp);
    current_scores = [current_scores,best_score_tmp];
    selected_feature_ind = [ind_to_remove,selected_feature_ind];
end


%% prepare submission files
load raw_data
for num_features = [20, 40, 60]
    num_features
    % compute score based on top features
    new_score = zeros(size(bdtnp_binary,1), size(scRNAseq_binarized_distMap,2));
    for i=1:size(bdtnp_binary,1)
        for j=1:size(scRNAseq_binarized_distMap,2)
            ref = bdtnp_binary(i,selected_feature_ind(1:num_features));
            new = scRNAseq_binarized_distMap(selected_feature_ind(1:num_features),j)';
            TP = sum(ref==1 & new==1);
            TN = sum(ref==0 & new==0);
            FP = sum(ref==0 & new==1);
            FN = sum(ref==1 & new==0);
            MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
            new_score(i,j) = MCC;
        end
    end

    
    % write file
    d_output = [];
    for i=1:10:num_features
        tmp = cell2mat(strcat(bdtnp_gene_names(selected_feature_ind(i:i+10-1)),','));
        d_output = [d_output, ',', tmp(1:end-1), char(13), char(10)];
    end
    for j=1:1297
        [Y,I] = sort(new_score(:,j),'descend');
        tmp = num2str([j,I(1:10)'],'%d,');
        tmp(tmp==' ')=[];
        if j~=1297
            d_output = [d_output, tmp(1:end-1), char(13), char(10)];
        else
            d_output = [d_output, tmp(1:end-1)];
        end
    end
    
    fid = fopen(fullfile('results_graph',num2str(num_features,'%dgenes.csv')),'w'); d = fwrite(fid,d_output,'char');  fclose(fid);
end




%% functions needed

