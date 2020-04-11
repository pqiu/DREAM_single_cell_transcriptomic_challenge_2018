%% load needed data
load raw_data scRNAseq_counts scRNAseq_normalized scRNAseq_gene_names bdtnp_gene_names scRNAseq_binarized_distMap

load('step_04_results_update.mat', 'selected_feature_ind')
submitted_ordered_features = bdtnp_gene_names(selected_feature_ind);

[~,IA,IB] = intersect(scRNAseq_gene_names, bdtnp_gene_names);
[~,I] = sort(IB);
IA = IA(I);
useful_genes_ind84 = IA;


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
useful_genes_ind = union(find(err_diff>mu+sigma+5),predictors_ind);
useful_genes_ind321= useful_genes_ind; 


%% PCA
data = per_gene_normalization(scRNAseq_normalized(useful_genes_ind,:));
gene_names = scRNAseq_gene_names(useful_genes_ind,:);
[num_significant_PCs, max_random_variance] = myPCA_significant_PC_by_permutation(data);
[PCA_cell_embeddings,PCA_gene_loadings,PCA_variances] = myPCA(data, num_significant_PCs);



%% kNN
[adj,adj2] = knn_graph(PCA_cell_embeddings,5);



%% iterative feature selection / elimination

useful_genes_ind = useful_genes_ind84;
[scores, scores2] = scoring_SPD_progression_similarity_continuous(adj, scRNAseq_normalized(useful_genes_ind,:));
current_scores = sum(scores);
selected_feature_ind = [];

for iter = 1:length(useful_genes_ind)
    scores_tmp = zeros(1,length(useful_genes_ind))+Inf;
    for i=1:length(useful_genes_ind)
        [iter,i]
        if ismember(useful_genes_ind(i), selected_feature_ind)
            continue;
        end
        scores_tmp(i) = sum(scoring_SPD_progression_similarity_continuous(adj, scRNAseq_normalized(setdiff(useful_genes_ind,[useful_genes_ind(i),selected_feature_ind]),:)));
    end
    [best_score_tmp,ind_to_remove] = min(scores_tmp);
    current_scores = [current_scores,best_score_tmp];
    selected_feature_ind = [useful_genes_ind(ind_to_remove),selected_feature_ind];
end

selected_feature_ind84 = selected_feature_ind;
[~, a, b] = intersect(scRNAseq_gene_names(selected_feature_ind84), submitted_ordered_features);
figure(1)
plot(a,b,'o'); % if good correlation, the continuous version of feature selection is the same as the MCC version



%%
useful_genes_ind = useful_genes_ind321;
[scores, scores2] = scoring_SPD_progression_similarity_continuous(adj, scRNAseq_normalized(useful_genes_ind,:));
current_scores = sum(scores);
selected_feature_ind = [];

for iter = 1:length(useful_genes_ind)
    scores_tmp = zeros(1,length(useful_genes_ind))+Inf;
    for i=1:length(useful_genes_ind)
        [iter,i]
        if ismember(useful_genes_ind(i), selected_feature_ind)
            continue;
        end
        scores_tmp(i) = sum(scoring_SPD_progression_similarity_continuous(adj, scRNAseq_normalized(setdiff(useful_genes_ind,[useful_genes_ind(i),selected_feature_ind]),:)));
    end
    [best_score_tmp,ind_to_remove] = min(scores_tmp);
    current_scores = [current_scores,best_score_tmp];
    selected_feature_ind = [useful_genes_ind(ind_to_remove),selected_feature_ind];
end

selected_feature_ind321 = selected_feature_ind;
[~, a, b] = intersect(scRNAseq_gene_names(selected_feature_ind321), submitted_ordered_features);
figure(2)
plot(a,b,'o'); % if good correlation, the continuous version of feature selection is the same as the MCC version


%% 
save submit_step_03_post_conference_rank_321_results

%%
load('ground_true_distMap_score.mat')
load('raw_data.mat')
load('submit_step_03_post_conference_rank_321_results.mat','selected_feature_ind321','selected_feature_ind84')
[~,I] = sort(score_ground_truth,'descend');
scRNAseq_positions = geometry_positions(I(1,:),:);

%%
for i =1:321
    gene = scRNAseq_gene_names(selected_feature_ind321(i));
    subplot(1,2,1);
    scatter(scRNAseq_positions(:,1),scRNAseq_positions(:,3),30,scRNAseq_normalized(ismember(scRNAseq_gene_names,gene),:),'fill')
    title(gene{1})
    if ismember(gene,bdtnp_gene_names)
        subplot(1,2,2);
        scatter(geometry_positions(:,1),geometry_positions(:,3),30,bdtnp_raw(:,ismember(bdtnp_gene_names,gene)),'fill')
    else
        subplot(1,2,2);
        plot(0)
    end    
    pause
end