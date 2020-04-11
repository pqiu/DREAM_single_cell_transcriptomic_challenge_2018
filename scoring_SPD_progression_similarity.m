function scores = scoring_SPD_progression_similarity(adj, scRNAseq_binarized_distMap)


tic
fprintf('scoring_SPD_progression_similarity ... %4d',0);
for i=1:size(scRNAseq_binarized_distMap,2)
    fprintf('\b\b\b\b%4d',i)
    ref = scRNAseq_binarized_distMap(:,i);
    TP = sum(repmat(ref,1,size(scRNAseq_binarized_distMap,2)) & scRNAseq_binarized_distMap);
    TN = sum(repmat(~ref,1,size(scRNAseq_binarized_distMap,2)) & ~scRNAseq_binarized_distMap);
    FP = sum(repmat(~ref,1,size(scRNAseq_binarized_distMap,2)) & scRNAseq_binarized_distMap);
    FN = sum(repmat(ref,1,size(scRNAseq_binarized_distMap,2)) & ~scRNAseq_binarized_distMap);
    mcc_score(i,:) = (TP.*TN - FP.*FN)./sqrt((TP+FP).*(TP+FN).*(TN+FP).*(TN+FN));
end
toc

% tic
% fprintf('scoring_SPD_progression_similarity ... %4d',0);
% for i=1:size(scRNAseq_binarized_distMap,2)
%     fprintf('\b\b\b\b%4d',i)
%     for j=i:size(scRNAseq_binarized_distMap,2)
%         ref = scRNAseq_binarized_distMap(:,i);
%         new = scRNAseq_binarized_distMap(:,j);
%         TP = sum(ref==1 & new==1);
%         TN = sum(ref==0 & new==0);
%         FP = sum(ref==0 & new==1);
%         FN = sum(ref==1 & new==0);
%         MCC = (TP*TN - FP*FN)/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
%         mcc_score(i,j) = MCC;
%         mcc_score(j,i) = MCC;
%     end
% end
% toc

[Y,I] = sort(mcc_score,2,'descend');

scores = zeros(1,size(scRNAseq_binarized_distMap,2));
for i=1:size(scRNAseq_binarized_distMap,2)
    scores(i) = sum(find(ismember(I(i,:),find(adj(i,:)==1))));
end

