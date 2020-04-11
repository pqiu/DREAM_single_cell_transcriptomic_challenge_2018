%%
s = importdata('download\geometry.txt');
geometry_positions = s.data;
scatter3(geometry_positions(:,1),geometry_positions(:,2),geometry_positions(:,3),10,geometry_positions(:,2),'fill')
% a little bit strange that a few points have negative values for the 2nd
% column
save raw_data geometry_positions

%%
s = importdata('download\bdtnp.txt');
bdtnp_raw = s.data;
bdtnp_gene_names = s.textdata;
s = importdata('download\binarized_bdtnp.csv');
bdtnp_binary = s.data;
% for i=1:84
%     subplot(1,2,1); scatter3(geometry_positions(:,1),geometry_positions(:,2),geometry_positions(:,3),10,bdtnp_raw(:,i),'fill')
%     title([num2str(i), ' - ', bdtnp_gene_names{i},' - raw'])
%     subplot(1,2,2); scatter3(geometry_positions(:,1),geometry_positions(:,2),geometry_positions(:,3),10,bdtnp_binary(:,i),'fill')
%     title([num2str(i), ' - ', bdtnp_gene_names{i},' - binary'])
%     pause
% end
% based on the figures it seems that the bdtnp_binary data is not very
% useful, and is in some sense more noisy than the bdtnp_raw data
save raw_data bdtnp_raw bdtnp_binary bdtnp_gene_names -append


%%
s = importdata('download\dge_raw.txt');
scRNAseq_gene_names = s.textdata;
scRNAseq_counts = s.data;

s = importdata('download\dge_normalized.txt');
scRNAseq_cell_names = regexp(s.textdata{1,:},char(9),'split');
scRNAseq_normalized = s.data;
plot(scRNAseq_counts(:,12), scRNAseq_normalized(:,12),'o')

s = importdata('download\dge_binarized_distMap.csv');
scRNAseq_binarized_distMap = s.data;
tmp = regexp(s.textdata{1},'"','split');
isequal(scRNAseq_cell_names,tmp(2:2:end))
isequal(s.textdata(2:end),bdtnp_gene_names')

save raw_data scRNAseq_gene_names  scRNAseq_cell_names scRNAseq_counts scRNAseq_normalized scRNAseq_binarized_distMap -append
