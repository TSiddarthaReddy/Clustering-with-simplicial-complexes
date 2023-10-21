clc;
clear all;
close all;


%%  Loading the data

load('Karatey_adjacency.mat');
load('Karatey_labels.mat');
addpath("Helper_Functions")

%% Listing the triangles based on Motif_Adjacency

W_motif = MotifAdjacency(adj_matrix, 'M4');
W_simplicial = full(W_motif);
[ind1,ind2] = find(W_motif>0);

% Indices of the edges that are to be removed as explained in the paper this is
% exhaustive search method 
ind_randsmple = [55,54];
W_simplicial(ind1(ind_randsmple),ind2(ind_randsmple)) = 0;
W_simplicial(ind2(ind_randsmple),ind1(ind_randsmple)) = 0;

%[LCC, lcc_inds] = LargestConnectedComponent(W);
[cluster_simplicial_aux, condv_simp_aux, condc_sim_aux] = SpectralPartitioning(W_simplicial);


% %% Check with Laplacian eigen vectors
% 
% Laplacian_motif = diag(sum(full(W_motif),1)) - full(W_motif);
% 
% Laplacian_graph = diag(sum(adj_matrix)) - adj_matrix;
% 
% Laplacian_simplicial = diag(sum(full(W_simplicial),1)) - full(W_simplicial);
% 
% [V,lambda] = eig(Laplacian_simplicial);

%%
clusters_simplicial  = zeros(length(ground_truth),1);

clusters_simplicial(cluster_simplicial_aux) = 1;

%miss_classfied_nodes_Simplicial = nnz(clusters_simplicial-ground_truth');

%B = setdiff(1:34,cluster_simplicial_aux);
NMI_Simplicial = nmi(clusters_simplicial,ground_truth');
%% For motifs

[cluster_motif_aux, condv_motif_aux, condc_motif_aux] = SpectralPartitioning(W_motif);

clusters_motif  = zeros(length(ground_truth),1);

clusters_motif(cluster_motif_aux) = 1;

%miss_classfied_nodes_motif = nnz(clusters_motif-ground_truth');

%B = setdiff(1:34,cluster_motif_aux);
NMI_motif = nmi(clusters_motif,ground_truth');
%% For graphs

[cluster_graph_aux, condv_graph_aux, condc_graph_aux] = SpectralPartitioning(adj_matrix);

clusters_graph  = zeros(length(ground_truth),1);
clusters_graph(cluster_graph_aux) = 1;
NMI_graph = nmi(clusters_graph,ground_truth');
%% Plotting the clusters in zarachy 

G = graph(adj_matrix);

P= plot(G);


x_coords = P.XData;

y_coords = P.YData;

G_gsp = gsp_graph(adj_matrix);

G_gsp.coords = [x_coords',y_coords'];

%gsp_plot_graph(G_gsp)
plot(G)
G_gsp.plottintg.vertex_size = 100;
figure;
gsp_plot_signal(G_gsp,clusters_simplicial)
colorbar off
title('clusters using simplicial adjacency')
figure;


gsp_plot_signal(G_gsp,clusters_motif)
colorbar off
title('clusters using motif adjacency')
figure;

gsp_plot_signal(G_gsp,ground_truth');
colorbar off
title('Groundtruth')
