clc;
clear all;
close all;

addpath("Helper_Functions")
%% 

load('football.txt');
load('Ground_truth_football.mat');
s= football(:,1)+1;
t= football(:,2)+1;
G = digraph(s,t);
A = adjacency(G);
Adjacency_football = full(A)+full(A');

%% Using the laplacian


L = diag(sum(full(Adjacency_football),1)) - full(Adjacency_football);
L_normalised  = diag(sum(full(Adjacency_football),1))^(-1/2)*L*diag(sum(full(Adjacency_football),1))^(-1/2);
[V,lambda] = eig(L_normalised);
%opts = statset('Display','final');
[cluster_indices_graph, ctrs] = kmeans(V(:,2:13), 12, ...
'Replicates',200);  

%% Evaluating the NMI

NMI_graph = nmi(cluster_indices_graph,ground_truth_football);


%% Using the motif laplacian
W_motif = MotifAdjacency(full(Adjacency_football), 'M4');

%W_simplicial = full(W);
%[ind1,ind2] = find(W>0);

L_motif = diag(sum(full(W_motif),1)) - full(W_motif);
L_motif_normalised = diag(sum(full(W_motif),1))^(-1/2)*L_motif*diag(sum(full(W_motif),1))^(-1/2);
[V_mot,lambda_mot] = eig(L_motif_normalised);
%opts = statset('Display','final');
[cluster_indices_motif, ctrs] = kmeans(V_mot(:,2:13), 12, ...
'Replicates',100); 

NMI_motif = nmi(cluster_indices_motif,ground_truth_football);


%%
W = MotifAdjacency(Adjacency_football, 'M4');

W_simplicial = full(W);
[ind1,ind2] = find(W>0);
%montecarlo = 100;
%store = cell(montecarlo,1);
%for montecarlo = 1 : 100
%ind_randsmple = randsample(1:length(ind1),5);
%ind_randsmple = [163 586 25 501 294];
%ind_randsmple = [1,160,678,257,195];
%ind_randsmple = [67,101];
ind_randsmple= [858	690	744	780	81];
W_simplicial(ind1(ind_randsmple),ind2(ind_randsmple)) = 0;
W_simplicial(ind2(ind_randsmple),ind1(ind_randsmple)) = 0;

L_simplicial = diag(sum(full(W_simplicial),1)) - full(W_simplicial);
L_simplicial_normalised = diag(sum(full(W_simplicial),1))^(-1/2)*L_simplicial* diag(sum(full(W_simplicial),1))^(-1/2);
[V_sim,lambda_sim] = eig(L_simplicial_normalised);
%opts = statset('Display','final');
[cluster_indices_simplicial, ctrs] = kmeans(V_sim(:,2:13), 12, ...
'Replicates',400); 
%print off;
NMI_simplicial= nmi(cluster_indices_simplicial,ground_truth_football);
%store{montecarlo} = ind_randsmple;
%W_simplicial = full(W);
%end


%% Plotting

G = graph(Adjacency_football);

P= plot(G);


x_coords = P.XData;

y_coords = P.YData;

G_gsp = gsp_graph(Adjacency_football);

G_gsp.coords = [x_coords',y_coords'];

%gsp_plot_graph(G_gsp)
plot(G)
G_gsp.plottintg.vertex_size = 100;
figure;
gsp_plot_signal(G_gsp,cluster_indices_simplicial)

figure;

gsp_plot_signal(G_gsp,cluster_indices_motif);

figure;

gsp_plot_signal(G_gsp,ground_truth_football);