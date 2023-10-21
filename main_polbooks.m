clc;
clear all;
close all;
%%

addpath("Helper_Functions")
load('Ground_truth_polbooks.mat')
S = fileread('polbooks.gml');
nodes = regexp(S, 'node.*?id (?<id>\d+).*?label\s*"(?<label>[^"]*)"', 'names');
edges = regexp(S, 'edge.*?source\s*(?<source>\d+).*?target\s*(?<target>\d+)', 'names');
ground_truth = regexp(S, 'node.*?value\s*"(?<value>[^"]*)"', 'names');
all_ids = {nodes.id};
all_names = {nodes.label};
all_sources = {edges.source};
all_targets = {edges.target};
[source_found, s] = ismember(all_sourcesecut ...
    all_ids);
nfidx = find(~source_found);
if ~isempty(nfidx)
error('Source ids not found in node list, starting from "%s"', edges(nfidx(1).source));
end
[target_found, t] = ismember(all_targets, all_ids);
nfidx = find(~target_found);
if ~isempty(nfidx)
error('Target ids not found in node list, starting from "%s"', edges(nfidx(1).target));
end
EdgeTable = table([s.', t.'], ones(length(s),1), 'VariableNames', {'EndNodes' 'Weight'});
NodeTable = table(all_names.', 'VariableNames',{'Name'});
G = graph(EdgeTable,NodeTable);
addpath('C:\Users\Smart\OneDrive - Indian Institute of Science\Desktop\Conference_papers\Eusipco_conference\Matlabcodes\Helper_functions\higher-order-organization-matlab-master\')
%%
S_graph = G.Edges(:,1);

T_graph = G.Edges(:,2);

s_g = zeros(length(all_sources),1);
t_g = zeros(length(all_targets),1);
for i = 1 : length(all_sources)
    s_g(i) = str2num(all_sources{i});
    t_g(i) = str2num(all_targets{i});
end


G1 = digraph(s_g'+1,t_g'+1);

A = adjacency(G1);

%Adjacency_polbooks = abs(diag(full(A)+full(A'))- (full(A)+full(A')));
%Adjacency_polbooks = A;
Adjacency_polbooks = abs((full(A)+full(A')));

%%

W_motif = MotifAdjacency(Adjacency_polbooks, 'M4');

W_simplicial = full(W_motif);
[ind1,ind2] = find(W_motif>0);
%ind_randsmple = randsample(1:length(ind1),4);
%ind_randsmple = [89,657,270,410];
ind_randsmple= [244	416	93	697];
%ind_randsmple= [ 110 574 244 93];

W_simplicial(ind1(ind_randsmple),ind2(ind_randsmple)) = 0;
W_simplicial(ind2(ind_randsmple),ind1(ind_randsmple)) = 0;

simplicial_Laplacian_polbooks = diag(sum(full(W_simplicial),1)) - full(W_simplicial);

%Motif_Laplacian_polbooks = diag(sum(full(W),1)) - full(W);

%[V,lambda] = eig(Motif_Laplacian_polbooks);
simplicial_Laplacian_polbooks_normalised =  (diag(sum(full(W_simplicial),1))+10*eye(105,105))^(-1/2)*simplicial_Laplacian_polbooks* (diag(sum(full(W_simplicial),1))+10*eye(105,105))^(-1/2);

[V_Simpl,lambda_Simpl] = eig(simplicial_Laplacian_polbooks_normalised);

opts = statset('Display','final');
[cluster_indices_simplicies, ctrs_simplicial] = kmeans(V_Simpl(:,2:4), 3, ...
'Replicates',200);

%%


motif_Laplacian_polbooks = diag(sum(full(W_motif),1)) - full(W_motif);

%Motif_Laplacian_polbooks = diag(sum(full(W),1)) - full(W);

%[V,lambda] = eig(Motif_Laplacian_polbooks);
motif_Laplacian_polbooks_normalised =  (diag(sum(full(W_motif),1))+15*eye(105,105))^(-1/2)*motif_Laplacian_polbooks* (diag(sum(full(W_motif),1))+15*eye(105,105))^(-1/2);
[V_motif,lambda_motif] = eig(motif_Laplacian_polbooks_normalised);

opts = statset('Display','final');
[cluster_indices_motif, ctrs_motif] = kmeans(V_motif(:,2:4), 3, ...
'Replicates',200);

%% 

nmi_simplicial = nmi(cluster_indices_simplicies,Ground_truth_labels_final);

nmi_motif = nmi(cluster_indices_motif,Ground_truth_labels_final);


%% Plotting the graphs

G = graph(Adjacency_polbooks);

P= plot(G);


x_coords = P.XData;

y_coords = P.YData;

G_gsp = gsp_graph(Adjacency_polbooks);

G_gsp.coords = [x_coords',y_coords'];

%gsp_plot_graph(G_gsp)
plot(G)
G_gsp.plottintg.vertex_size = 100;
figure;
gsp_plot_signal(G_gsp,cluster_indices_simplicies)

figure;

gsp_plot_signal(G_gsp,cluster_indices_motif);


%% For graphs

Laplacian_polbooks = diag(sum(full(Adjacency_polbooks),1)) - full(Adjacency_polbooks);

Laplacian_polbooks_normalised = diag(sum(full(Adjacency_polbooks),1))^(-1/2)*Laplacian_polbooks*diag(sum(full(Adjacency_polbooks),1))^(-1/2);

%Motif_Laplacian_polbooks = diag(sum(full(W),1)) - full(W);

%[V,lambda] = eig(Motif_Laplacian_polbooks);
[V_graphs,lambda_graph] = eig(Laplacian_polbooks_normalised);

%opts = statset('Display','final');
[cluster_indices_graph, ctrs_graph] = kmeans(V_graphs(:,2:4), 3, ...
'Replicates',200);


nmi_graph = nmi(cluster_indices_graph,Ground_truth_labels_final);


figure;

gsp_plot_signal(G_gsp,cluster_indices_graph);




%% Drawing the bar plot


