clc;
clear all;

% This code helps to reproduce the results on synthetic data in the paper
%clustering with simplicial complexes.

%The cluster indices of respective methods are given by simplicial cluster,
%motif cluster, graph cluster.

addpath("Helper_Functions")

%% Motif_Laplacian

B1 = [1 1 0 1 0 0 0 0 0 0 0 0 0;0 0 0 1 1 1 0 1 0 0 0 0 0;1 0 1 0 0 0 0 0 0 0 0 0 0;0 1 1 0 1 0 1 0 0 0 0 0 0;0 0 0 0 0 0 0 1 1 1 1 0 0;...
    0 0 0 0 0 0 0 0 0 0 1 1 0;0 0 0 0 0 1 1 0 1 0 0 0 1;0 0 0 0 0 0 0 0 0 1 0 1 1];

%Generated considering both open and closed triangles
B2 = [1 0 0 0 0 0;1 1 0 0 0 0;1 0 0 0 0 0;0 1 0 0 0 0;0 1 1 0 0 0;0 0 1 1 0 0;0 0 1 0 0 0;0 0 0 1 0 0;0 0 0 1 1 0;0 0 0 0 1 1;0 0 0 0 0 1;...
      0 0 0 0 0 1;0 0 0 0 1 0];

Motif_Laplacian = B1*(B2*B2')*B1'/4;

Motif_Adjacency =  abs(diag(diag(Motif_Laplacian)) - Motif_Laplacian);

[Motif_cluster,Motif_condv,Motif_condc,Motif_order] = SpectralPartitioning(Motif_Adjacency);


%% Graph Laplacian

Graph_Laplacian = B1*B1';

Graph_Adjacency =  abs(diag(diag(Graph_Laplacian)) - Graph_Laplacian);

[Graph_cluster,Graph_condv,Graph_condc,Graph_order] = SpectralPartitioning(Graph_Adjacency);


%% Simplicial Laplcaian

B1 = [1 1 0 1 0 0 0 0 0 0 0 0 0;0 0 0 1 1 1 0 1 0 0 0 0 0;1 0 1 0 0 0 0 0 0 0 0 0 0;0 1 1 0 1 0 1 0 0 0 0 0 0;0 0 0 0 0 0 0 1 1 1 1 0 0;...
    0 0 0 0 0 0 0 0 0 0 1 1 0;0 0 0 0 0 1 1 0 1 0 0 0 1;0 0 0 0 0 0 0 0 0 1 0 1 1];

B2 = [1 0 0 0 0;1 1 0 0 0;1 0 0 0 0;0 1 0 0 0;0 1 1 0 0;0 0 1 1 0;0 0 1 0 0;0 0 0 1 0;0 0 0 1 0;0 0 0 0 1;0 0 0 0 1;...
      0 0 0 0 1;0 0 0 0 0];

Simplicial_Laplacian = B1*(B2*B2')*B1'/4;

Simplicial_Adjacency =  abs(diag(diag(Simplicial_Laplacian)) - Simplicial_Laplacian);

[simplicial_cluster,simplicial_condv,simplicial_condc,simplicial_order] = SpectralPartitioning(Simplicial_Adjacency);

%% NMI evaluation

groundtruths = [1,1,1,1,2,2,1,2];

cluster_simplicial = 2*ones(8,1);

cluster_simplicial(simplicial_cluster) = 1;

Nmi_simplicial  = nmi(groundtruths,cluster_simplicial);
cluster_motif = 2*ones(8,1);

cluster_motif(Motif_cluster) = 1;

Nmi_motif  = nmi(groundtruths,cluster_motif);

cluster_graph = 2*ones(8,1);

cluster_graph(Graph_cluster) = 1;

Nmi_graph  = nmi(groundtruths,cluster_graph);


% %% Histogram plots of all the methods
% 
% 
% %X = {'Synthetic data','Karate network data','Polbooks','Football'};
% X = 1:4
% y = [0.56 0.56 1; 0.83 0.68 0.83; 0.57 0.55 0.59;0.90 0.89 0.91]';
% %y = [0.56 0.56 1; 0.83 0.68 0.83; 0.57 0.55 0.59;0.90 0.89 0.91]';
% barWidth = 0.5;
% %y = y';
% b = bar(X,y,barWidth);
% 
% xtips1 = b(1).XEndPoints;
% ytips1 = b(1).YEndPoints;
% labels1 = string(b(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% 
% xtips2 = b(2).XEndPoints;
% ytips2 = b(2).YEndPoints;
% labels2 = string(b(2).YData);
% text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% 
% xtips3 = b(3).XEndPoints;
% ytips3 = b(3).YEndPoints;
% labels3= string(b(3).YData);
% text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% 
% % xtips4= b(4).XEndPoints;
% % ytips4 = b(4).YEndPoints;
% % labels4 = string(b(4).YData);
% % text(xtips4,ytips4,labels4,'HorizontalAlignment','center',...
% %     'VerticalAlignment','bottom')
% 
% ylabel('NMI')
% 
% 
% 
% %%
% 
% components = {'Synthetic data','Karate network data','Polbooks','Football'};
% x = [0.56 0.56 1; 0.83 0.68 0.83; 0.57 0.55 0.59;0.90 0.89 0.91];
% figure
% y = bar(x,'grouped');
% xtips1 = y(1).XEndPoints;
% ytips1 = y(1).YEndPoints;
% labels1 = string(y(1).YData);
% text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% xtips2 = y(2).XEndPoints;
% ytips2 = y(2).YEndPoints;
% labels2 = string(y(2).YData);
% text(xtips2,ytips2,labels2,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% xtips3 = y(3).XEndPoints;
% ytips3 = y(3).YEndPoints;
% labels3 = string(y(3).YData);
% text(xtips3,ytips3,labels3,'HorizontalAlignment','center',...
%     'VerticalAlignment','bottom')
% % xtips4 = y(4).XEndPoints;
% % ytips4 = y(4).YEndPoints;
% % labels4 = string(y(4).YData);
% % text(xtips4,ytips4,labels4,'HorizontalAlignment','center',...
% %     'VerticalAlignment','bottom')
% xticklabels(components);
% grid on
% xlabel ('Design Name','fontweight','bold','FontSize',12);
% ylabel ('Penetration Level (%)','fontweight','bold','FontSize',12);
% legend('DG', 'PV');
% ylim([min(ylim) 75])
%




