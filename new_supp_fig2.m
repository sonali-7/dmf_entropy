%% Supplementary fig 2
% Centrality measures vs Delta H

basefold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/';
filename = 'dmf_pcb_lsd_nm_i_v5.mat';
load('Deco90_Core_IgProfile.mat', 'core')
k_core = core;
simdata = load([basefold,filename]);
condnames = {'pcb','lsd'};

figfold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/figures/supp2/';

load('SC_and_5ht2a_receptors.mat')
pla = squeeze(simdata.reg_ent(:,:,1));
lsd = squeeze(simdata.reg_ent(:,:,2));
mean_pla = mean(pla,2);
mean_lsd = mean(lsd,2);
std_lsd = std(lsd,0,2);
std_pla= std(pla,0,2);
delta_h = (lsd - pla)./pla;
mean_delta_h = mean(delta_h,2); 
std_delta_h = std(delta_h,0,2); 
%% Centrality measures
sc90 = sc90./max(sc90(:))*0.2;
stren = sum(sc90)./2;
sc90g = graph(sc90);
g_deg = centrality(sc90g,'degree');
g_cc = centrality(sc90g,'closeness');
g_bc = centrality(sc90g,'betweenness');
g_prc = centrality(sc90g,'pagerank');
g_eivc = centrality(sc90g,'eigenvector');
comc = sum(expm(sc90));
% sgc = subgraph_centrality(sc90);

%% Plotting
figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','visible','on');
subplot(2,4,1)
errorbar(g_deg,mean_delta_h,std_delta_h,'ko','markerfacecolor','r');
xlabel('Degree')
grid on;
ylabel('\Delta H')

subplot(2,4,2)
errorbar(stren,mean_delta_h,std_delta_h,'ko','markerfacecolor','r');
xlabel('Strength')
grid on;

subplot(2,4,3)
errorbar(g_bc,mean_delta_h,std_delta_h,'ko','markerfacecolor','r');
xlabel('Betweeness')
grid on;

subplot(2,4,4)
errorbar(g_eivc,mean_delta_h,std_delta_h,'ko','markerfacecolor','r');
xlabel('Eigenvector Centrality')
grid on;

subplot(2,4,5)
errorbar(g_cc,mean_delta_h,std_delta_h,'ko','markerfacecolor','r');
xlabel('Closeness')
grid on;
ylabel('\Delta H')

subplot(2,4,6)
errorbar(comc,mean_delta_h,std_delta_h,'ko','markerfacecolor','r');
xlabel('Communicability')
grid on;

subplot(2,4,7)
errorbar(g_prc,mean_delta_h,std_delta_h,'ko','markerfacecolor','r');
xlabel('PageRank')
grid on;

subplot(2,4,8)
errorbar(k_core ,mean_delta_h,std_delta_h,'ko','markerfacecolor','r');
xlabel('S-core')
grid on;


print(gcf,'-dpdf',[figfold,'centrality_vs_dh.pdf'],'-r300')
%% Plotting
figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','visible','on');
subplot(2,4,1)
errorbar(g_deg,mean_pla,std_pla,'ko','markerfacecolor','r');
xlabel('Degree')
grid on;
ylabel('\Delta H')

subplot(2,4,2)
errorbar(stren,mean_pla,std_pla,'ko','markerfacecolor','r');
xlabel('Strength')
grid on;

subplot(2,4,3)
errorbar(g_bc,mean_pla,std_pla,'ko','markerfacecolor','r');
xlabel('Betweeness')
grid on;

subplot(2,4,4)
errorbar(g_eivc,mean_pla,std_pla,'ko','markerfacecolor','r');
xlabel('Eigenvector Centrality')
grid on;

subplot(2,4,5)
errorbar(g_cc,mean_pla,std_pla,'ko','markerfacecolor','r');
xlabel('Closeness')
grid on;
ylabel('\Delta H')

subplot(2,4,6)
errorbar(comc,mean_pla,std_pla,'ko','markerfacecolor','r');
xlabel('Communicability')
grid on;

subplot(2,4,7)
errorbar(g_prc,mean_pla,std_pla,'ko','markerfacecolor','r');
xlabel('PageRank')
grid on;

subplot(2,4,8)
% errorbar(sgc,mean_delta_h,std_delta_h,'ko','markerfacecolor','r');
errorbar(k_core,mean_delta_h,std_delta_h,'ko','markerfacecolor','r');
xlabel('S-core')
grid on;

print(gcf,'-dpdf',[figfold,'centrality_vs_pla_h.pdf'],'-r300')
%% Plotting centrality vs Receptors
rec2a =receptors;
figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','visible','on');
subplot(2,4,1)
plot(g_deg,rec2a,'ko','markerfacecolor','r');
xlabel('Degree')
grid on;
ylabel('5HT2A-R Density')

subplot(2,4,2)
plot(stren,rec2a,'ko','markerfacecolor','r');
xlabel('Strength')
grid on;

subplot(2,4,3)
plot(g_bc,rec2a,'ko','markerfacecolor','r');
xlabel('Betweeness')
grid on;

subplot(2,4,4)
plot(g_eivc,rec2a,'ko','markerfacecolor','r');
xlabel('Eigenvector Centrality')
grid on;

subplot(2,4,5)
plot(g_cc,rec2a,'ko','markerfacecolor','r');
xlabel('Closeness')
grid on;
ylabel('5HT2A-R Density')

subplot(2,4,6)
plot(comc,rec2a,'ko','markerfacecolor','r');
xlabel('Communicability')
grid on;

subplot(2,4,7)
plot(g_prc,rec2a,'ko','markerfacecolor','r');
xlabel('PageRank')
grid on;

subplot(2,4,8)
plot(k_core,rec2a,'ko','markerfacecolor','r');
xlabel('S-core')
grid on;

print(gcf,'-dpdf',[figfold,'centrality_vs_receptors.pdf'],'-r300')
%%

legend([p1(1) p2],{'All','Mean'},'location','southwest')