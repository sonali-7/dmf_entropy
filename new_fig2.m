% NEW FIGUREs OF SCIREP: corrigendum
%% % Paper fig # 2
% 1.- H_NM vs H_rs
% 2.- bar \delta H vs regions sorted in ascending order
% 3.- \delta H vs Receptor density

% Loading data
% basefold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/';
basefold = 'E:\Matlab\fastdmf-master\newSciRep\';
simdata = load([basefold,'dmf_pcb_lsd_nm_i_v5.mat']);
N = length(simdata.params.C);
nnms = 2;
figfold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/figures/fig2/';
load('Structural.mat','labelsymm')
dmf_labels = arrayfun(@(x) labelsymm(x,:),1:90,'uni',0)';
%% Creating anatomical labels
sort_ids = cat(2,[1:2:N],fliplr([2:2:N])); % sorting id to deco to aa
an_regs = {'Cingulate','Frontal','Limbic','Occipital','Parietal',...
    'Sensorimotor','Subcortical','Temporal'};
an_ngs = length(an_regs);
an_regs_id = {[31:36],[3:16 19:28],[37:42],[43:56],[59:70],...
    [1 2 17 18 57 58],[71:78],[29:30, 79:90]};
an_reg_id = zeros(N,1);
for g=1:an_ngs
    an_reg_id(an_regs_id{g}) = g;
end
dmf_an_reg_id=an_reg_id(sort_ids);

%% 1.- Plotting  H_NM vs H_rs
anat_cols = flipud(othercolor('Dark25',an_ngs+1));
anat_cols =anat_cols(2:end,:);
mean_reg_ent = squeeze(mean(simdata.reg_ent,2));
std_reg_ent = squeeze(std(simdata.reg_ent,0,2));
pla_reg_ent = squeeze(simdata.reg_ent(:,:,1));
nm_reg_ent = squeeze(simdata.reg_ent(:,:,2));
delta_h = (nm_reg_ent - pla_reg_ent)./pla_reg_ent;
mean_delta_h = mean(delta_h,2);
std_delta_h = std(delta_h,0,2);
[sort_delta_h,dh_sort_id] = sort(mean_delta_h);
sort_std_dh = std_delta_h(dh_sort_id);


min_ent = min(simdata.reg_ent(:));
max_ent = max(simdata.reg_ent(:));
thislims = [2 2.6];

%%
figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','visible','on');
subplot(121)

plot(mean_reg_ent(:,1),mean_reg_ent(:,2),'o','color','k','markersize',10,...
    'markerfacecolor',[0.7 0 0]);hold on

plot([min_ent-0.1 max_ent+0.1],[min_ent-0.1 max_ent+0.1],'k--');
xlim(thislims );
ylim(thislims );
xlabel('Resting State Region Entropy')
ylabel('5HT2A-R Region Entropy')
box off
grid on;
axis square

subplot(122)
for g=1:an_ngs
    plot(mean_reg_ent(dmf_an_reg_id==g,1),mean_reg_ent(dmf_an_reg_id==g,2),'o','color','k','markersize',10,...
        'markerfacecolor',anat_cols(g,:));hold on
    
end
plot([min_ent-0.1 max_ent+0.1],[min_ent-0.1 max_ent+0.1],'k--');
xlim(thislims );
ylim(thislims );
xlabel('Resting State Region Entropy')
ylabel('5HT2A-R Region Entropy')
box off
grid on;
axis square;
% alpha 0.5
legend(an_regs)

print(gcf,'-dpdf',[figfold,'scatter.pdf'],'-r300')

%% 1.1 Histograms of regions entropy for both conditions
% Statistical test
mean_brain_ent = squeeze(mean(simdata.reg_ent,2));
mean_brain_ent_p_val = signrank(mean_reg_ent(:,1),mean_reg_ent(:,2),'tail','left');
%% Cohens d for each simulation
nreps = simdata.nreps;
cohen_d = zeros(nreps,1);
tic
for r=1:nreps
    cohen_d(r) = computeCohen_d(simdata.reg_ent(:,r,2),simdata.reg_ent(:,r,1));
end
toc

%% Building and Plotting histogram
edges = linspace(min_ent-0.1,max_ent+0.1,100);
% edges = linspace(thislims(1),thislims(2),100);
[pla_nel,pla_e] =  histcounts(pla_reg_ent(:),edges,'Normalization','probability');
[nm_nel,nm_e] =  histcounts(nm_reg_ent(:),edges,'Normalization','probability');

%
figure
a1=area(pla_e(2:end),smooth(pla_nel,5),'linestyle','none','FaceAlpha',0.7,'FaceColor',[0 0 0.7]);hold on
a2=area(nm_e(2:end),smooth(nm_nel,5),'linestyle','none','FaceAlpha',0.7,'FaceColor',[0.7 0 0]);hold on
plot([mean(pla_reg_ent(:)) mean(pla_reg_ent(:))],[0 0.14],'--','color',[0 0 0.7],'linewidth',2);hold on
plot([mean(nm_reg_ent(:)) mean(nm_reg_ent(:))],[0 0.14],'--','color',[0.7 0 0],'linewidth',2);hold on
xlim([min_ent-0.1 max_ent+0.1]);
xlim([2 2.7]);
ylim([0 0.07])
grid on
xlabel('Region Entropy')
ylabel('Normalized Frequency')
box off
legend([a1,a2], {'Resting State','5HT2A-R'},'location','northeast')
axis square

print(gcf,'-dpdf',[figfold,'histogram.pdf'],'-r300')

%% 2.- Averaging between hemispheres
labelsymm_h = cellfun(@(x) x(3:end),dmf_labels(1:45),'uni',0);
homo_d_h = (delta_h(1:45,:) + flipud(delta_h(46:end,:)))./2;
h_mean_delta_h = mean(homo_d_h,2);
h_std_delta_h = std(homo_d_h,0,2);
[h_sort_delta_h,h_dh_sort_id] = sort(h_mean_delta_h);
h_sort_std_dh = std_delta_h(h_dh_sort_id);

figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape'); % Maximize figure.
bar(1:N/2,h_sort_delta_h,'edgecolor','none','facecolor',[0.6 0 0]);hold on
errorbar(1:N/2,h_sort_delta_h,[],h_sort_std_dh,'linestyle','none',...
    'linewidth',2,'color',[0.3 0.3 0.3]);
set(gca,'xtick',1:45,'xticklabel',labelsymm_h(h_dh_sort_id))
xtickangle(90)
grid on;
box off
set(gca,'Xgrid','off')
ylabel('Entropy Change (\Delta H)')
ylim([0 0.1])
print(gcf,'-dpdf',[figfold,'homotopic_average_bar_plot_sorted_delta_h.pdf'],'-r300')

%% new figure 2C -> cerebro con Delta H coloreado
load('aal_node_pos.mat')
node_pos = v;
sort_ids = cat(2,[1:2:N],fliplr([2:2:N]));
delta_h_aal = mean_delta_h;
delta_h_aal(sort_ids) = mean_delta_h;
sel_nodes = [65,19,43];
% sel_nodes = [10,11,12];
figure;
method  = 'euclidean';
% 'raycast', 'euclidean' or 'spheres'
% mcolors = flipud(othercolor('Spectral4',N));
% mcolors = flipud(othercolor('RdBu11',N));
mcolors  = othercolor('Reds9',5000);
% mcolors = spring(N);
meshclims = [0 0.1];
nodecols=cat(1,[0 0 0.7],[0.5 0.5 0.5],[0.7 0 0]);
aa=atemplate('overlay',delta_h_aal,'meshcols',mcolors,'method',method,...
    'hemi','both','nocolbar','meshclims',meshclims);%,...
%     'nodes',sel_nodes,'nodecols',nodecols);
%
material shiny
hold on;
%% set views:
v=cell(3,1);
v{1} = [270 0]; % L
v{2} = [0  90]; % Topo
v{3} = [90  0]; % R
zview={'L','T','R'};
nviews = length(v);
hemis = {'both','left','right'};
nhemis = length(hemis);
% Extracting zenital plot of brain

for i=1:nviews % views
    for h=1:nhemis
        figure();
        atemplate('overlay',delta_h_aal,'method',method,'meshclims',meshclims,...
            'hemi',hemis{h},'nocolbar','meshcols',mcolors); % g
        hold on;        
        material shiny
        view(gca,v{i});
        cl =camlight('headlight');
        set(cl,'Color',[0.3 0.3 0.3])
        print(gcf,'-dpng',[figfold,'delta_h_',hemis{h},'_',zview{i},'.png'],'-r300')        
        close gcf
        
    end
    
    
end

