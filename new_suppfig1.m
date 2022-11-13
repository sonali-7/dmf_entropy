% New Supllementary Figure 1: Functional and anatomical grouping
% Loading data
% basefold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/';
basefold = 'E:\Matlab\fastdmf-master\newSciRep\';
simdata = load([basefold,'dmf_pcb_lsd_nm_i_v5.mat']);
N = length(simdata.params.C);
nnms = 2;
figfold = 'E:\Matlab\fastdmf-master\newSciRep\figures\supp1\';
load('Structural.mat','labelsymm')
dmf_labels = arrayfun(@(x) labelsymm(x,:),1:90,'uni',0)';
load([basefold,'yeo72aal.mat'])
nrsns = length(rsns_labs);
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
dmf_rsns_reg_id = rsns_ids(:,sort_ids);

%% Preparing data
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


%% Plotting deltaH per grouping
msize = 30;
mcols = [0.6 0 0];
% figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape'); % Maximize figure.
subplot(211)
% plot([0 0],[0 an_ngs],'k--','linewidth',2);hold on
for g=1:an_ngs
    nregs = sum(dmf_an_reg_id==g);
%     errorbar(mean_delta_h(dmf_an_reg_id==g),ones(nregs,1)*(g),[],[],...
%         std_delta_h(dmf_an_reg_id==g),std_delta_h(dmf_an_reg_id==g),...
%     'color',mcols*1.3,'CapSize',10,'linestyle','none');hold on
    plot(mean_delta_h(dmf_an_reg_id==g),ones(nregs,1)*(g),...
    '.','markersize',msize,'color',mcols);hold on

end
set(gca,'ytick',1:an_ngs,'yticklabel',[an_regs])
xlim([0 0.1])
grid on
box off;
xlabel('\Delta H')
ylim([0.5 an_ngs+0.5])

subplot(212)
% plot([0 0],[0 nrsns],'k--','linewidth',2);hold on
for g=1:nrsns 
    nregs = sum(dmf_rsns_reg_id(g,:)>0);
%     errorbar(mean_delta_h(rsns_ids(g,:)),ones(nregs,1)*(g),[],[],...
%         std_delta_h(rsns_ids(g,:)),std_delta_h(rsns_ids(g,:)),...
%     'color',mcols*1.3,'CapSize',10,'linestyle','none');hold on
    plot(mean_delta_h(dmf_rsns_reg_id(g,:)),ones(nregs,1)*(g),...
    '.','markersize',msize,'color',mcols);hold on

end
set(gca,'ytick',1:nrsns,'yticklabel',rsns_labs)
xlim([0 0.1])
grid on
box off;
xlabel('\Delta H')
ylim([0.5 nrsns+0.5])

print(gcf,'-dpdf',[figfold,'anat_and_rsns_grouping_dh.pdf'],'-r300')