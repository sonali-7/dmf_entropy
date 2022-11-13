% New Fig3 for SciRep: Corrigendum
% basefold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/';
% basefold = 'E:\Matlab\fastdmf-master\newSciRep\';
basefold = 'E:\Matlab\fastdmf-master\newSciRep\';
% filename = 'dmf_pcb_lsd_music_reg_ents.mat';
filename = 'dmf_pcb_lsd_nm_i_v5.mat';
load('Deco90_Core_IgProfile.mat', 'ignitionProfile', 'core')
load('Structural.mat','labelsymm')
dmf_labels = arrayfun(@(x) labelsymm(x,:),1:90,'uni',0)';
load('yeo72aal.mat')
simdata = load([basefold,filename]);
N=90;
condnames = {'pcb','lsd'};
stren = sum(simdata.params.C)./2;receptors = simdata.params.receptors;
comc = sum(expm(simdata.params.C));
sort_ids = cat(2,[1:2:N],fliplr([2:2:N])); % sorting id to deco to aa
rsns_ids = rsns_ids(:,sort_ids);
nreps = simdata.nreps;

%% Grouping AAL regions
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

%% Using Quadratic for S and linear for R
% --------------------------Selection of connectivity feature
sel_con = stren;
% --------------------------Preparing data
mean_ent= squeeze(mean(simdata.reg_ent,2));
deltah= (simdata.reg_ent(:,:,2)-simdata.reg_ent(:,:,1))./simdata.reg_ent(:,:,1);
mean_dh = mean(deltah,2);
std_dh =  std(deltah,0,2);

%% LME formula and data
varnames = {'dH';'C';'R';'G'};
lme_formula = 'dH ~R:G + C^2:G '; % 

tab = table(mean_dh,sel_con',receptors,'VariableNames',varnames(1:3));

% ------------------------- Anatomical prior
sel_anat = dmf_an_reg_id==4 ; % Occipital 
not_sel = not(sel_anat);
% --------------------------functions and options
ft = fittype('b*x^2+c*x+a');
options = fitoptions(ft);
options.Upper = [0 inf inf];
options.Lower = [-inf -inf -inf];
options.Robust = 'on';

ft2 = fittype('c*x+a');
options2 = fitoptions(ft2);
options2.Upper = [inf inf];
options2.Lower = [0 -inf];
options2.Robust = 'on';

% Residuals factor to assign data to groups
sel_resfac =1;
% --------------------------- Fitting
[fit1] = fit(sel_con(sel_anat)',mean_dh(sel_anat),ft,options);
pred_dh_s = fit1(sel_con);

[fit2] = fit(receptors(not_sel),mean_dh(not_sel),ft2,options2);
pred_dh_r = fit2(receptors);
% -------------------------- Comparison between residuals
res_s = (pred_dh_s - mean_dh);
res_r = (pred_dh_r - mean_dh);
is_s = abs(res_s)<(sel_resfac.*abs(res_r)) | sel_anat & not_sel; % lower residuals for s, belong to s
is_r = not(is_s);
% --------------------------- Searching the group that maximizes the R² of a LME with 3-way separation

nreps = 50000;
r_regs = find(is_r);
nselregs = length(r_regs);

size_offset = 5;
minpts = size_offset:nselregs-5;
nmpt = numel(minpts);
max_sel_ids = cell(nmpt,1);
best_rsq_per_s = zeros(nmpt,1);
lme_rsq = zeros(nmpt,1);
tic
parfor s=1:nmpt
    subset_ids = cell(nreps,1);
    subset_rsq = zeros(nreps,1);
    this_size = minpts(s);
    for r=1:nreps
        % Random subset size        
        rand_ids = r_regs(randperm(nselregs,this_size)); % random ids from r_regs
        subset_ids{r} = rand_ids;
        
        % Computing correlation
        subset_rsq(r) = corr2(receptors(rand_ids),mean_dh(rand_ids));
        
    end
    % keeping ids of best    
    [best_rsq_per_s(s),max_rsq_id] = max(subset_rsq);
    max_sel_ids{s} = subset_ids{max_rsq_id};
    
    % Fitting LME    
    group_var = zeros(N,1);
    group_var(is_s)=1; % 1 for S, 0 for Null, -1 for R
    group_var(max_sel_ids{s})=-1;
    cat_var = categorical(group_var);
    sel_tab = [tab,table(cat_var,'VariableNames',{'G'})];
    % sel_tab = tab;
    try
        lme_fit= fitlme(sel_tab,lme_formula,'DummyVarCoding','full');
        lme_rsq(s)=lme_fit.Rsquared.Ordinary;
    end
    
    
end
toc
%

figure
plot(minpts,lme_rsq,'o-')
% Extarcting best
nreps = simdata.nreps;
[best_rsq, best_rsq_id]= max(lme_rsq);
this_sel_ids = max_sel_ids{best_rsq_id};

% Fitting LME to each simulation
group_var = zeros(N,1);
group_var(is_s)=1; % 1 for S, 0 for Null, -1 for R
group_var(this_sel_ids)=-1;
cat_var = categorical(group_var);
lme_rsq_rep = zeros(nreps,1);


tic
for r=1:nreps
    tab = table(deltah(:,r),sel_con',receptors,'VariableNames',varnames(1:3));
    sel_tab = [tab,table(cat_var,'VariableNames',{'G'})];
    lme_fit= fitlme(sel_tab,lme_formula,'DummyVarCoding','full');
    lme_rsq_rep(r)=lme_fit.Rsquared.Ordinary;
end
toc
%% Fitting using averages
lme_formula = 'dH ~R:G + C^2:G'; % 
strenxx = linspace(0,0.5,1000);
tab = table(mean_dh,sel_con',receptors,'VariableNames',varnames(1:3));
sel_tab = [tab,table(cat_var,'VariableNames',{'G'})];
lme_fit= fitlme(sel_tab,lme_formula,'DummyVarCoding','full');
% Fitting to whole dataset and to grouped variables
nreps = simdata.nreps;
ngs = 3;
gvals = [-1 0 1];
q_s_dh_rsq = zeros(nreps,1);
l_r_dh_rsq = q_s_dh_rsq;
q_s_dh_rsq_g = zeros(nreps,ngs);
l_r_dh_rsq_g = q_s_dh_rsq_g ;

tic
parfor r=1:nreps
    % Fit options
    ft = fittype('b*x^2+c*x+a');
    options = fitoptions(ft);
    options.Upper = [0 inf inf];
    options.Lower = [-inf -inf -inf];
    options.Robust = 'off';
    
    ft2 = fittype('c*x+a');
    options2 = fitoptions(ft2);
    options2.Upper = [inf inf];
    options2.Lower = [0 -inf];
    options2.Robust = 'off';
    % Fitting to strength
    [~,q_s_error] = fit(sel_con',deltah(:,r),ft,options);
    q_s_dh_rsq(r) =q_s_error.rsquare;
    % Fitting to receptors
    [~,l_r_error] = fit(receptors,deltah(:,r),ft2,options2);
    l_r_dh_rsq(r)=l_r_error.rsquare;
    for g=1:ngs        
        % Fitting to strength
        [~,q_s_error] = fit(sel_con(group_var==gvals(g))',deltah(group_var==gvals(g),r),ft,options);
        q_s_dh_rsq_g(r,g) =q_s_error.rsquare;
        % Fitting to receptors
        [~,l_r_error] = fit(receptors(group_var==gvals(g)),deltah(group_var==gvals(g),r),ft2,options2);
        l_r_dh_rsq_g(r,g)=l_r_error.rsquare;
    end
end
toc

% Fitting to averages for plotting
q_s_mean_dh_fit = cell(4,1); % whole and grouping
l_r_mean_dh_fit = q_s_mean_dh_fit;
% Fitting to strength
q_s_mean_dh_fit{1} = fit(sel_con',mean_dh,ft,options);
% Fitting to receptors
l_r_mean_dh_fit{1} = fit(receptors,mean_dh,ft2,options2);
for g=1:ngs
    % Fitting to strength
    q_s_mean_dh_fit{g+1} = fit(sel_con(group_var==gvals(g))',mean_dh(group_var==gvals(g)),ft,options);
    
    % Fitting to receptors
    l_r_mean_dh_fit{g+1} = fit(receptors(group_var==gvals(g)),mean_dh(group_var==gvals(g)),ft2,options2);
    
end

%
s_xx = linspace(0,max(sel_con)*1.1,1000);
r_xx = linspace(0,1,1000);
s_yy = feval(q_s_mean_dh_fit{1},s_xx);
r_yy = feval(l_r_mean_dh_fit{1},r_xx);

s_yy_1 = feval(q_s_mean_dh_fit{4},s_xx);
s_yy_0 = feval(q_s_mean_dh_fit{3},s_xx);
r_yy_m_1 = feval(l_r_mean_dh_fit{2},r_xx);


%% Plottin deltaH vs Stren and Receptors and the best line
% figfold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/figures/fig3/';
figfold = 'E:\Matlab\fastdmf-master\newSciRep\figures\fig3\';
msize = 30;

figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape'); % Maximize figure.
subplot(1,2,1)
errorbar(sel_con(group_var==1),mean_dh(group_var==1),std_dh(group_var==1),'b.','markersize',msize);hold on 
errorbar(sel_con(group_var==-1),mean_dh(group_var==-1),std_dh(group_var==-1),'r.','markersize',msize);hold on 
errorbar(sel_con(group_var==0),mean_dh(group_var==0),std_dh(group_var==0),'c.','markersize',msize);hold on 
% Plotting lines
p1=plot(s_xx,s_yy,'k-'); % whole data
p2=plot(s_xx,s_yy_1,'b-'); % s1-group
p3=plot(s_xx,s_yy_0,'c-'); % s2-group
xlim([0.99*min(sel_con),max(sel_con)*1.01])
grid on
axis square
xlabel('Strength')
ylabel('\Delta h')
ylim([0 0.1])

legend([p1 p3 p2],{['All, R²=',num2str(mean(q_s_dh_rsq)),' \pm ',num2str(std(q_s_dh_rsq))],...
    ['S_2, R²=',num2str(mean(q_s_dh_rsq_g(:,2))),' \pm ',num2str(std(q_s_dh_rsq_g(:,2)))],...
    ['S_1, R²=',num2str(mean(q_s_dh_rsq_g(:,3))),' \pm ',num2str(std(q_s_dh_rsq_g(:,3)))]},'location','northwest')


%
subplot(1,2,2)
errorbar(receptors(group_var==1),mean_dh(group_var==1),std_dh(group_var==1),'b.','markersize',msize);hold on 
errorbar(receptors(group_var==0),mean_dh(group_var==0),std_dh(group_var==0),'c.','markersize',msize);hold on 
errorbar(receptors(group_var==-1),mean_dh(group_var==-1),std_dh(group_var==-1),'r.','markersize',msize);hold on 
% Plotting lines
p1=plot(r_xx,r_yy,'k-'); % whole data
p2=plot(r_xx,r_yy_m_1,'r-'); % r-group
xlim([0,1])
grid on
axis square
xlabel('Receptor Density')
ylabel('\Delta h')
ylim([0 0.1])

legend([p1 p2],{['All, R²=',num2str(mean(l_r_dh_rsq)),' \pm ',num2str(std(l_r_dh_rsq))],...
    ['R_1, R²=',num2str(mean(l_r_dh_rsq_g(:,1))),' \pm ',num2str(std(l_r_dh_rsq_g(:,1)))]},...
    'location','northwest')
%
print(gcf,'-dpdf',[figfold,'lme.pdf'],'-r300')


%% Printing groups
g_reg_names=arrayfun(@(x) dmf_labels(group_var==x),[1 0 -1],'uni',0);
%% Plotting brains with regions: test with one side
load('aal_node_pos.mat')
sort_ids = cat(2,[1:2:N],fliplr([2:2:N]));
node_pert = group_var;
node_pert_aal = node_pert;
node_pert_aal(sort_ids) = node_pert;
figure;
method  = 'euclidean';
mcolors = (othercolor('RdBu11',N));
meshclims = [-0.9 0.9];
nodecols=cat(1,[0 0 0.7],[0.5 0.5 0.5],[0.7 0 0]);
aa=atemplate('overlay',node_pert_aal,'meshcols',mcolors,'method',method,...
    'hemi','both','nocolbar','meshclims',meshclims);%,...
brighten(-.6)

% set views:
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
        atemplate('overlay',node_pert_aal,'meshcols',mcolors,'method',method,...
            'hemi',hemis{h},'nocolbar','meshclims',meshclims);  
        brighten(-.6)
        view(gca,v{i});        
        cl =camlight('headlight');
        set(cl,'Color',[0.4 0.4 0.4])
        print(gcf,'-dpng',[figfold,'brain_three_way_sep_',hemis{h},'_',zview{i},'.png'],'-r300')        
        close gcf
        
    end
    
    
end

%%
msize = 50;

% subplot(1,3,1)
% plot(sel_con(group_var==1),mean_ent(group_var==1,1),'b.','markersize',msize);hold on 
% plot(sel_con(group_var==-1),mean_ent(group_var==-1,1),'r.','markersize',msize);hold on 
% plot(sel_con(group_var==0),mean_ent(group_var==0,1),'k.','markersize',msize/2);hold on 
% grid on
% axis square

subplot(1,3,1)
plot(mean_ent(group_var==1,1),mean_ent(group_var==1,2),'b.','markersize',msize);hold on 
plot(mean_ent(group_var==-1,1),mean_ent(group_var==-1,2),'r.','markersize',msize);hold on 
plot(mean_ent(group_var==0,1),mean_ent(group_var==0,2),'k.','markersize',msize/2);hold on 
refline(1,0)
grid on
axis square

subplot(1,3,2)
plot(sel_con(group_var==1),mean_dh(group_var==1),'b.','markersize',msize);hold on 
plot(sel_con(group_var==-1),mean_dh(group_var==-1),'r.','markersize',msize);hold on 
plot(sel_con(group_var==0),mean_dh(group_var==0),'k.','markersize',msize/2);hold on 
plot(sel_con,predict(lme_fit,sel_tab),'g.','markersize',msize/3);hold on 

grid on
axis square

subplot(1,3,3)
% plot(receptors,predict(lme_fit,sel_tab),'g.');hold on 
plot(receptors(group_var==1),mean_dh(group_var==1),'b.','markersize',msize);hold on 
plot(receptors(group_var==-1),mean_dh(group_var==-1),'r.','markersize',msize);hold on 
plot(receptors(group_var==0),mean_dh(group_var==0),'k.','markersize',msize/2);hold on 
plot(receptors,predict(lme_fit,sel_tab),'g.','markersize',msize/3);hold on 
grid on
axis square
%

%% Using homotopic averages
homo_d_h = mean((deltah(1:45,:) + flipud(deltah(46:end,:)))./2,2);
homo_sel_con = (sel_con(1:45) + flipud(sel_con(46:end)))./2;
homo_rec = (receptors(1:45) + flipud(receptors(46:end)))./2;

figure
subplot(121)
plot(homo_sel_con,homo_d_h,'o');

subplot(122)
plot(homo_rec,homo_d_h,'o');

