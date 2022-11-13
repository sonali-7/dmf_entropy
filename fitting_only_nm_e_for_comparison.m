%% Running DMF with and without neuromodulation
% basefold = '/media/ruben/ssd240/Matlab/cb-neuromod-master/';
basefold = 'E:\Matlab\cb-neuromod-master\';
PARAMETERS_DIR=fullfile(dotenv.read().PROJECT_DIR, 'parameters')
load(fullfile(PARAMETERS_DIR,'fc_fcd_bold_sig_pcb_lsd.mat'),'fcd','fc','tr','flp','fhi',...
    'wsize','overlap','condnames','sel_conds')
load(fullfile(PARAMETERS_DIR,'SC_and_5ht2a_receptors.mat'))
C = sc90./max(sc90(:))*0.2;
stren = sum(C)./2;
nsubs = size(fc,3);
N = length(C);
nwins = size(fcd,1);
nconds=2;
dmf_labels = load('Structural.mat','labelsymm');
dmf_labels = arrayfun(@(x) string(dmf_labels.labelsymm(x,:)),1:N,'uni',0)';
%% Computing corr(SC,FC) per subject
% KS-Dist between all pairs of subjects
% FC corr between all pairs of subject
ave_fc = squeeze(mean(fc,3));
vec_sc = squareform(C);
emp_sc_fc_corr = zeros(nsubs,nconds);
%

fcd_ks_per_sub = cell(nconds,1);
fc_corr_per_sub = fcd_ks_per_sub;

fcd_ks_2_all = zeros(nsubs,nconds);
fc_corr_2_all = fcd_ks_2_all;

for c=1:nconds
    this_fcd = fcd(:,:,:,c);
    aux_fcd_ks_per_sub = zeros(nsubs,nsubs);
    aux_fc_corr_per_sub = aux_fcd_ks_per_sub;
    vec_fc = squareform(ave_fc(:,:,c) - eye(N));
    for s=1:nsubs

        aux_fc = squareform(fc(:,:,s,c) - eye(N));
        emp_sc_fc_corr(s,c) = corr2(aux_fc,vec_sc);
        aux_fcd = squareform(this_fcd(:,:,s) - eye(nwins));

        [~,~,fcd_ks_2_all(s,c)] = kstest2(reshape(this_fcd(:,:,s),[numel(this_fcd(:,:,s)),1]),this_fcd(:));
%         fc_corr_2_all(s,c) = corr2(fc(:,:,s,c),ave_fc(:,:,c));
        fc_corr_2_all(s,c) = mean((aux_fc-vec_fc).^2); % MSE FC


        for s2=s:nsubs
            aux_fc2 = squareform(fc(:,:,s2,c) - eye(N));
%             aux_fc_corr_per_sub(s,s2) = corr2(aux_fc(:),aux_fc2(:));
            aux_fc_corr_per_sub(s,s2) = mean((aux_fc(:)-aux_fc2(:)).^2); % MSE FC

            aux_fcd2 = squareform(this_fcd(:,:,s2) - eye(nwins));
            [~,~,aux_fcd_ks_per_sub(s,s2)] = kstest2(aux_fcd(:),aux_fcd2(:));
        end

    end
    isub_sub = find(triu(ones(nsubs),1));
    fc_corr_per_sub{c} = aux_fc_corr_per_sub(isub_sub);
    fcd_ks_per_sub{c} = aux_fcd_ks_per_sub(isub_sub);
end

%% Preparing parameters
% dmf parameters
[ params ] = dyn_fic_DefaultParams('C',C);
params.burnout = 10;
params.flp = flp;
params.fhi = fhi;
params.wsize = wsize;
params.overlap = overlap;
params.TR = tr;
params.batch_size = 50000;
params.receptors = receptors;
% no optimization of dynamic FIC
params.lrj = 0;
params.taoj = Inf;

%% Running optimizer: Optimizing only nm_e with fix G
% Bayes optimization options
sel_c = 2;
savefold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/';
savename = [savefold,'bayesopt_dmf_checkpoint_G_only_nm_e_',condnames{sel_c},'_G2.4_v2.mat'];

bo_opts = {'IsObjectiveDeterministic',false,'UseParallel',true,...
        'MinWorkerUtilization',10,...
        'AcquisitionFunctionName','expected-improvement-plus',...
        'MaxObjectiveEvaluations',1e16,...
        'ParallelMethod','clipped-model-prediction',...
        'GPActiveSetSize',300,'ExplorationRatio',0.5,'MaxTime',3600*7*24,...
        'OutputFcn',@saveToFile,...
        'SaveFileName',savename};


%% Optimizable parameters
G = 2.42;
nm_e = [0 0.1];
nm_i = 0;

T = 510;
[opt_fc_error,opt_fcd_ks,opt_pars,bayesopt_out] = fit_fc_fcd_dmf_neuromod(...
    T,ave_fc(:,:,sel_c),fcd(:,:,:,sel_c),G,nm_e,nm_i,params,bo_opts);

%% Checking optimization results
% load([savefold,'bayesopt_dmf_checkpoint_',condnames{1},'.mat'])
load(savename)
[best_pars,est_min_ks] = bestPoint(BayesoptResults,'Criterion','min-mean')
%% Extracting firing rates and data to plot objective function and others
% Firing rate and entropy
par_names = {'G','alpha'};
sel_id = 702;
gvals = BayesoptResults.XTrace.G(sel_id:end);
alphavals = BayesoptResults.XTrace.alpha(sel_id:end);
ksvals = BayesoptResults.ObjectiveTrace(sel_id:end);
data_trace =BayesoptResults.UserDataTrace(sel_id:end);
mean_fr = cellfun(@(x) mean(x{1}),data_trace);
mean_ent= cellfun(@(x) mean(x{2}),data_trace);
% Fitting to rates
fit3d_frs = fit([gvals,alphavals],mean_fr,'loess','Span',0.5);
% Fitting to entropie
fit3d_ent = fit([gvals,alphavals],mean_ent,'loess','Span',0.5);

nxx = 100;
gvec = linspace(G(1),G(2),nxx);
alphavec = linspace(fic_alpha(1),fic_alpha(2),nxx);
[xx,yy] = meshgrid(gvec,alphavec);

% extracts rates
frs_vs_pars = feval(fit3d_frs,[xx(:),yy(:)]);
this_frs= reshape(frs_vs_pars,[nxx,nxx]);

% extracts entropies
ent_vs_pars = feval(fit3d_ent,[xx(:),yy(:)]);
this_ent= reshape(ent_vs_pars,[nxx,nxx]);

Xtable = table(xx(:),yy(:),'VariableNames',par_names);
[this_objective,this_sigma] = predictObjective(BayesoptResults,Xtable);
this_objective = reshape(this_objective,[nxx,nxx]);


%%
cmap = flipud(othercolor('YlGnBu5',256));
figfold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/figures/';
figname = [condnames{sel_c},'_bayesopt_objfunc_G_and_alpha_v4'];
kslims = [0 0.7];
kslevels = linspace(kslims(1),kslims(2),11);
frlims = [0 10];
frlevels = linspace(frlims(1),frlims(2),11);
entlims = [1.8 2.8];
entlevels = linspace(entlims(1),entlims(2),11);

figure('units','normalized','outerposition',[0 0 1 1],'paperpositionmode','auto')
subplot(221)
% imagesc(gvec,alphavec,this_objective,kslims);hold on
contourf(gvec,alphavec,this_objective,kslevels,'showtext','on');hold on
set(gca,'ydir','normal')
pp=plot(best_pars.G(1),best_pars.alpha(1),'r*','markersize',15);
legend(pp,['G = ',num2str(best_pars.G(1)),...
    ', \alpha = ',num2str(best_pars.alpha(1))],'fontsize',11)
grid on
cb = colorbar;
cb.Label.String = 'mean K-S (FCD_{emp},FCD_{dmf})';
colormap(flipud(cmap))
ylabel('\alpha')
xlabel('G')
axis square
title('KS FCD')

subplot(223)
[cc,hh] = contourf(gvec,alphavec,this_frs,frlevels,'showtext','on');hold on
s=contourdata(cc);
set(gca,'ydir','normal')
pp=plot(best_pars.G(1),best_pars.alpha(1),'r*','markersize',15);
grid on
cb = colorbar;
cb.Label.String = 'mean E Firing Rate (Hz))';
colormap(flipud(cmap))
caxis(frlims)
ylabel('\alpha')
xlabel('G')
axis square
title('Brain Average E Firing Rate (Hz)')

subplot(222)
% [cc,hh] = contourf(gvec,alphavec,this_objective,kslevels,'showtext','on');hold on
imagesc(gvec,alphavec,this_objective,kslims);hold on
for l=1:length(s)
    if s(l).level==3
        pp=plot(s(l).xdata,s(l).ydata,'color',[0 0.4 0.2],'linewidth',2);
    else
        plot(s(l).xdata,s(l).ydata,'k');
    end
end
set(gca,'ydir','normal')
plot(best_pars.G(1),best_pars.alpha(1),'r*','markersize',15);
grid on
cb = colorbar;
cb.Label.String = 'mean K-S (FCD_{emp},FCD_{dmf})';
colormap(flipud(cmap))
caxis(kslims)
ylabel('\alpha')
xlabel('G')
axis square
title('KS FCD')
legend(pp,'3 Hz')


subplot(224)
contourf(gvec,alphavec,this_ent,entlevels,'showtext','on');hold on
set(gca,'ydir','normal')
pp=plot(best_pars.G(1),best_pars.alpha(1),'r*','markersize',15);
grid on
cb = colorbar;
cb.Label.String = 'mean E Entropy (nats))';
colormap(flipud(cmap))
caxis(entlims)
ylabel('\alpha')
xlabel('G')
axis square
title('Brain Average E Entropy (Hz)')

%
print(gcf,'-dpng',[figfold,figname,'.png'],'-r300')
print(gcf,'-dpdf',[figfold,figname,'.pdf'],'-r300')
%% Running simulation at minimum to check ks vals for PCB
gamma_ent_fun = @(a) a(1) + log(a(2)) + log(gamma(a(1))) + (1-a(1))*psi(a(1));
nreps = 60;
selG = 2.42;
sel_alpha = 1.49;
selpars = params;
stren = sum(selpars .C)./2;
N = length(stren);
isubfc = find(tril(ones(N),-1));

selpars.G = selG;
selpars.J = sel_alpha*selpars.G*stren' + 1; % updates it
selpars.batch_size = 10000;
selpars.wgaine = 0.0;
selpars.wgaini = 0.0;

sel_ks_fcd = zeros(nreps,nconds);
sel_fc_mse = sel_ks_fcd;
reg_fr = zeros(N,nreps);
reg_ent = reg_fr;
bold_sigs = cell(nreps,1);

% T = 110;
T = 510;
nsteps = T.*(1000); % number of DMF timepoints
%
init1 = tic;
parfor r=1:nreps
    % Simulating
    initic=tic;
    [rates,bold] = dyn_fic_DMF(selpars, nsteps,'both'); % runs simulation
    rates = rates(:,(selpars.burnout*1000*2):end);
    reg_fr(:,r) = mean(rates,2);
    for n=1:N
        gamma_pars = gamfit(rates(n,:));
        reg_ent(n,r) = gamma_ent_fun(gamma_pars);
    end
    bold = bold(:,selpars.burnout:end); % remove initial transient
    bold(isnan(bold))=0;
    bold(isinf(bold(:)))=max(bold(~isinf(bold(:))));
    % Filtering and computing FC
    filt_bold = filter_bold(bold',selpars.flp,selpars.fhi,selpars.TR);
    sim_fc = corrcoef(filt_bold);
    % FCD
    sim_fcd = compute_fcd(filt_bold,selpars.wsize,selpars.overlap,isubfc);
    sim_fcd(isnan(sim_fcd))=0;
    sim_fcd = corrcoef(sim_fcd);

    aux_ks_fcd = zeros(nconds,1);

    for c=1:nconds
        this_fc = ave_fc(:,:,c);
        this_fcd= fcd(:,:,:,c);
        % Computing FC error: 1-Corrrelation between FC's
        sel_fc_mse(r,c) = mean((sim_fc(isubfc)-this_fc(isubfc)).^2); % MSE FC

        % Computing KS FCD for each condition and MSE on FC
        [~,~,aux_ks_fcd(c)] = kstest2(sim_fcd(:),this_fcd(:));
    end
    bold_sigs{r} = filt_bold;
    sel_ks_fcd(r,:) = aux_ks_fcd;

    endtime=toc(initic);
    disp(['Simulation ',num2str(r),', PCB K-S = ',num2str(aux_ks_fcd(1)),...
        ', LSD K-S = ',num2str(aux_ks_fcd(2)),'. Time = ',num2str(endtime),' seconds'])


end
toc(init1)
%
mean_fr = mean(reg_fr,2);
mean_ent = mean(reg_ent,2);

%% Running optimizer to optimize both E and I neuromodulation
sel_c = 2;
savefold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/';
savename = [savefold,'bayesopt_dmf_checkpoint_nm_e_nm_i_',condnames{sel_c},'.mat'];
bo_opts = {'IsObjectiveDeterministic',false,'UseParallel',true,...
        'MinWorkerUtilization',30,...
        'AcquisitionFunctionName','expected-improvement-plus',...
        'MaxObjectiveEvaluations',1e16,...
        'ParallelMethod','clipped-model-prediction',...
        'GPActiveSetSize',300,'ExplorationRatio',0.5,'MaxTime',3600*7*24   ,...
        'OutputFcn',@saveToFile,...
        'SaveFileName',savename};


% Optimizable parametersse
thispars  = params;
selG = 2.4;
sel_alpha = 1.5;

thispars.G = selG;
thispars.J = sel_alpha*thispars.G*stren' + 1; %
G = [];
nm_e = [0 0.2];
nm_i = [0 0.2];
% nm_i = [];


T = 510;
[opt_fc_error,opt_fcd_ks,opt_pars,bayesopt_out] = fit_fc_fcd_dmf_neuromod(...
    T,ave_fc(:,:,sel_c),fcd(:,:,:,sel_c),G,nm_e,nm_i,thispars,bo_opts);
%% Checking optimization results
% load(savename)
[best_pars,est_min_ks] = bestPoint(BayesoptResults,'Criterion','min-mean')
%% Checking objective function
par_names = {'nm_e','nm_i'};
nxx = 1000;
nm_e = linspace(0,0.25,nxx);
nm_i = linspace(0,0.25,nxx);
[xx,yy] = meshgrid(nm_e,nm_i);
Xtable = table(xx(:),yy(:),'VariableNames',par_names);
[this_objective,this_sigma] = predictObjective(BayesoptResults,Xtable);
this_objective = reshape(this_objective,[nxx,nxx]);

cmap = flipud(othercolor('YlGnBu5',256));
%% Selected points
sel_nm = 0.024;
Xtable = table(sel_nm ,sel_nm ,'VariableNames',par_names);
[sel_objective,sel_sigma] = predictObjective(BayesoptResults,Xtable)

%% Plotting objective function
selG = 2.4;
sel_alpha = 1.5;
% figfold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/figures/';
figfold = 'E:\Matlab\fastdmf-master\newSciRep\';
figname = ['bayeopt_objfunc_nm_e_nm_i_G_',num2str(selG),...
    '_alpha_',num2str(sel_alpha)];

kslims = [0 0.7];
kslevels = linspace(kslims(1),kslims(2),11);

figure('units','normalized','outerposition',[0 0 1 1],'paperpositionmode','auto')
subplot(121)
imagesc(nm_e,nm_i,this_objective,kslims);hold on
set(gca,'ydir','normal')
% plot(best_pars.nm_e(1),best_pars.nm_i(1),'r*','markersize',15)
% legend(['nm_e = ',num2str(best_pars.nm_e(1)),...
%     ', nm_i = ',num2str(best_pars.nm_i(1))],'fontsize',11)
plot(sel_nm,sel_nm,'r*','markersize',15);
legend(['nm_e = ',num2str(sel_nm),...
    ', nm_i = ',num2str(sel_nm)],'fontsize',11)

grid on


cb = colorbar;
cb.Label.String = 'mean K-S (FCD_{emp},FCD_{dmf})';
title(['G=',num2str(selG),', \alpha_{fic} =',num2str(sel_alpha)])
colormap(flipud(cmap))
ylabel('nm_i')
xlabel('nm_e')
axis square

subplot(122)
contourf(nm_e,nm_i,this_objective,kslevels,'showtext','on');hold on
set(gca,'ydir','normal')
% pp=plot(best_pars.nm_e(1),best_pars.nm_i(1),'r*','markersize',15);
% legend(pp,['nm_e = ',num2str(best_pars.nm_e(1)),...
%     ', nm_i = ',num2str(best_pars.nm_i(1))],'fontsize',11)
pp=plot(sel_nm,sel_nm,'r*','markersize',15);
legend(pp,['nm_e = ',num2str(sel_nm),...
    ', nm_i = ',num2str(sel_nm)],'fontsize',11)
grid on
cb = colorbar;
cb.Label.String = 'mean K-S (FCD_{emp},FCD_{dmf})';
colormap(flipud(cmap))
caxis(kslims)
ylabel('nm_i')
xlabel('nm_e')
axis square
title(['K-S = ',num2str(sel_objective),' \pm ',num2str(sel_sigma)])


%
print(gcf,'-dpng',[figfold,figname,'.png'],'-r300')
print(gcf,'-dpdf',[figfold,figname,'.pdf'],'-r300','-painters')

%% Running model for both conditions with optimal parameters
selG = 2.4;
sel_alpha = 1.5;
sel_wgaine = [0 0.024];
sel_wgaini = [0 0.024];

% Default pars
[ params ] = dyn_fic_DefaultParams('C',C);
params.burnout = 10;
params.flp = flp;
params.fhi = fhi;
params.wsize = wsize;
params.overlap = overlap;
params.TR = tr;
params.batch_size = 50000;
params.receptors = receptors;
params.lrj = 0;
params.taoj = Inf;
params.G = selG;
params.J = sel_alpha*params.G*stren' + 1; % updates

% Creating parameters for simulations
nmods = 2;
nreps = 300;
iniconds = randperm(1000,nreps);
parlist = cell(nmods,nreps);
for m=1:nmods
    thispars = params;
    thispars.wgaine = sel_wgaine(m);
    thispars.wgaini = sel_wgaini(m);
    for r=1:nreps
        thispars.seed = iniconds(r);
        parlist{m,r} = thispars;
    end
end
%

gamma_ent_fun = @(a) a(1) + log(a(2)) + log(gamma(a(1))) + (1-a(1))*psi(a(1));

stren = sum(params .C)./2;
N = length(stren);
isubfc = find(tril(ones(N),-1));

sel_ks_fcd = zeros(nreps,nconds,nmods);
sel_fc_mse = sel_ks_fcd;
reg_fr = zeros(N,nreps,nmods);
reg_ent = reg_fr;
bold_sigs = cell(nreps,nmods);

% T = 110;
T = 510;
nsteps = T.*(1000); % number of DMF timepoints
%
init1 = tic;
for m=1:nmods
    parfor r=1:nreps
        selpars = parlist{m,r};
        % Simulating
        initic=tic;
        [rates,bold] = dyn_fic_DMF(selpars, nsteps,'both'); % runs simulation
        rates = rates(:,(selpars.burnout*1000*2):end);
        reg_fr(:,r,m) = mean(rates,2);
        for n=1:N
            gamma_pars = gamfit(rates(n,:));
            reg_ent(n,r,m) = gamma_ent_fun(gamma_pars);
        end
        bold = bold(:,selpars.burnout:end); % remove initial transient
        bold(isnan(bold))=0;
        bold(isinf(bold(:)))=max(bold(~isinf(bold(:))));
        % Filtering and computing FC
        filt_bold = filter_bold(bold',selpars.flp,selpars.fhi,selpars.TR);
        sim_fc = corrcoef(filt_bold);
        % FCD
        sim_fcd = compute_fcd(filt_bold,selpars.wsize,selpars.overlap,isubfc);
        sim_fcd(isnan(sim_fcd))=0;
        sim_fcd = corrcoef(sim_fcd);

        aux_ks_fcd = zeros(nconds,1);

        for c=1:nconds
            this_fc = ave_fc(:,:,c);
            this_fcd= fcd(:,:,:,c);
            % Computing FC error: 1-Corrrelation between FC's
            sel_fc_mse(r,c,m) = mean((sim_fc(isubfc)-this_fc(isubfc)).^2); % MSE FC

            % Computing KS FCD for each condition and MSE on FC
            [~,~,aux_ks_fcd(c)] = kstest2(sim_fcd(:),this_fcd(:));
        end
        bold_sigs{r,m} = filt_bold;
        sel_ks_fcd(r,:,m) = aux_ks_fcd;

        endtime=toc(initic);
%         disp(['model = ',condnames{m}])
%         disp(['Simulation ',num2str(r),', PCB K-S = ',num2str(aux_ks_fcd(1)),...
%             ', LSD K-S = ',num2str(aux_ks_fcd(2)),'. Time = ',num2str(endtime),' seconds'])


    end
end
toc(init1)
%%
save([savefold,'dmf_pcb_lsd_nm_i_v5.mat'],'selG','sel_alpha','sel_wgaini','sel_wgaine','iniconds',...
    'sel_ks_fcd','sel_fc_mse','reg_fr','reg_ent','bold_sigs','params','nreps');

%% Preparing data to plot
entlims = [1.8 2.8];
frlims = [0 10];
condnames = {'pcb','lsd'};
mean_fr = squeeze(mean(reg_fr,2));
mean_ent= squeeze(mean(reg_ent,2));
mean_ks_fcd = squeeze(mean(sel_ks_fcd));
mean_mse_fc = squeeze(mean(sel_fc_mse));
% Plotting Firing Rates, Entropy, and FCD dist
% (condition,model)
figname =['ks_optimal_models_pcb_lsd_v4'];
figure('units','normalized','outerposition',[0 0 1 1],'paperpositionmode','auto')
subplot(241)
plot([1 2],[sel_ks_fcd(:,1,1) sel_ks_fcd(:,2,2)],'o-','color',[0.5 0.5 0.5]);hold on
plot([1 2],[mean_ks_fcd(1,1) mean_ks_fcd(2,2)],'ko-','linewidth',2);hold on
xlim([0.8 2.2])
set(gca,'xtick',[1 2],'xticklabel',condnames)
ylabel('KS FCD')
title('Using Respective Optimal Model')
set(gca,'ygrid','on')
ylim([0 0.5]);
axis square

subplot(242)
plot([1 2],[sel_ks_fcd(:,1,1) sel_ks_fcd(:,2,1)],'o-','color',[0.5 0.5 0.5]);hold on
plot([1 2],[mean_ks_fcd(1,1) mean_ks_fcd(2,1)],'ko-','linewidth',2);hold on
xlim([0.8 2.2])
set(gca,'xtick',[1 2],'xticklabel',condnames)
ylabel('KS FCD')
title('Using PCB Optimal Model')
set(gca,'ygrid','on')
ylim([0 0.5]);
axis square

subplot(243)
plot([1 2],[sel_ks_fcd(:,1,2) sel_ks_fcd(:,2,2)],'o-','color',[0.5 0.5 0.5]);hold on
plot([1 2],[mean_ks_fcd(1,2) mean_ks_fcd(2,2)],'ko-','linewidth',2);hold on
xlim([0.8 2.2])
set(gca,'xtick',[1 2],'xticklabel',condnames)
ylabel('KS FCD')
title('Using LSD Optimal Model')
set(gca,'ygrid','on')
ylim([0 0.5]);
axis square

deltah= (reg_ent(:,:,2)-reg_ent(:,:,1))./reg_ent(:,:,1);

subplot(245)
plot(mean_ent(:,1),mean_ent(:,2),'.','markersize',10);hold on
plot(entlims,entlims,'k--')
xlabel('Region Entropy PCB (nats)')
ylabel('Region Entropy 5HT2A (nats)')
grid on
axis square
xlim([2 2.6])
ylim([2 2.6])
title(['nm_e=',num2str(sel_wgaine(2)),', nm_i=',num2str(sel_wgaini(2))])

subplot(246)
plot(mean_fr(:,1),mean_fr(:,2),'.','markersize',10);hold on
plot(frlims,frlims,'k--')
xlabel('Region Firing Rate PCB (Hz)')
ylabel('Region Firing Rate 5HT2A (Hz)')
grid on
axis square
xlim([3 4])
ylim([3 4])

subplot(247)
plot(stren,mean(deltah,2),'.','markersize',10);hold on
xlabel('Region Strength')
ylabel('Region \Delta Entropy')
grid on
axis square

subplot(248)
plot(receptors,mean(deltah,2),'.','markersize',10);hold on
xlabel('Region 5-HT2AR Density')
ylabel('Region \Delta Entropy')
grid on
axis square

%
print(gcf,'-dpng',[figfold,figname,'.png'],'-r300')
print(gcf,'-dpdf',[figfold,figname,'.pdf'],'-r300','-painters')

%% Running optimizer with Fix G and Alpha and neuromodulating FIC
sel_c = 2;
savefold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/';
savename = [savefold,'bayesopt_dmf_checkpoint_fix_G_alpha_fic_nm_',condnames{sel_c},'.mat'];
bo_opts = {'IsObjectiveDeterministic',false,'UseParallel',true,...
        'MinWorkerUtilization',20,...
        'AcquisitionFunctionName','expected-improvement-plus',...
        'MaxObjectiveEvaluations',1e16,...
        'ParallelMethod','clipped-model-prediction',...
        'GPActiveSetSize',300,'ExplorationRatio',0.6,'MaxTime',3600*7*24   ,...
        'OutputFcn',@saveToFile,...
        'SaveFileName',savename};


% Optimizable parameters
sel_alpha = 1.49;
thispars = params;
thispars.G = 2.42;
thispars.J = sel_alpha*thispars.G*stren' + 1; %

G = [];
fic_alpha = [];
nm = [-0.2 0.2];


T = 510;
[opt_fc_error,opt_fcd_ks,opt_pars,bayesopt_out] = fit_fc_fcd_dmf_nm_fic(...
    T,ave_fc(:,:,sel_c),fcd(:,:,:,sel_c),G,fic_alpha,nm,thispars,bo_opts);

%% Checking optimization results
load(savename)
[best_pars,est_min_ks] = bestPoint(BayesoptResults,'Criterion','min-mean')

%%
gamma_ent_fun = @(a) a(1) + log(a(2)) + log(gamma(a(1))) + (1-a(1))*psi(a(1));
nreps = 60;
selG = 2.42;
sel_alpha = 1.49;

selpars = params;
stren = sum(selpars .C)./2;
N = length(stren);
isubfc = find(tril(ones(N),-1));

selpars.G = selG;
selpars.J = sel_alpha*selpars.G*stren' + 1; % updates it
selpars.batch_size = 10000;
selpars.wgaine = 0.02;
selpars.wgaini = 0.02;

sel_ks_fcd = zeros(nreps,nconds);
sel_fc_mse = sel_ks_fcd;
reg_fr = zeros(N,nreps);
reg_ent = reg_fr;
bold_sigs = cell(nreps,1);

% T = 110;
T = 510;
nsteps = T.*(1000); % number of DMF timepoints
%
init1 = tic;
parfor r=1:nreps
    % Simulating
    initic=tic;
    [rates,bold] = dyn_fic_DMF(selpars, nsteps,'both'); % runs simulation
    rates = rates(:,(selpars.burnout*1000*2):end);
    reg_fr(:,r) = mean(rates,2);
    for n=1:N
        gamma_pars = gamfit(rates(n,:));
        reg_ent(n,r) = gamma_ent_fun(gamma_pars);
    end
    bold = bold(:,selpars.burnout:end); % remove initial transient
    bold(isnan(bold))=0;
    bold(isinf(bold(:)))=max(bold(~isinf(bold(:))));
    % Filtering and computing FC
    filt_bold = filter_bold(bold',selpars.flp,selpars.fhi,selpars.TR);
    sim_fc = corrcoef(filt_bold);
    % FCD
    sim_fcd = compute_fcd(filt_bold,selpars.wsize,selpars.overlap,isubfc);
    sim_fcd(isnan(sim_fcd))=0;
    sim_fcd = corrcoef(sim_fcd);

    aux_ks_fcd = zeros(nconds,1);

    for c=1:nconds
        this_fc = ave_fc(:,:,c);
        this_fcd= fcd(:,:,:,c);
        % Computing FC error: 1-Corrrelation between FC's
        sel_fc_mse(r,c) = mean((sim_fc(isubfc)-this_fc(isubfc)).^2); % MSE FC

        % Computing KS FCD for each condition and MSE on FC
        [~,~,aux_ks_fcd(c)] = kstest2(sim_fcd(:),this_fcd(:));
    end
    bold_sigs{r} = filt_bold;
    sel_ks_fcd(r,:) = aux_ks_fcd;

    endtime=toc(initic);
    disp(['Simulation ',num2str(r),', PCB K-S = ',num2str(aux_ks_fcd(1)),...
        ', LSD K-S = ',num2str(aux_ks_fcd(2)),'. Time = ',num2str(endtime),' seconds'])


end
toc(init1)
%
mean_fr = mean(reg_fr,2);
mean_ent = mean(reg_ent,2);

figure
subplot(121)
plot(mean_fr)

subplot(122)
plot(mean_ent)



%% Running optimizer with FIX G and neuromodulation
% Bayes optimization options
% best parameters is 0.026 for lsd
savefold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/';
prev_opt = load([savefold,'bayesopt_dmf_checkpoint_nm_linked_',condnames{2},'.mat']);
sel_c = 2;
bo_opts = {'IsObjectiveDeterministic',false,'UseParallel',true,...
        'MinWorkerUtilization',16,...
        'AcquisitionFunctionName','expected-improvement-plus',...
        'MaxObjectiveEvaluations',1e16,...
        'ParallelMethod','clipped-model-prediction',...
        'GPActiveSetSize',300,'ExplorationRatio',0.3,'MaxTime',3600*7*24   ,...
        'OutputFcn',@saveToFile,...
        'SaveFileName',[savefold,'bayesopt_dmf_checkpoint_nm_linked_',condnames{sel_c},'_v2.mat'],...
        'InitialX',prev_opt.BayesoptResults.XTrace,...
        'InitialObjective',prev_opt.BayesoptResults.ObjectiveTrace,...
        'InitialErrorValues',prev_opt.BayesoptResults.ErrorTrace};


% Optimizable parametersse
selG = 2.4;
sel_alpha=1.5;
selpars=params;
selpars.G = selG;
selpars.J = sel_alpha*selpars.G*stren' + 1; % updates it
G = [];
nm_e = [0 0.2];
% nm_i = [0 0.4];
nm_i = [];


T = 510;
[opt_fc_error,opt_fcd_ks,opt_pars,bayesopt_out] = fit_fc_fcd_dmf_neuromod(...
    T,ave_fc(:,:,sel_c),fcd(:,:,:,sel_c),G,nm_e,nm_i,selpars,bo_opts);

%% Checking optimization results
% load([savefold,'bayesopt_dmf_checkpoint_',condnames{sel_c},'.mat'])
load([savefold,'bayesopt_dmf_checkpoint_nm_linked_',condnames{sel_c},'_v2.mat'])
[best_pars,est_min_ks] = bestPoint(BayesoptResults,'Criterion','min-mean')
%% Plotting objective function
xx = linspace(0,0.4,100);
Xtable = table(xx(:),'VariableNames',{'nm_e'});
[this_objective,this_sigma] = predictObjective(BayesoptResults,Xtable);

figure
shadedErrorBar(xx,this_objective,this_sigma)

%%
figfold = [savefold,'figures'];
figure
figname = ['bayesopt_iterations_pcb_music_only_G'];
plot(BayesoptResults.EstimatedObjectiveMinimumTrace);hold on
plot(BayesoptResults.ObjectiveMinimumTrace)
% plot(est_bayesopt_trace(sel_ids));hold on
% plot(obs_bayesopt_trace(sel_ids))
xlim([1 350])
axis square
xlabel('Iterations')
ylabel('KS(FCDemp,FCDsim)')
legend({'Estimated Minimum','Observed Minimum'})
grid on
%
print(gcf,'-dpng',[figfold,figname,'.png'],'-r300')
print(gcf,'-dpdf',[figfold,figname,'.pdf'],'-r300')


%% Running simulation at minimum to check ks vals for PCB
gamma_ent_fun = @(a) a(1) + log(a(2)) + log(gamma(a(1))) + (1-a(1))*psi(a(1));
nreps = 90;
selG = 2.5;
selpars = params;
stren = sum(selpars .C)./2;
N = length(stren);
isubfc = find(tril(ones(N),-1));

selpars.G = selG;
selpars.J = 1.5*selpars.G*stren' + 1; % updates it
selpars.batch_size = 10000;
selpars.wgaine = 0.065;
selpars.wgaini = 0.065;

sel_ks_fcd = zeros(nreps,nconds);
sel_fc_mse = sel_ks_fcd;
reg_fr = zeros(N,nreps);
reg_ent = reg_fr;
bold_sigs = cell(nreps,1);

T = 510;
nsteps = T.*(1000); % number of DMF timepoints
%
init1 = tic;
parfor r=1:nreps
    % Simulating
    initic=tic;
    [rates,bold] = dyn_fic_DMF(selpars, nsteps,'both'); % runs simulation
    rates = rates(:,(selpars.burnout*1000*2):end);
    reg_fr(:,r) = mean(rates,2);
    for n=1:N
        gamma_pars = gamfit(rates(n,:));
        reg_ent(n,r) = gamma_ent_fun(gamma_pars);
    end
    bold = bold(:,selpars.burnout:end); % remove initial transient
    bold(isnan(bold))=0;
    bold(isinf(bold(:)))=max(bold(~isinf(bold(:))));
    % Filtering and computing FC
    filt_bold = filter_bold(bold',selpars.flp,selpars.fhi,selpars.TR);
    sim_fc = corrcoef(filt_bold);
    % FCD
    sim_fcd = compute_fcd(filt_bold,selpars.wsize,selpars.overlap,isubfc);
    sim_fcd(isnan(sim_fcd))=0;
    sim_fcd = corrcoef(sim_fcd);

    aux_ks_fcd = zeros(nconds,1);

    for c=1:nconds
        this_fc = ave_fc(:,:,c);
        this_fcd= fcd(:,:,:,c);
        % Computing FC error: 1-Corrrelation between FC's
        sel_fc_mse(r,c) = mean((sim_fc(isubfc)-this_fc(isubfc)).^2); % MSE FC

        % Computing KS FCD for each condition and MSE on FC
        [~,~,aux_ks_fcd(c)] = kstest2(sim_fcd(:),this_fcd(:));
    end
    bold_sigs{r} = filt_bold;
    sel_ks_fcd(r,:) = aux_ks_fcd;

    endtime=toc(initic);
    disp(['Simulation ',num2str(r),', PCB K-S = ',num2str(aux_ks_fcd(1)),...
        ', LSD K-S = ',num2str(aux_ks_fcd(2)),'. Time = ',num2str(endtime),' seconds'])


end
toc(init1)
%
% save([savefold,'checking_bayesopt_pcb_music_only_G.mat'],'bold_sigs','sel_ks_fcd',...
%     'sel_fc_corr','reg_fr','selpars');
save([savefold,'checking_bayesopt_lsd_music_fixG_nme_nmi_linked.mat'],'bold_sigs','sel_ks_fcd',...
    'sel_fc_corr','reg_fr','selpars');


%%
mean_fr = mean(reg_fr,2);
mean_ent = mean(reg_ent,2);

