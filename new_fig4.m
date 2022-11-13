%% New figure 4 for scientific reports: Connectome controls
load('SC_and_5ht2a_receptors.mat')
C = sc90./max(sc90(:))*0.2;
stren = sum(C)./2;
basefold = '/media/ruben/ssd240/Matlab/cb-neuromod-master/';
load([basefold,'fc_fcd_bold_sig_pcb_lsd.mat'],'fcd','fc','tr','flp','fhi',...
    'wsize','overlap','condnames','sel_conds')
ave_fc = squeeze(mean(fc,3));
vec_sc = squareform(C);

% Creating random connectivty matrix with the same connection density and link distribution
% Uses Brain Connectivity Toolbox
c_triu = squareform(C);
c_triu_rand = c_triu(randperm(length(c_triu)));
c_rand = squareform(c_triu_rand);
% Creating random matrix with preserved degree distribution
[c_dpr,eff]=randmio_und(C, 10);
% Creating random matrix that preservs degree and strength distribution
[c_dspr,R] = null_model_und_sign(C,0,0.1);
C_cond = cat(3,C,c_rand,c_dpr,c_dspr);
nnets = size(C_cond,3);

% DMF parameters
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
% params.J = sel_alpha.*stren.*params.G + 1;


% Creating parameters for simulations
nmods = 2;
nreps = 120;
iniconds = randperm(1000,nreps);
parlist = cell(nnets,nmods,nreps);
for c=1:nnets
    thisstren = sum(C_cond(:,:,c))./2;
    thispars = params;
    thispars.C = C_cond(:,:,c);
    thispars.J = sel_alpha.*thisstren.*thispars.G + 1;
    for m=1:nmods
        
        thispars.wgaine = sel_wgaine(m);
        thispars.wgaini = sel_wgaini(m);
        for r=1:nreps
            thispars.seed = iniconds(r);
            parlist{c,m,r} = thispars;
        end
    end
end

% Gamma entropy function
nconds=2;
gamma_ent_fun = @(a) a(1) + log(a(2)) + log(gamma(a(1))) + (1-a(1))*psi(a(1));

stren = sum(params.C)./2;
N = length(stren);
isubfc = find(tril(ones(N),-1));

sel_ks_fcd = zeros(nnets,nreps,nconds,nmods);
sel_fc_mse = sel_ks_fcd;
reg_fr = zeros(N,nreps,nnets,nmods);
reg_ent = reg_fr;
bold_sigs = cell(nreps,nnets,nmods);

T = 510;
% T = 110;
nsteps = T.*(1000); % number of DMF timepoints
%
init1 = tic;
for w=1:nnets
    for m=1:nmods
        initic=tic;
        thispars = parlist(w,m,:);
        parfor r=1:nreps
            selpars = thispars{r};
            % Simulating
            
            [rates,bold] = dyn_fic_DMF(selpars, nsteps,'both'); % runs simulation
            rates = rates(:,(selpars.burnout*1000*2):end);
            reg_fr(:,r,w,m) = mean(rates,2);
            for n=1:N
                gamma_pars = gamfit(rates(n,:));
                reg_ent(n,r,w,m) = gamma_ent_fun(gamma_pars);
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
                sel_fc_mse(w,r,c,m) = mean((sim_fc(isubfc)-this_fc(isubfc)).^2); % MSE FC
                
                % Computing KS FCD for each condition and MSE on FC
                [~,~,aux_ks_fcd(c)] = kstest2(sim_fcd(:),this_fcd(:));
            end
            bold_sigs{r,w,m} = filt_bold;
            sel_ks_fcd(w,r,:,m) = aux_ks_fcd;
            
        end
        toc(initic)
    end
end
toc(init1)

savefold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/';
save([savefold,'dmf_connectome_controls.mat'],'selG','sel_alpha','sel_wgaini','sel_wgaine','iniconds',...
    'sel_ks_fcd','sel_fc_mse','reg_fr','reg_ent','bold_sigs','params','nreps','C_cond');


%% Computing correlation to Human Delta H
deltah = squeeze((reg_ent(:,:,:,2)-reg_ent(:,:,:,1))./reg_ent(:,:,:,1));
mean_pcb = squeeze(mean(reg_ent(:,:,:,1),2));
mean_lsd= squeeze(mean(reg_ent(:,:,:,2),2));
r_sq_rnd = zeros(nreps,1);
r_sq_dpr = r_sq_rnd;
r_sq_dspr = r_sq_rnd;
sel_delta_ent = squeeze(deltah(:,:,1));
for r=1:nreps
    r_sq_rnd(r) = corr(sel_delta_ent(:,r),deltah(:,r,2));
    r_sq_dpr(r) = corr(sel_delta_ent(:,r),deltah(:,r,3));
    r_sq_dspr(r) = corr(sel_delta_ent(:,r),deltah(:,r,4));
end
r_sq_2_plot = cat(2,r_sq_rnd,r_sq_dpr,r_sq_dspr);

%% Plotting deltaH_hum vs deltaH_rnd
figfold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/figures/fig4/';
cmap = othercolor('Greys4',256);
net_labs = {'Human','Rand','DPR','DSPR'};
nnets = length(net_labs);
p_lims = [0 0.1];
ticks = linspace(0,0.1,11);

%% Plotting connectomes
figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','visible','on');
for n=1:nnets
    subplot(2,4,n)
    imagesc(C_cond(:,:,n),[0 0.2])
    colormap(cmap)
    axis square
    title(net_labs{n})
    set(gca,'xtick','','ytick','')
    
end
for n=1:3
    subplot(2,4,n+4+1)
    yy = squeeze(deltah(:,:,n+1));
%     plot(sel_delta_ent(:),yy(:),'k.');hold on
    p1=plot(mean(sel_delta_ent,2),mean(yy,2),'o','color','k','markersize',10,...
    'markerfacecolor',[0.7 0 0]);hold on
    plot(p_lims,p_lims,'k--');hold on
    set(gca,'xtick',ticks,'ytick',ticks)
    axis square;
    ylim(p_lims)
    xlim(p_lims)
    grid on;
    box off;
    legend(p1,['R = ',num2str(mean(r_sq_2_plot(:,n)).*1000./1000),' \pm ',...
        num2str(std(r_sq_2_plot(:,n)).*1000./1000)])
%         legend(p1,['R = ',num2str(mean(r_sq_2_plot(:,n)).*1000./1000),' \pm ',...
%         num2str(std(r_sq_2_plot(:,n)).*1000./1000)])
    xlabel(['\Delta H ',net_labs{1}])
    ylabel(['\Delta H ',net_labs{n+1}])
end

print(gcf,'-dpdf',[figfold,'delta_h_vs_connectomes.pdf'],'-r300')
%% Plotting PCB vs 2A

figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','visible','on');

for n=1:nnets
    subplot(1,4,n)
    plot(mean_pcb(:,n),mean_lsd(:,n),'k.','markersize',msize);hold on
    xlim([2 2.6])
    refline(1,0);
    title(net_labs{n})
    grid on    
    xlabel('PCB Entropy (nats)')        
    ylabel('2A Entropy')    
    axis square
end

print(gcf,'-dpdf',[figfold,'scatter_h_vs_connectomes.pdf'],'-r300')

%% Plotting the effect on different randomizations with the same initial condition
nrands = 30;
sel_ini = 1;
thisparlist = cell(nmods,nrands);

for r=1:nrands
    for m=1:nmods
        thispars = params;
        thispars.wgaine = sel_wgaine(m);
        thispars.wgaini = sel_wgaini(m);
        [c_dspr,R] = null_model_und_sign(C,0,0.1);
        thisstren = sum(c_dspr)./2;
        
        thispars.C = c_dspr;
        thispars.J = sel_alpha.*thisstren.*thispars.G + 1;
        thispars.seed = iniconds(sel_ini);
        thisparlist{m,r} = thispars;
    end
end

%% Simulating
thisreg_ent = zeros(N,nrands,nmods);
% T = 510;
T = 110;
nsteps = T.*(1000); % number of DMF timepoints
parfor r=1:nrands
    for m=1:nmods
        
        [rates,bold] = dyn_fic_DMF(thisparlist{m,r}, nsteps,'both'); % runs simulation
        rates = rates(:,(thisparlist{m,r}.burnout*1000*2):end);
        
        for n=1:N
            gamma_pars = gamfit(rates(n,:));
            thisreg_ent(n,r,m) = gamma_ent_fun(gamma_pars);
        end
    end
end

%%
thisdeltah = squeeze((thisreg_ent(:,:,2)-thisreg_ent(:,:,1))./thisreg_ent(:,:,1));
sel_delta_ent = squeeze(deltah(:,:,1));
r_sq_dspr_r = zeros(nrands,1);
for r=1:nrands
    r_sq_dspr_r(r) = corr(sel_delta_ent(:,sel_ini),thisdeltah(:,r));
end
