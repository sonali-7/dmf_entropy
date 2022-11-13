%% New supplementary figure 4 SciRep: Receptor controls
PARAMETERS_DIR=dotenv.read().PROJECT_DIR, 'parameters'
load(fullfile(PARAMETERS_DIR, 'SC_and_5ht2a_receptors.mat'))
C = sc90./max(sc90(:))*0.2;
N= length(C);
stren = sum(C)./2;
basefold = '/media/ruben/ssd240/Matlab/cb-neuromod-master/';
load([PARAMETERS_DIR,'fc_fcd_bold_sig_pcb_lsd.mat'],'fcd','fc','tr','flp','fhi',...
    'wsize','overlap','condnames','sel_conds')
ave_fc = squeeze(mean(fc,3));
vec_sc = squareform(C);
% Loading previous simulation initial condition to avoid re calculation of
% pcb and 2A condition
basefold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/';
prevdata = load([basefold,'dmf_pcb_lsd_nm_i_v5.mat']);
%% Loading receptors
load('mean5T_all.mat','symm_mean5HT2A','symm_mean5HT1A','symm_mean5HT1B','symm_mean5HT4','symm_mean5HTT') % loads receptor density
symm_mean5HT2A = symm_mean5HT2A./max(symm_mean5HT2A);
rec_uni = ones(size(symm_mean5HT2A))*mean(symm_mean5HT2A);
rec_rand = symm_mean5HT2A(randperm(N));
all_receptors = cat(2,symm_mean5HT1A./max(symm_mean5HT1A),...
    symm_mean5HT1B./max(symm_mean5HT1B),symm_mean5HT4./max(symm_mean5HT4),...
    symm_mean5HTT./max(symm_mean5HTT),rec_uni,rec_rand);

nmaps = size(all_receptors,2);
%% Preparing parameters
nreps = prevdata.nreps;
parlist = cell(nmaps,nreps);

for m=1:nmaps
    thispars = prevdata.params;
    thispars.receptors = all_receptors(:,m);
    for r=1:nreps
        thispars.wgaine = prevdata.sel_wgaine(2);
        thispars.wgaini = prevdata.sel_wgaini(2);
        thispars.seed = prevdata.iniconds(r);
        parlist{m,r} = thispars;
    end
end

%% Running simulation
nreps = 30;
nconds=2;
gamma_ent_fun = @(a) a(1) + log(a(2)) + log(gamma(a(1))) + (1-a(1))*psi(a(1));

isubfc = find(tril(ones(N),-1));

sel_ks_fcd = zeros(nreps,nmaps,nconds);
sel_fc_mse = sel_ks_fcd;
reg_fr = zeros(N,nreps,nmaps);
reg_ent = reg_fr;
bold_sigs = cell(nreps,nmaps);

T = 510;
% T = 110;
nsteps = T.*(1000); % number of DMF timepoints
%
parfor r=1:nreps
    tic
    for m=1:nmaps
        selpars = parlist{m,r};
        % Simulating

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
            sel_fc_mse(r,m) = mean((sim_fc(isubfc)-this_fc(isubfc)).^2); % MSE FC

            % Computing KS FCD for each condition and MSE on FC
            [~,~,aux_ks_fcd(c)] = kstest2(sim_fcd(:),this_fcd(:));
        end
        bold_sigs{r,m} = filt_bold;
        sel_ks_fcd(r,m,:) = aux_ks_fcd;
    end
    toc
end

savefold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/';
save([savefold,'dmf_receptor_controls.mat'],'iniconds',...
    'sel_ks_fcd','sel_fc_mse','reg_fr','reg_ent','bold_sigs','nreps');
%% Checking entropies
figfold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/figures/supp4/';
map_name = {'1A','2B','4','TT','uni','rand'};

pcb_ent = prevdata.reg_ent(:,1:nreps,1);
lsd_ent = prevdata.reg_ent(:,1:nreps,2);
delta_h_2a = (lsd_ent - pcb_ent)./pcb_ent;
mean_dh_2a = mean(delta_h_2a,2);

mean_ent = squeeze(mean(reg_ent,2));
deltah = (reg_ent -pcb_ent)./pcb_ent;
mean_dh = squeeze(mean(deltah,2));


%% Plotting scatter PCB vs LSD
msize = 20;
figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','visible','on');

for m=1:nmaps
    subplot(2,3,m)
    plot(mean(pcb_ent,2),mean(lsd_ent,2),'k.','markersize',msize);hold on
    plot(mean(pcb_ent,2),mean_ent(:,m),'r.','markersize',msize);hold on
    xlim([2 2.6])
    refline(1,0);
    title(map_name{m})
    grid on
    if m>3
        xlabel('PCB Entropy (nats)')
    end
    if m==1 || m==4
        ylabel('NM entropy')
    end
    legend('2A',map_name{m},'location','northwest')

end
%
print(gcf,'-dpdf',[figfold,'receptors_maps_vs_entropy.pdf'],'-r300')

%%
figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','visible','on');

for m=1:nmaps
    subplot(2,3,m)
    plot(stren,mean_dh_2a,'k.','markersize',msize);hold on
    plot(stren,mean_dh(:,m),'r.','markersize',msize);hold on
    xlim([0.1 0.5])

    title(map_name{m})
    grid on
    if m>3
        xlabel('Strength')
    end
    if m==1 || m==4
        ylabel('\Delta H')
    end
    legend('2A',map_name{m},'location','northwest')

end
print(gcf,'-dpdf',[figfold,'receptors_maps_vs_stren_vs_dh.pdf'],'-r300')
%% Plotting KS-FCD for receptors
xx= 1:3:(nmaps*3);
figure('units','normalized','outerposition',[0 0 1 1],'PaperOrientation','landscape','visible','on');
s1=plotSpread(prevdata.sel_ks_fcd(1:nreps,1,2),'xValues',-2,'showMM',5,'distributionColors','b');hold on
set(s1{2},'color','b')
s2=plotSpread(prevdata.sel_ks_fcd(1:nreps,2,2),'xValues',-1,'showMM',5,'distributionColors','r');hold on

s1=plotSpread(sel_ks_fcd(:,:,1),'xValues',xx,'showMM',5,'distributionColors','b');hold on
set(s1{2},'color','b')
s2=plotSpread(sel_ks_fcd(:,:,2),'xValues',xx+1,'showMM',5,'distributionColors','r');
set(gca,'xtick',[-1.5 xx+0.5],'xticklabel',['2A',map_name])
xlim([-3 nmaps*3])
legend([s1{1}(1) s2{1}(1)],'pcb','lsd')
ylabel('KS FCD')
xlabel('5-HT Receptor Map')
set(gca,'ygrid','on')

print(gcf,'-dpdf',[figfold,'receptors_maps_vs_ksfcd.pdf'],'-r300')
