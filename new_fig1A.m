%% New SciRep Paper fig # 1A
load(fullfile(dotenv.read().PROJECT_DIR, 'parameters', 'SC_and_5ht2a_receptors.mat'))
C = sc90./max(sc90(:))*0.2;
stren = sum(C)./2;
N = length(C);
selG = 2.4;
sel_alpha = 1.5;
sel_wgaine = [0 0.024];
sel_wgaini = [0 0.024];

% Default pars
[ params ] = dyn_fic_DefaultParams('C',C);
params.burnout = 10;
params.batch_size = 50000;
params.receptors = receptors;
params.lrj = 0;
params.taoj = Inf;
params.G = selG;
params.J = sel_alpha*params.G*stren' + 1; % updates
params.seed = 1;
gamma_ent_fun = @(a) a(1) + log(a(2)) + log(gamma(a(1))) + (1-a(1))*psi(a(1));
%% Running model
% T = 510;
T = 110;
nsteps = T.*(1000); % number of DMF timepoints
tic
reg_ent = zeros(N,1);
gamma_pars = zeros(2,N);
[rates,bold] = dyn_fic_DMF(params, nsteps,'both'); % runs simulation
rates = rates(:,(params.burnout*1000*2):end);
for n=1:N
    gamma_pars(:,n) = gamfit(rates(n,:));
    reg_ent(n,1) = gamma_ent_fun(gamma_pars(:,n));
end
toc
%% Plotting the brain

sel_nodes_id = [24,71,68,45];

%% 2 & 3.- Plot firing rates of highlighted regions and Plot firing rates distribution, gamma fit and entropy
xx=linspace(0,20,10000);
tini = 10000;
tend = 20000;
xx_fr = (tini:tend)*params.dt/1000;
xx_fr = xx_fr - min(xx_fr);
figure
for i=1:4

    subplot(4,2,2*i-1)
    plot(xx_fr,rates(sel_nodes_id(i),tini:tend),'color',[0 0 0.75])
    set(gca,'xtick',0:0.25:1)
    if i<4
        set(gca,'xticklabel','')
    else
        xlabel('Time (s)')
    end
    ylabel('Firing Rate (Hz)')
    ylim([0 15])
    xlim([0 max(xx_fr)])
    box off


    gpdf = gampdf(xx,gamma_pars(1,sel_nodes_id(i)),gamma_pars(2,sel_nodes_id(i)));
    subplot(4,2,2*i)
    p1=plot(xx,gpdf,'linewidth',3,'color',[0 0 0.4]);hold on
    histogram(rates(sel_nodes_id(i),:),'Normalization','pdf',...
        'FaceColor',[0 0 0.7],'EdgeColor','none','FaceAlpha',0.6);hold on
    xlim([0 15])
    ylim([0 0.26])
    box off
    if i<4
        set(gca,'xticklabel','')
    else
        xlabel('Firing Rate (Hz)')
    end
    ylabel ('PDF')
    legend(['H(X) = ',num2str(reg_ent(sel_nodes_id(i)))])
end

%% 4.- Topographical Distribution of Entropy.
load('dmf_simulation_fig2_v3.mat')
overlay = zeros(N,1);
ave_reg_ent = mean(reg_ent,3);
overlay(sort_ids) = ave_reg_ent;
cmap = othercolor('YlOrRd9',5000);
figure;
aa=atemplate('network',sc_thr>0,'overlay',overlay,'method',method,...
    'hemi','both','nocolbar','meshcols',cmap); % g
%
% set(aa.overlay.cb,'ColorMap','custom')
% set(aa.overlay.cb,'CLim',[min(reg_ent) max(reg_ent)])