% New figure of SciRep 2020.
% FIGURE 1C brains colored by entropy

basefold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/';
simdata = load([basefold,'dmf_pcb_lsd_nm_i_v5.mat']);
N = length(simdata.params.C);
nnms = 2;

%% Plotting both brains
sort_ids = cat(2,[1:2:N],fliplr([2:2:N]));
ave_reg_ent = squeeze(mean(simdata.reg_ent,2));
min_ent = min(ave_reg_ent(:));
max_ent = max(ave_reg_ent(:));
meshclims = [min_ent max_ent];
overlay_pla = zeros(N,1);
overlay_nm=overlay_pla;
overlay_pla(sort_ids) = ave_reg_ent(:,1);
overlay_nm(sort_ids) = ave_reg_ent(:,2);
overlay_all = zeros(N,2);
overlay_all(sort_ids,:) = ave_reg_ent;

cmap = othercolor('YlOrRd9',5000);
%% Setting the views and fig properties
v{1} = [270 0]; % L
v{2} = [0  90]; % Topo
v{3} = [90  0]; % R
zview={'L','T','R'};
nviews = length(v);
hemis = {'both','left','right'};
nhemis = length(hemis);
conds = {'pla','nm'};
method  = 'euclidean';  
figfold = '/media/ruben/ssd240/Matlab/fastdmf-master/newSciRep/figures/brains/';
meshclims = [2 2.6];
%% Extracting zenital plot of brain
for nm=1:nnms
    for i=1:nviews % views
        for h=1:nhemis
            figure();
            atemplate('overlay',overlay_all(:,nm),'method',method,'meshclims',meshclims,...
                'hemi',hemis{h},'nocolbar','meshcols',cmap); % g
            view(gca,v{i});
            cl =camlight('headlight');
            material shiny
            set(cl,'Color',[0.3 0.3 0.3])
            alpha 0.95
            print(gcf,'-dpng',[figfold,conds{nm},'_',hemis{h},'_',zview{i},'.png'],'-r300')
            %             print(gcf,'-dpdf',[figfold,conds{nm},'_',hemis{h},'_',zview{i},'.pdf'],'-bestfit','-r300')            
            close gcf
            
        end
    end
end
