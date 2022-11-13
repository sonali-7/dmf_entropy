%% Workflow for computing FC and FCD for LSD and Rest data Deco et al 2018.
% Loading data
% 1,2,3 LSD rest1, LSD music, LSD rest2
% 4,5,6 PCB rest1, PCB music, PCB rest2
basefold = '/media/ruben/ssd240/Matlab/cb-neuromod-master/';
load([basefold,'LSDnew.mat'],'tc_aal')
% load([basefold,'SC_and_5ht2a_receptors.mat'])

%% Selecting only data with 200 or more datapoints
sel_conds = [5,2]; % PCB, LSD with music
bold_sigs = tc_aal(:,sel_conds);
subt = cellfun(@(x) size(x,2),bold_sigs);
nsubs = length(subt);
nconds = 2;
condnames = {'pcb','lsd'};
% sewl_by_ 
%% Setting filter for DMf-Simulated data
tr=2;
fnq=1/(2*tr);                 % Nyquist frequency
flp = 0.01;                   % lowpass frequency of filter
fhi = 0.1;                    % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter
nzeros = 40;

%% Computing FC and FCD for each subject
N=90;
isubfc = find(tril(ones(N),-1));
wsize = 30;
overlap = 28;
win_start = 0:wsize-overlap:subt(1)-wsize-1;
nwins = length(win_start);
isubfcd = find(tril(ones(nwins),-1));
nints = length(isubfcd);

fc = zeros(N,N,nsubs,nconds);
fcd = zeros(nwins,nwins,nsubs,nconds);
var_fcd = zeros(nsubs,nconds);
for c=1:nconds
    
    for s=1:nsubs
        
        fc(:,:,s,c) = corrcoef(bold_sigs{s,c}');
        aux_fcd = compute_fcd(bold_sigs{s,c}',wsize,overlap,isubfc);
        aux_fcd = corrcoef(aux_fcd);
        fcd(:,:,s,c) = aux_fcd;
        var_fcd(s,c) = var(aux_fcd(isubfcd));
        
    end
    
%     ave_fc = mean(fc,3) - eye(N);
%     vec_ave_fc = squareform(ave_fc);
%     vec_emp_fcd = squeeze(fcd(:,:,:));
end

% Saving FC and data parameters
save([basefold,'fc_fcd_bold_sig_pcb_lsd.mat'],'fcd','fc','tr','flp','fhi',...
    'wsize','overlap','condnames','sel_conds')