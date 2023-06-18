run setsubj;

datapath = sprintf('%s/data.orig', rootpath);
% datapath = sprintf('%s/data.sfeq', rootpath);

gavgpath = sprintf('%s/grandavg',datapath);

%% Perform ERP analysis

postfix = '_ica_pruned_ar';

for isub = 1 : length(subj)
    
    % Subject path name
    subpath = fullfile(datapath, subj{isub});
    
    % Load dataset
    filename = [subj{isub} postfix '.set'];
    EEG = pop_loadset('filename', filename, 'filepath', subpath);
    
    % Averaging. Only good trials.  Include standard deviation.  Save to disk.
    ERP = pop_averager( EEG, 'Criterion', 'good', 'SEM', 'on');
    ERP.erpname = [subj{isub} '_ERPs'];  % name for erpset menu
    pop_savemyerp(ERP, 'erpname', ERP.erpname, 'filename', [ERP.erpname '.erp'], 'filepath', subpath, 'warning', 'off');

    % Bin Operations. 
    fname = [subj{isub} '_ERPs.erp'];  % Re-create filename for unfiltered ERP
    ERP = pop_loaderp( 'filename', fname, 'filepath', subpath);   % Load the file
    
    ERP = pop_binoperator( ERP, bineq_fn );
    ERP.erpname = [ERP.erpname '_diff'];  % name for erpset menu
    pop_savemyerp(ERP, 'erpname', ERP.erpname, 'filename', [ERP.erpname '.erp'], 'filepath', subpath, 'warning', 'off');

    % Filtering ERP. Channels = 1 to 64; No high-pass;
    % Lowpass cutoff at 30 Hz; Order of the filter = 2.
    % Type of filter = "Butterworth"; Do not remove DC offset
    fprintf('\n\n\n**** %s: Low-pass filtering ERP at 30 Hz ****\n\n\n', subj{isub});              
    ERP = pop_filterp( ERP, 1:ERP.nchan, 'Cutoff', 30, 'Design', 'butter', 'Filter', 'lowpass', 'Order',2 );
    ERP.erpname = [ERP.erpname '_30Hz'];  % name for erpset menu
    pop_savemyerp(ERP, 'erpname', ERP.erpname, 'filename', [ERP.erpname '.erp'], 'filepath', subpath, 'warning', 'off');
    
end

% Use checkindividuals.m to check individual ERP plots


%% Grand Average

% Subset of subjects with rejection rates less than 20%
subj_tmp = subj; % subj(rejrates < 25);

% Initialize the ALLERP structure and CURRENTERP
ALLERP = buildERPstruct([]);
CURRENTERP = 0;

% Kind of ERP data
postfix = '_ERPs_diff_30Hz';

for isub = 1 : length(subj_tmp)

    subpath = fullfile(datapath, subj_tmp{isub});
    fname = [subj_tmp{isub} postfix '.erp'];
    ERP = pop_loaderp( 'filename', fname, 'filepath', subpath);
    fprintf('Loaded %s\n', fname);
    fprintf(' Channels = %g\n', ERP.nchan);
    
    CURRENTERP = CURRENTERP + 1;
    ALLERP(CURRENTERP) = ERP;
end

if ~exist(gavgpath,'dir')
    mkdir(gavgpath);
end

% Make a grand average. The final ERP from each subject was saved in ALLERP
ERP = pop_gaverager( ALLERP, 'Criterion', 100, 'ERPindex', 1:length(ALLERP) ); %'Criterion', 100 is left over from a previous version
% name for erpset menu
ERP.erpname = ['grandavg' postfix];
ERP = pop_savemyerp(ERP, 'filename', [ERP.erpname '.erp'], 'filepath', gavgpath, 'warning', 'off');
CURRENTERP = CURRENTERP + 1;
ALLERP(CURRENTERP) = ERP;



%% Waveforms

run setsubj;

datapath = sprintf('%s/data.orig', rootpath);
% datapath = sprintf('%s/data.sfeq', rootpath);

cmap = gray(30);
set(0, 'DefaultAxesColorOrder', flipud(cmap([5 15 25],:)));

% load grandavg LME results
gavgpath = fullfile(datapath,'grandavg');
gERP = pop_loaderp('filename','grandavg_ERPs_diff_30Hz.erp', 'filepath',gavgpath);

h = figure('Position', [100 100 500 800]);

subplot(2,1,1);
% chan = [1 19];
chan = 17;
typenum = 1:3;

% Figure parameters
ylim = [-5 8];
xlim = [-100 500];
ytick = [-4 -2 0 2 4 6];
xtick = [-100:100:500];
xticklabel = {'-100 ms','','','','','',''};
yticklabel = {'','','','+2\muV','',''};

PlotAxisAtOrigin(gERP.times, squeeze(mean(gERP.bindata(chan,:,typenum),1)), ...
    'XLim', xlim, 'YLim', ylim, 'XTick', xtick, 'YTick', ytick, ...
    'XTickLabel', xticklabel, 'YTickLabel', yticklabel,'LineWidth',1);
legend('N low','N med','N high','Location','SouthEast');
title(sprintf('Ch %d',chan));

subplot(2,1,2);
% chan = 50;
chan = 20;
typenum = 1:3;

% Figure parameters
ylim = [-3 10];
xlim = [-100 500];
ytick = [-2 0 2 4 6 8];
xtick = [-100:100:500];
xticklabel = {'-100 ms','','','','','',''};
yticklabel = {'','','+2\muV','','',''};

PlotAxisAtOrigin(gERP.times, squeeze(mean(gERP.bindata(chan,:,typenum),1)), ...
    'XLim', xlim, 'YLim', ylim, 'XTick', xtick, 'YTick', ytick, ...
    'XTickLabel', xticklabel, 'YTickLabel', yticklabel,'LineWidth',1);
legend('N low','N med','N high','Location','SouthEast');
title(sprintf('Ch %d',chan));

[pth,dname,ext] = fileparts(datapath);
jp_print2pdf(h,fullfile(rootpath,'figs',sprintf('Waveforms_%s.pdf',[dname,ext])),'orientation','portrait','width',5);

