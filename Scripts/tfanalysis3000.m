clear;
datapath = 'C:\Users\Owner\Documents\CoDeNeuroLab\Data\data.orig';
subj = ["3000" "3001" "3002" "3003" "3004" "3005" "3006" "3007", ... 
    "3008" "3009" "3010" "3011" "3023" "3031"];

nsub = length(subj);
chan = [17];
timesout = 236;

FreqPow = struct([]);
number = [10 11 13 19 20 22 25 26 27; 
        4 5 7 12 14 16 21 23 24; 1 2 3 6 8 9 15 17 18];
sizes = [1 2 3 4 5 7 10 11 13; 
        6 8 9 12 14 16 19 20 22; 15 17 18 21 23 24 25 26 27];
spacing = [1 4 6 10 12 15 19 21 25; 
        2 5 8 11 14 17 20 23 26; 3 7 9 13 16 18 22 24 27];
    
icond = cat(3, number, sizes, spacing);
[numRows, numCols, numDepth] = size(icond);



alphaSumCond = zeros(3, 3, 236);
betaSumCond = zeros(3, 3, 236);
thetaSumCond = zeros(3, 3, 236);
gammaSumCond = zeros(3, 3, 236);

%isub = 1;
for isub = 1:nsub
    % Load data
    subpath = fullfile(datapath);
    EEG = pop_loadset('filename', sprintf('%s_ica_pruned_ar.set',subj(isub)), 'filepath', subpath);
    
    FreqPow(isub).cond = struct([]);
    
 %%
    for icond2 = 1:numDepth
        
        for icond3 = 1:numRows
            % Event code for each epoch (in the "eventbini" field)

            if iscell(EEG.epoch(1).eventbini)
                events_epochs = cellfun(@(x) x{1}, {EEG.epoch(:).eventbini});
            else
                events_epochs = cellfun(@(x) x(1), num2cell([EEG.epoch(:).eventbini]));
            end
            % Index for this condition and rejected trials
            idx = ismember(events_epochs,icond(icond2,icond3));
            if ~isempty(EEG.reject.rejmanual)
                idx_r = EEG.reject.rejmanual > 0;
            else
                idx_r = false(size(idx));
            end
    
           
            FreqPow(isub).cond(icond2, icond3).theta  = nan(length(chan),timesout);
            FreqPow(isub).cond(icond2, icond3).alpha  = nan(length(chan),timesout);
            FreqPow(isub).cond(icond2, icond3).beta  = nan(length(chan),timesout);
            FreqPow(isub).cond(icond2, icond3).gamma  = nan(length(chan),timesout);
            %[FreqPow(isub).cond(icond(1,:)).hgamma] = deal(nan(length(chan),timesout));
            FreqPow(isub).cond(icond2, icond3).ersp   = nan(51,236,length(chan));
          
            for ichan = 1 %: length(chan)
                % Current channel and condition
                thisdata = EEG.data(chan(1),:,~idx_r & idx);
    
                h = figure;
       
                [ersp,itc,powbase,times,freqs,erspboot,itcboot,tfdata] = newtimef(thisdata,EEG.pnts,[EEG.xmin EEG.xmax]*1000,...
                    EEG.srate,0,'timesout',timesout,'padratio',8);
                close(h);
                
             
                % Theta (4-8); Alpha (9-12); Beta (13-30); Gamma (30-80); High-gamma (80-200)           
                FreqPow(isub).cond(icond2, icond3).theta(ichan,:)  = mean(ersp(freqs > 4 & freqs < 8,:));
                FreqPow(isub).cond(icond2, icond3).alpha(ichan,:)  = mean(ersp(freqs > 9 & freqs < 12,:));
                FreqPow(isub).cond(icond2, icond3).beta(ichan,:)   = mean(ersp(freqs > 13 & freqs < 30,:));
                FreqPow(isub).cond(icond2, icond3).gamma(ichan,:)  = mean(ersp(freqs > 30 & freqs < 80,:));
                %FreqPow(isub).cond(cond_idx).hgamma(ichan,:) = mean(ersp(freqs > 80 & freqs < 200,:));
                FreqPow(isub).cond(icond2, icond3).ersp(:,:,ichan) = ersp;                        
                
            end
           
        end
    
    end

    SubFreqPow = FreqPow(isub);
    
    % Save preprocessed data
    save(fullfile(subpath,[num2str(subj(isub)) '_freqpow.mat']),'SubFreqPow','times','freqs');
    fprintf('### Subject %s saved. ###\n',subj(isub));

   
    %%

    % Alpha Power

    alphaTempCond = zeros(3, 3, 236);
    betaTempCond = zeros(3, 3, 236);
    thetaTempCond = zeros(3, 3, 236);
    gammaTempCond = zeros(3, 3, 236);
    for i = 1:3
        if isub == 1
            alphaSumCond(i, 1, :) = FreqPow(isub).cond(i, 1).alpha(ichan,:);
            alphaSumCond(i, 2, :) = FreqPow(isub).cond(i, 2).alpha(ichan,:);
            alphaSumCond(i, 3, :) = FreqPow(isub).cond(i, 3).alpha(ichan,:);
            betaSumCond(i, 1, :) = FreqPow(isub).cond(i, 1).beta(ichan,:);
            betaSumCond(i, 2, :) = FreqPow(isub).cond(i, 2).beta(ichan,:);
            betaSumCond(i, 3, :) = FreqPow(isub).cond(i, 3).beta(ichan,:);
            thetaSumCond(i, 1, :) = FreqPow(isub).cond(i, 1).theta(ichan,:);
            thetaSumCond(i, 2, :) = FreqPow(isub).cond(i, 2).theta(ichan,:);
            thetaSumCond(i, 3, :) = FreqPow(isub).cond(i, 3).theta(ichan,:);
            gammaSumCond(i, 1, :) = FreqPow(isub).cond(i, 1).gamma(ichan,:);
            gammaSumCond(i, 2, :) = FreqPow(isub).cond(i, 2).gamma(ichan,:);
            gammaSumCond(i, 3, :) = FreqPow(isub).cond(i, 3).gamma(ichan,:);
        else
            alphaTempCond(i, 1, :) = FreqPow(isub).cond(i, 1).alpha(ichan,:);
            alphaTempCond(i, 2, :) = FreqPow(isub).cond(i, 2).alpha(ichan,:);
            alphaTempCond(i, 3, :) = FreqPow(isub).cond(i, 3).alpha(ichan,:);
            betaTempCond(i, 1, :) = FreqPow(isub).cond(i, 1).beta(ichan,:);
            betaTempCond(i, 2, :) = FreqPow(isub).cond(i, 2).beta(ichan,:);
            betaTempCond(i, 3, :) = FreqPow(isub).cond(i, 3).beta(ichan,:);
            thetaTempCond(i, 1, :) = FreqPow(isub).cond(i, 1).theta(ichan,:);
            thetaTempCond(i, 2, :) = FreqPow(isub).cond(i, 2).theta(ichan,:);
            thetaTempCond(i, 3, :) = FreqPow(isub).cond(i, 3).theta(ichan,:);
            gammaTempCond(i, 1, :) = FreqPow(isub).cond(i, 1).gamma(ichan,:);
            gammaTempCond(i, 2, :) = FreqPow(isub).cond(i, 2).gamma(ichan,:);
            gammaTempCond(i, 3, :) = FreqPow(isub).cond(i, 3).gamma(ichan,:);
            alphaSumCond(i, 1, :) = alphaSumCond(i, 1, :) + alphaTempCond(i, 1, :);
            alphaSumCond(i, 2, :) = alphaSumCond(i, 2, :) + alphaTempCond(i, 2, :);
            alphaSumCond(i, 3, :) = alphaSumCond(i, 3, :) + alphaTempCond(i, 3, :);
            betaSumCond(i, 1, :) = betaSumCond(i, 1, :) + betaTempCond(i, 1, :);
            betaSumCond(i, 2, :) = betaSumCond(i, 2, :) + betaTempCond(i, 2, :);
            betaSumCond(i, 3, :) = betaSumCond(i, 3, :) + betaTempCond(i, 3, :);
            thetaSumCond(i, 1, :) = thetaSumCond(i, 1, :) + thetaTempCond(i, 1, :);
            thetaSumCond(i, 2, :) = thetaSumCond(i, 2, :) + thetaTempCond(i, 2, :);
            thetaSumCond(i, 3, :) = thetaSumCond(i, 3, :) + thetaTempCond(i, 3, :);
            gammaSumCond(i, 1, :) = gammaSumCond(i, 1, :) + gammaTempCond(i, 1, :);
            gammaSumCond(i, 2, :) = gammaSumCond(i, 2, :) + gammaTempCond(i, 2, :);
            gammaSumCond(i, 3, :) = gammaSumCond(i, 3, :) + gammaTempCond(i, 3, :);
        end

    end
    %%
end

% Alpha
titles = ["Number" "Size" "Spacing"];
figure;
t = tiledlayout(3, 1);
title(t, 'Average Signals Across All Subjects');
for i = 1:3
    nexttile;
    alphaAvgCond1 = squeeze(alphaSumCond(i, 1, :)) / nsub;
    alphaAvgCond2 = squeeze(alphaSumCond(i, 2, :)) / nsub;
    alphaAvgCond3 = squeeze(alphaSumCond(i, 3, :)) / nsub;
    dat = [alphaAvgCond1'; alphaAvgCond2'; alphaAvgCond3'];
    plot(times, dat);
    title(titles(i));
end    
legend({'low', 'medium', 'high'}, 'Location','northeast');
xlabel(t, 'Time (ms)');
ylabel(t, 'Alpha Power');

% Beta
figure;
t = tiledlayout(3, 1);
title(t, 'Average Signals Across All Subjects');
for i = 1:3
    nexttile;
    betaAvgCond1 = squeeze(betaSumCond(i, 1, :)) / nsub;
    betaAvgCond2 = squeeze(betaSumCond(i, 2, :)) / nsub;
    betaAvgCond3 = squeeze(betaSumCond(i, 3, :)) / nsub;
    dat = [betaAvgCond1'; betaAvgCond2'; betaAvgCond3'];
    plot(times, dat);
    title(titles(i));
end    
legend({'low', 'medium', 'high'}, 'Location','northeast');
xlabel(t, 'Time (ms)');
ylabel(t, 'Beta Power');

% Theta
figure;
t = tiledlayout(3, 1);
title(t, 'Average Signals Across All Subjects');
for i = 1:3
    nexttile;
    thetaAvgCond1 = squeeze(thetaSumCond(i, 1, :)) / nsub;
    thetaAvgCond2 = squeeze(thetaSumCond(i, 2, :)) / nsub;
    thetaAvgCond3 = squeeze(thetaSumCond(i, 3, :)) / nsub;
    dat = [thetaAvgCond1'; thetaAvgCond2'; thetaAvgCond3'];
    plot(times, dat);
    title(titles(i));
end    
legend({'low', 'medium', 'high'}, 'Location','northeast');
xlabel(t, 'Time (ms)');
ylabel(t, 'Theta Power');

% Gamma 
figure;
t = tiledlayout(3, 1);
title(t, 'Average Signals Across All Subjects');
for i = 1:3
    nexttile;
    gammaAvgCond1 = squeeze(gammaSumCond(i, 1, :)) / nsub;
    gammaAvgCond2 = squeeze(gammaSumCond(i, 2, :)) / nsub;
    gammaAvgCond3 = squeeze(gammaSumCond(i, 3, :)) / nsub;
    dat = [gammaAvgCond1'; gammaAvgCond2'; gammaAvgCond3'];
    plot(times, dat);
    title(titles(i));
end    
legend({'low', 'medium', 'high'}, 'Location','northeast');
xlabel(t, 'Time (ms)');
ylabel(t, 'Gamma Power');
    






% Calculate average ERSP across all subjects
ERSP_avg = zeros(numDepth, numRows, 51, timesout);

for icond2 = 1:numDepth
    for icond3 = 1:numRows
        ERSP_sum = zeros(51, timesout);
        
        for isub = 1:nsub
            ERSP_sum = ERSP_sum + FreqPow(isub).cond(icond2, icond3).ersp(:,:,1);
        end
        
        ERSP_avg(icond2, icond3, :, :) = ERSP_sum / nsub;
    end
end

% Plot heat maps
titles = ["Number" "Size" "Spacing"];
cond_labels = {'low', 'medium', 'high'};

for icond2 = 1:numDepth
    fig = figure; % Create a figure for each icond2 (Number, Size, Spacing)
    t = tiledlayout(fig, numRows, 1); % Create a tiled layout with numRows rows and 1 column
    title(t, titles(icond2));
    
    for icond3 = 1:numRows
        nexttile; % Move to the next tile in the tiled layout
        imagesc(times, freqs, squeeze(ERSP_avg(icond2, icond3, :, :)));
        axis xy;
        colorbar;
        clim([-2 2]);
        title(sprintf('%s', cond_labels{icond3}));
    end
    xlabel(t, 'Time (ms)');
    ylabel(t, 'Frequency (Hz)');
end
