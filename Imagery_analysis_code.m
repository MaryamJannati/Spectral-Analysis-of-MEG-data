addpath('C:\Users\Dan\Documents\Master thesis Trento\Codes\no fidelity codes\fieldtrip-20220407');
datadir = 'C:\Users\Dan\Documents\Master thesis Trento\Codes\scoutexport-fidelity';
savedir= 'C:\Users\Dan\Documents\Master thesis Trento\Codes\fidelity codes\TFRhanns';
Ffiles = dir(fullfile(datadir, '*Face*.mat'));
Pfiles = dir(fullfile(datadir, '*Place*.mat'));

%%%%
% 241-244 = none,poor,fair,vivid
% 231-234 = Face,Face,Place,Place
%%%%
cfg_freq              = [];
cfg_freq.output       = 'pow';
cfg_freq.channel      = 'all';
cfg_freq.method       = 'mtmconvol';
cfg_freq.taper        = 'hanning';
cfg_freq.foi          = 1:0.2:15;                         % analysis 1 to 15 Hz in steps of 0.2 Hz
cfg_freq.t_ftimwin    = ones(length(cfg_freq.foi),1).*0.3;   % length of time window = 0.3 sec
cfg_freq.toi          = -0.2:0.050: 6;                  % time window "slides" from -0.2 to 6 sec in steps of 0.05
cfg_freq.pad          = 10;%'maxperlen';
%%%
% For the HIGH FREQUENCY analysis, update the cfg_freq as below:

% cfg_freq.foi          = 15:1:100;                         % analysis 15 to 100 Hz in steps of 1 Hz
% cfg_freq.t_ftimwin    = ones(length(cfg_freq.foi),1).*0.3;   % length of time window = 0.3 sec

%% frequency alasysis FACE (ft_freqanalysis): 
   for iCond = 1:numel (Ffiles)  
       fa = fullfile(datadir, Ffiles(iCond).name);
       load(fa);               
       for i = 1:length(finalStruct.trial)
           if finalStruct.confirmation(i) > 232 || finalStruct.fidelity(i) < 243 || isnan(finalStruct.confirmation(i)) || isnan(finalStruct.fidelity(i))
              finalStruct.trial{i}(:,:) = NaN;
           end
       end        
       TFRhann = ft_freqanalysis(cfg_freq, finalStruct);        
       filename = [Ffiles(iCond).name(1:5) 'face_lowFreq'];  % For the HIGH FREQUENCY analysis, replace 'lowFreq' with 'highFreq' in this line    
       save (fullfile(savedir, [filename '_TFRhann']), 'TFRhann');
   end
%save (fullfile(savedir, ('MRCN_face_highFreq_TFRhann')), 'TFRhann'); % in
%case of single subject saving
%% frequency alasysis PLACE (ft_freqanalysis):   
   for iCond = 1:numel (Pfiles)  
       Pl = fullfile(datadir, Pfiles(iCond).name);
       load(Pl);
       for i = 1:length(finalStruct.trial)
           if finalStruct.confirmation(i) < 233 || finalStruct.fidelity(i) < 243 || isnan(finalStruct.confirmation(i)) || isnan(finalStruct.fidelity(i))
              finalStruct.trial{i}(:,:) = NaN;
           end
       end   
       TFRhann = ft_freqanalysis(cfg_freq, finalStruct);
       filename = [Pfiles(iCond).name(1:5) 'place_lowFreq'];  % For the HIGH FREQUENCY analysis, replace 'lowFreq' with 'highFreq' in this line 
       save (fullfile(savedir, [filename '_TFRhann']), 'TFRhann');
   end
%%%%
parcelsR = TFRhann.label(181:360); % defining  each hemisphere's parcels to plot
parcelsL = TFRhann.label (1:180);

%%%%%%%%%%%%
%% Plotting (ft_singleplotTFR):
%  Flowfiles = dir(fullfile(savedir, '*face_lowFreq_TFRhann*')); %% dirs of different conditions to plot
%  Plowfiles = dir(fullfile(savedir, '*place_lowFreq_TFRhann*'));
  Fhighfiles = dir(fullfile(savedir, '*face_highFreq_TFRhann*'));
%  Phighfiles = dir(fullfile(savedir, '*place_highFreq_TFRhann*'));

%% Right Hemisphere
figure;
screenSize = get(groot, 'screensize');
set(gcf, 'Position', screenSize);
for i = 1:numel(Fhighfiles) % replace 'Flowfiles' with 'PlowFiles', 'Fhighfiles', or 'Phighfiles'
     TFRhann = fullfile(savedir, Fhighfiles(i).name); % replace 'Flowfiles' with 'PlowFiles', 'Fhighfiles', or 'Phighfiles'
     load (TFRhann); 
     filename = [Fhighfiles(i).name(1:10) 'highFreq']; % replace  'Flowfiles'    with 'PlowFiles',     'Fhighfiles',    or 'Phighfiles'
     % for the HIGH FREQUENCY analysis, replace 'lowFreq' with 'highFreq'
     % for the place condition, replace (1:10) with (1:11)                                                       
        for j = 1:numel(parcelsR) 
            cfg = [];
            cfg.baseline     = [-0.2 0];
            cfg.baselinetype ='relchange';
            cfg.parameter    = 'powspctrm';
            cfg.masknans     = 'no';
            cfg.layout       = 'vertical';
            cfg.channel      = parcelsR{j};                                
            cfg.figure       = gcf;
            subplot(12,15,j)
            ft_singleplotTFR(cfg, TFRhann);                                     
        end
        saveas(gcf, fullfile(savedir, [filename '_ftplotTFR_RH.jpg']));
end

%% Left Hemisphere
% figure;
% screenSize = get(groot, 'screensize');
% set(gcf, 'Position', screenSize);
% for i = 1:numel(Flowfiles) % replace 'Flowfiles' with 'PlowFiles', 'Fhighfiles', or 'Phighfiles'
%      TFRhann = fullfile(savedir, Flowfiles(i).name); % replace 'Flowfiles' with 'PlowFiles', 'Fhighfiles', or 'Phighfiles'
%      load (TFRhann); 
%      filename = [Flowfiles(i).name(1:10) 'lowFreq']; % replace  'Flowfiles'    with 'PlowFiles',     'Fhighfiles',    or 'Phighfiles'
     % for the HIGH FREQUENCY analysis, replace 'lowFreq' with 'highFreq'
     % for the place condition, replace (1:10) with (1:11)   
%         for j = 1:numel(parcelsL)            
%            cfg = [];
%             cfg.baseline     = [-0.2 0];
%             cfg.baselinetype ='relchange';
%             cfg.parameter    = 'powspctrm';
%             cfg.masknans     = 'no';
%             cfg.layout       = 'vertical';
%             cfg.channel      = parcelsL{j};                                
%             cfg.figure       = gcf;
%             subplot(12,15,j)
%             ft_singleplotTFR(cfg, TFRhann)                                     
%         end
%         saveas(gcf, fullfile(savedir, [filename '_ftplotTFR_LH.jpg']));
%  end

%% Grand Average %%
%% Paths:
addpath('C:\Users\Dan\Documents\Master thesis Trento\Codes\no fidelity codes\fieldtrip-20220407');
subdirg= 'C:\Users\Dan\Documents\Master thesis Trento\Codes\fidelity codes\TFRhanns';
savedirg =  'C:\Users\Dan\Documents\Master thesis Trento\Codes\fidelity codes\grand average';

%% Data directories:
FLfiles = dir(fullfile(subdirg, '*face_lowFreq_TFRhann*'));
FHfiles = dir(fullfile(subdirg, '*face_highFreq_TFRhann*'));
PLfiles = dir(fullfile(subdirg, '*place_lowFreq_TFRhann*'));
PHfiles = dir(fullfile(subdirg, '*place_highFreq_TFRhann*'));
%%
% Concatenating the TFRhann files to the 4th dimension for the FACE/LOWFREQ data: 
facelowFreq_cat = [];
for i = 1:numel(FLfiles)    
        FL=fullfile(subdirg, FLfiles(i).name);
        load (FL);         
        facelowFreq_cat = cat(4, facelowFreq_cat, TFRhann.powspctrm);
end     
        save (fullfile(savedirg, 'facelowFreq_cat'), 'facelowFreq_cat');
% Computing the Grand Average for the condition : face low frequency        
facelowFreq_GA = squeeze(nanmean(facelowFreq_cat, 4));
facelowFreq_GA_TFR = TFRhann;
facelowFreq_GA_TFR.powspctrm = facelowFreq_GA;
save (fullfile(savedirg, 'facelowFreq_GA_TFR'), 'facelowFreq_GA_TFR');
%%        
% Concatenating the TFRhann files to the 4th dimension for the PLACE/LOWFREQ data:
placelowFreq_cat = [];
for i = 1:numel(PLfiles)    
        PL=fullfile(subdirg, PLfiles(i).name);
        load (PL);         
        placelowFreq_cat = cat(4, placelowFreq_cat, TFRhann.powspctrm);
end     
        save (fullfile(savedirg, 'placelowFreq_cat'), 'placelowFreq_cat');
% Computing the Grand Average for the condition : place low frequency         
placelowFreq_GA = squeeze(nanmean(placelowFreq_cat, 4));
placelowFreq_GA_TFR = TFRhann;
placelowFreq_GA_TFR.powspctrm = placelowFreq_GA;
save (fullfile(savedirg, 'placelowFreq_GA_TFR'), 'placelowFreq_GA_TFR');  
%%
% Concatenating the TFRhann files to the 4th dimension for the FACE/HIGHFREQ data: 
facehighFreq_cat = [];
for i = 1:numel(FHfiles)    
        FH=fullfile(subdirg, FHfiles(i).name);
        load (FH);         
        facehighFreq_cat = cat(4, facehighFreq_cat, TFRhann.powspctrm);
end     
        save (fullfile(savedirg, 'facehighFreq_cat'), 'facehighFreq_cat');
% Computing the Grand Average for the condition : face high frequency       
facehighFreq_GA = squeeze(nanmean(facehighFreq_cat, 4));
facehighFreq_GA_TFR = TFRhann;
facehighFreq_GA_TFR.powspctrm = facehighFreq_GA;
save (fullfile(savedirg, 'facehighFreq_GA_TFR'), 'facehighFreq_GA_TFR');
%%
% Concatenating the TFRhann files to the 4th dimension for the PLACE/HIGHFREQ data:
placehighFreq_cat = [];
for i = 1:numel(PHfiles)    
        PH=fullfile(subdirg, PHfiles(i).name);
        load (PH);         
        placehighFreq_cat = cat(4, placehighFreq_cat, TFRhann.powspctrm);
end     
        save (fullfile(savedirg, 'placehighFreq_cat'), 'placehighFreq_cat');
% Computing the Grand Average for the condition : place high frequency         
placehighFreq_GA = squeeze(nanmean(placehighFreq_cat, 4));
placehighFreq_GA_TFR = TFRhann;
placehighFreq_GA_TFR.powspctrm = placehighFreq_GA;
save (fullfile(savedirg, 'placehighFreq_GA_TFR'), 'placehighFreq_GA_TFR');
%%
parcelsR = TFRhann.label(181:360); % defining  each hemisphere's parcels to plot
parcelsL = TFRhann.label (1:180);
%% Plotting the spectogram of Grand Average for each condition (ft_singleplotTFR)
% Right hemisphere:
figure;
screenSize = get(groot, 'screensize');
set(gcf, 'Position', screenSize);
for i = 1:numel(parcelsR)                         
         cfg = [];                
         cfg.baseline     = [-0.2 0];
         cfg.baselinetype ='relchange';
         cfg.channel      = parcelsR{i};
         cfg.parameter    = 'powspctrm';
         cfg.masknans     = 'no';
         cfg.layout       = 'vertical';
         cfg.figure       = gcf;
         subplot (12,15,i)
         ft_singleplotTFR (cfg, facehighFreq_GA_TFR); % replace 'facelowFreq_GA_TFR' with 'placelowFreq_GA_TFR', or 'facehighFreq_GA_TFR', or 'facehighFreq_GA_TFR'                                                     
         caxis ([-0.4 0.4]);
end
         saveas(gcf, fullfile(savedirg,'facehighFreq_GA_RH.jpg')); % replace 'facelowFreq_GA_RH' with 'placelowFreq_GA_RH', or 'facehighFreq_GA_RH', or 'facehighFreq_GA_RH'
         
% % Left hemisphere:
% figure;
% screenSize = get(groot, 'screensize');
% set(gcf, 'Position', screenSize);
% for i = 1:numel(parcelsL)                         
%          cfg = [];                
%          cfg.baseline     = [-0.2 0];
%          cfg.baselinetype ='relchange';
%          cfg.channel      = parcelsL{i};
%          cfg.parameter    = 'powspctrm';
%          cfg.masknans     = 'no';
%          cfg.layout       = 'vertical';
%          cfg.figure       = gcf;
%          subplot (12,15,i)
%          ft_singleplotTFR(cfg, facelowFreq_GA_TFR); % replace 'facelowFreq_GA_TFR' with 'placelowFreq_GA_TFR', or 'facehighFreq_GA_TFR', or 'facehighFreq_GA_TFR'   
%          caxis ([-0.4 0.4]); 
% end
%          saveas(gcf, fullfile(savedirg,'facelowFreq_GA_LH.jpg'));  % replace 'facelowFreq_GA_LH' with 'placelowFreq_GA_LH', or 'facehighFreq_GA_LH', or 'facehighFreq_GA_LH'

           

%% Delta %%
%% Paths:
addpath('C:\Users\Dan\Documents\Master thesis Trento\Codes\no fidelity codes');
addpath('C:\Users\Dan\Documents\Master thesis Trento\Codes\no fidelity codes\fieldtrip-20220407');
subdird =  'C:\Users\Dan\Documents\Master thesis Trento\Codes\fidelity codes\grand average';
savedird =  'C:\Users\Dan\Documents\Master thesis Trento\Codes\fidelity codes\delta';
%%
% Computing the delta: LOWFREQ (FACE-PLACE)
lowFreq_delta_cat = facelowFreq_cat - placelowFreq_cat;
save (fullfile(savedird, 'lowFreq_delta_cat'), 'lowFreq_delta_cat');
lowFreq_delta_TFR = TFRhann; % load a low freq struct from TFRhanns
lowFreq_delta_TFR.powspctrm = lowFreq_delta_cat;
lowFreq_delta_TFR.dimord = 'chan_freq_time_subj';
save (fullfile(savedird, 'lowFreq_delta_TFR'), 'lowFreq_delta_TFR');
% Average of the low frequency delta:
lowFreq_delta_Avg = squeeze(nanmean(lowFreq_delta_cat,4));
lowFreq_delta_Avg_TFR = TFRhann; % load a low freq struct from TFRhanns
lowFreq_delta_Avg_TFR.powspctrm = lowFreq_delta_Avg;
save (fullfile(savedird, 'lowFreq_delta_Avg_TFR'), 'lowFreq_delta_Avg_TFR');
%%
% Computing the delta: HIGHFREQ (FACE-PLACE)
highFreq_delta_cat = facehighFreq_cat - placehighFreq_cat;
save (fullfile(savedird, 'highFreq_delta_cat'), 'highFreq_delta_cat');
highFreq_delta_TFR = TFRhann; % load a low freq struct from TFRhanns
highFreq_delta_TFR.powspctrm = highFreq_delta_cat;
highFreq_delta_TFR.dimord = 'chan_freq_time_subj';
save (fullfile(savedird, 'highFreq_delta_TFR'), 'highFreq_delta_TFR');
% Average of the high frequency delta:
highFreq_delta_Avg = squeeze(nanmean(highFreq_delta_cat,4));
highFreq_delta_Avg_TFR = TFRhann; % load a low freq struct from TFRhanns
highFreq_delta_Avg_TFR.powspctrm = highFreq_delta_Avg;
save (fullfile(savedird, 'highFreq_delta_Avg_TFR'), 'highFreq_delta_Avg_TFR');

%% Plotting the Avg Delta:

% Right hemisphere:
figure;
screenSize = get(groot, 'screensize');
set(gcf, 'Position', screenSize);
for i = 1:numel(parcelsR)                         
         cfg = [];                
         cfg.baseline     = 'no';
         %cfg.baselinetype ='relchange';
         cfg.channel      = parcelsR{i};
         cfg.parameter    = 'powspctrm';
         cfg.masknans     = 'no';
         cfg.layout       = 'vertical';
         cfg.figure       = gcf;
         subplot (12,15,i) 
         ft_singleplotTFR(cfg, highFreq_delta_Avg_TFR); % replace with 'lowFreq_delta_Avg_TFR'        
         caxis ([-1*10^-21 1*10^-21]);  % ([-1*10^-22 1*10^-22]) for low frequency      
         set(gcf, 'Colormap', rbw_db2); % load MyColormaps.mat
%   smoothed_TFR = smoothdata(lowFreq_delta_Avg_TFR.powspctrm,'gaussian', 2);  % uncomment for high frequency plotting
%   lowFreq_delta_Avg_TFR.powspctrm = smoothed_TFR;
end
         saveas(gcf, fullfile(savedird, 'highFreq_delta_Avg_RH21.jpg'));  % replace with 'highFreq_delta_Avg_RH' 
         
% % Left hemisphre:
% figure;
% screenSize = get(groot, 'screensize');
% set(gcf, 'Position', screenSize);
% for i = 1:numel(parcelsL)                
%          subplot (12,15,i) 
%          cfg = [];                
%          cfg.baseline     = 'no';
%          %cfg.baselinetype ='relchange';
%          cfg.channel      = parcelsL{i};
%          cfg.parameter    = 'powspctrm';
%          cfg.masknans     = 'no';
%          cfg.layout       = 'vertical';
%          cfg.figure       = gcf;
%          ft_singleplotTFR(cfg, lowFreq_delta_Avg_TFR); % replace with 'highFreq_delta_Avg_TFR' 
%          caxis ([-1*10^-21 1*10^-21]);
%          set(gcf, 'Colormap', rbw_db2);
%   smoothed_TFR = smoothdata(lowFreq_delta_Avg_TFR.powspctrm,'gaussian', 2);  % uncomment for high frequency plotting
%   lowFreq_delta_Avg_TFR.powspctrm = smoothed_TFR;
% end
%          saveas(gcf, fullfile(savedird, 'lowFreq_delta_Avg_LH.jpg')); % replace with 'highFreq_delta_Avg_LH' 
%%

%%%%%%%%%%%%%%%%%%%%
%% Statistics & brain plots %%
 savedirs =  'C:\Users\Dan\Documents\Master thesis Trento\Codes\fidelity codes\statistics';

load ('C:\Users\Dan\Documents\Master thesis Trento\Codes\fidelity codes\delta\lowFreq_delta_TFR.mat');
load ('C:\Users\Dan\Documents\Master thesis Trento\Codes\fidelity codes\delta\highFreq_delta_TFR.mat');
cd 'C:\Users\Dan\Documents\Master thesis Trento\Codes\no fidelity codes';
addpath ('C:\Users\Dan\Documents\Master thesis Trento\Codes\fidelity codes');
% Defining ROIs
ROIs ={'R_FFC_ROI R','R_H_ROI R','R_IFJa_ROI R','R_IFJp_ROI R','R_LO1_ROI R','R_LO2_ROI R','R_LO3_ROI R',...
      'R_PCV_ROI R','R_PHA1_ROI R','R_PHA2_ROI R','R_PHA3_ROI R','R_PIT_ROI R','R_POS1_ROI R','R_POS2_ROI R',...
      'R_RSC_ROI R','R_STSda_ROI R','R_STSdp_ROI R','R_STSva_ROI R','R_STSvp_ROI R'};
% region number: 243,250,251,252,263,264,265,276,286.287,288,291,293,294,303,307,308,309,310

    parcels = lowFreq_delta_TFR.label (1:360); %lowFreq_delta_TFR
% Find indices of elements in parcelsR that match elements in ROIs
[~, indices] = ismember (ROIs, parcels);

% Defining frequency band limits
bands = {    
    'Delta', [1, 3];
    'Theta' , [4, 7];
    'Alpha', [7, 13];   
    'Beta', [18, 25];
    'Gamma', [30, 100 ];
};
alpha = 0.05; % The critical value for statistical tests

for cond = 1:size (bands, 1)
    bandName = bands{cond, 1};
    bandLimits = bands{cond, 2};

    if cond == 1 || cond == 2 || cond == 3
         data = analyze_band (lowFreq_delta_TFR, bandLimits, ROIs, parcels); 
    else
         data = analyze_band (highFreq_delta_TFR, bandLimits, ROIs, parcels);  
    end
 rh_data = data (181:360, :);
 stat_results = table('Size', [length(rh_data), 7], ...
                         'VariableTypes', {'logical', 'double', 'double', 'double', 'logical', 'double', 'logical'}, ...
                         'VariableNames', {'H_shapiro', 'p_shapiro', 'w_shapiro', 'p', 'H', 'adj_p', 'h_fdr'});
 variableNames = {'H_shapiro', 'p_shapiro', 'w_shapiro', 'p', 'H', 'adj_p', 'h_fdr'};
for col = variableNames
    stat_results.(col{:}) = NaN(size(stat_results, 1), 1);
end
    for i = 1:size (rh_data,1) 
        if rh_data (i,1) ~= 0
           [H_shapiro, p_shapiro, W_shapiro] = swtest (rh_data (i, :));
          stat_results.H_shapiro(i) = H_shapiro;        
          stat_results.p_shapiro(i)  = p_shapiro;
          stat_results.w_shapiro(i)  = W_shapiro;       
        end
    end
 normality_ratio = sum (stat_results.H_shapiro == 0) / sum (stat_results.H_shapiro == 1);
 disp(['Normality Ratio:' num2str(normality_ratio)]);
if normality_ratio > 1
    for i = 1:size (rh_data, 1)
        if rh_data (i,1) ~= 0
        [H, p] = ttest (rh_data (i,:), 0,'alpha', alpha);
                 stat_results.p (i) = p;
                 stat_results.H (i) = H;             
        end
    end
else
    for i = 1:size (rh_data, 1)
        if rh_data (i,1) ~= 0
        [p, h] = signrank (rh_data (i,:), 0,'method', 'approximate','alpha', alpha); 
              stat_results.p (i) = p;
              stat_results.H (i) = h;
        end
    end
end
  p_values = stat_results.p (~isnan (stat_results.p));
  [h_fdr, crit_p, adj_ci_cvrg, adj_p] = fdr_bh (p_values, 0.05, 'pdep', 'yes');
  valid_rows = find (~isnan(stat_results.p));
  stat_results.adj_p (valid_rows) = adj_p;
  stat_results.h_fdr (valid_rows) = h_fdr;   
   for i = 1:size (rh_data, 1)
       if stat_results.h_fdr (i) == 0
           rh_data (i, :) = 0;               
       end
   end        
data (181:360, :) = rh_data;
vector = mean (data,2);

for i = 1:size (vector)
    if vector (i,1) ~= 0 
        disp (['Significant ROI after FDR correction: ', parcels{i}]);
    end
end

if sum (vector(:,1)) ~= 0
       brainplotter_db2(vector); 
       saveas (gcf, fullfile (savedirs, [bandName, '_brainplot.jpg']));
end
save(fullfile (savedirs, [bandName, '_stat_results']), 'stat_results');
end

%% Plotting the signifacnt areas
savedir =  'C:\Users\Dan\Desktop';
addpath('C:\Users\Dan\Documents\Master thesis Trento\Codes\no fidelity codes\fieldtrip-20220407');
load ('C:\Users\Dan\Documents\Master thesis Trento\Codes\fidelity codes\delta\lowFreq_delta_Avg_TFR.mat');
load ('C:\Users\Dan\Documents\Master thesis Trento\Codes\fidelity codes\grand average\facelowFreq_GA_TFR.mat')
load ('C:\Users\Dan\Documents\Master thesis Trento\Codes\fidelity codes\grand average\placelowFreq_GA_TFR.mat')
load 'C:\Users\Dan\Documents\Master thesis Trento\Codes\no fidelity codes\MyColormaps.mat';
  parcelsR = {'R_FFC_ROI R','R_H_ROI R','R_IFJa_ROI R','R_IFJp_ROI R','R_LO1_ROI R','R_LO2_ROI R','R_LO3_ROI R',...
      'R_PCV_ROI R','R_PHA1_ROI R','R_PHA2_ROI R','R_PHA3_ROI R','R_PIT_ROI R','R_POS1_ROI R','R_POS2_ROI R',...
      'R_RSC_ROI R','R_STSda_ROI R','R_STSdp_ROI R','R_STSva_ROI R','R_STSvp_ROI R'};
%parcelsR = { 'R_FFC_ROI R', 'R_H_ROI R'};

% spectral GA
figure;
set(gcf, 'Position');
for j = 1:numel(parcelsR) 
            cfg = [];
            cfg.baseline     = [-0.2 0];
            cfg.baseline     = 'yes';
            cfg.baselinetype ='relchange';
            cfg.parameter    = 'powspctrm';
            cfg.masknans     = 'no';
            cfg.layout       = 'vertical';
            cfg.channel      = parcelsR{j};            
            cfg.figure       = gcf;
            subplot(7,3,j)
            caxis ([-0.5 0.5]);
            ft_singleplotTFR(cfg, placehighFreq_GA_TFR);                                     
end
        saveas(gcf, fullfile(savedir, 'place_high_GA.jpg'));   %'place_sig_GA.jpg'

% delta

figure;
set(gcf, 'Position');
for i = 1:numel(parcelsR)                         
         cfg = [];                
         cfg.baseline     = 'no';
         %cfg.baselinetype ='relchange';
         cfg.channel      = parcelsR{i};
         cfg.parameter    = 'powspctrm';
         cfg.masknans     = 'no';
         cfg.layout       = 'vertical';
         cfg.figure       = gcf;
         subplot (2,1,i) 
         ft_singleplotTFR(cfg, lowFreq_delta_Avg_TFR); % replace with 'lowFreq_delta_Avg_TFR'        
         caxis ([-1*10^-21.3 1*10^-21.3]);  % ([-1*10^-22 1*10^-22]) for low frequency      
         set(gcf, 'Colormap', rbw_db2); % load MyColormaps.mat
end
         saveas(gcf, fullfile(savedir, 'delta_sig.jpg'));  
         %%%%%%%%%%%%%%%%%%%%%%%%
