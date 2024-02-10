%% This script creates TFRhann files
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
