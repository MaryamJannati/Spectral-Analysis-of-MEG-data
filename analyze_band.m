function data = analyze_band(tfr_data, bandLimits, ROIs, parcels)
% Regions of interest that are defined in the main code
% parcels = TFRstructure.label (1:360); 
    [~, indices] = ismember(ROIs, parcels);
    bandlimits = find(tfr_data.freq > bandLimits(1) & tfr_data.freq < bandLimits(2));
    time_window = find(tfr_data.time >= 0 & tfr_data.time <= 6);   
    data = squeeze(nanmean(nanmean(tfr_data.powspctrm(:, bandlimits, time_window, :), 2), 3));   
    data(setdiff(1:numel(parcels), indices)) = 0;
    for i = 1:size (data,1)
        if data(i,1)==0
           data(i,:)=0;
        end
    end
end

