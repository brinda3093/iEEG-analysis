% Edit date: 1-8-2020; Brinda Sevak
% Edit - Changed the 1st frequency band from 1-4 to 2-4 Hz for SWA type 1
% and SWA type 2

function [fft_bands, freq_bands, freq] = csc_calculate_freq_bands_last(fft_all, freq_range, options)
% concatenates all the ffts from the input into one large variable

% define the frequency bands of interest

    switch options.fft_bands
        case 1
    freq_bands =   [2,   4; 
                    4,   8;
                    1,   4;
                    8,  12;
                    12, 16;
                    15, 25;
                    25, 40]; % changed by Anna (in wispic 1-4 instead of 0.3-4, and 20-30 instead of 20-30)
                
freq = {'SW1','SW2','Delta(1-4)','Theta (4-8 Hz)','Alpha (8-12 Hz)','Spindles (12-16 Hz)','HiBeta (15-25 Hz)','Gamma (25-40 Hz)'};
                
%         case 2
%     freq_bands =   [0.5, 4 ;
%                     5,   8 ;
%                     8,  12 ;
%                     12, 16 ;
%                     16, 25]; 
    end % commented by Anna


% calculate the power of the specified frequencies
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

% pre-allocate the range
number_bands = size(freq_bands, 1);
fft_bands = zeros(size(fft_all, 1), size(fft_all, 3), number_bands);

% loop over each band and average the frequencies
for b = 1:number_bands
   range = freq_range >= freq_bands(b, 1)  &  freq_range <= freq_bands(b, 2);
   fft_bands(:, :, b) = squeeze(nanmean(fft_all(:, range, :), 2));
end

% save to an external file if requested
if options.save_file == 0 % ==0 added by Anna
    
    % append to the already saved file
    save(fullfile(options.save_path, options.save_name)', 'fft_bands', 'freq_bands', 'freq', '-append');
    
end

end
