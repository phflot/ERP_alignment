
clear;
close all;
clc;

load('erp.mat');
load('erp_reg_BM3D_1.mat');

EEG.srate = 200;
pnt_trial = size(erp);
EEG.pnts = pnt_trial(2);
EEG.trials = pnt_trial(1);
EEG.data(1,:,:) = single(erp');
EEG.times = 0:1/EEG.srate:0.7;
EEG.times = 1000*(EEG.times(1:end-1));

% wavelet parameters
min_freq = 2;
max_freq = 50;
num_frex = 30;

% other wavelet parameters
frequencies = logspace(log10(min_freq),log10(max_freq),num_frex);
time = 0:1/EEG.srate:0.7;
time = time(1:end-1);
half_of_wavelet_size = (length(time)-1)/2;

% FFT parameters (use next-power-of-2)
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
n_conv_pow2   = pow2(nextpow2(n_convolution));
wavelet_cycles= 4; 

% define baseline period (not useful with the synthetic data set, but with real ERP data)
baselinetime = [ -200 0 ]; % in ms

% convert baseline window time to indices
[~,baselineidx(1)]=min(abs(EEG.times-baselinetime(1)));
[~,baselineidx(2)]=min(abs(EEG.times-baselinetime(2)));


iterations= 10;
dbcorrect = false;

powerByTrialFreq = zeros(length(frequencies),EEG.trials);

start_time = 0; % in ms
end_time   = 700;

figure;
for input = 1:2

if input >1
    EEG.data(1,:,:) = single(Reg_BM3D_1');
end

    fft_data = fft(reshape(EEG.data(1,:,:),1,[]),n_conv_pow2);
    % fft_data = fft(reshape(EEG.data(strcmpi(chan2plot,{EEG.chanlocs.labels}),:,:),1,[]),n_conv_pow2);
    timeidx = dsearchn(EEG.times',[start_time end_time]');

    for fi=1:length(frequencies)

        % create wavelet and get its FFT
        wavelet = exp(2*1i*pi*frequencies(fi).*time) .* exp(-time.^2./(2*( wavelet_cycles /(2*pi*frequencies(fi)))^2))/frequencies(fi);
        fft_wavelet = fft(wavelet,n_conv_pow2);

        % run convolution
        convolution_result = ifft(fft_wavelet.*fft_data,n_conv_pow2) * sqrt(wavelet_cycles /(2*pi*frequencies(fi)));
        convolution_result = convolution_result(1:n_convolution);
        convolution_result = convolution_result(half_of_wavelet_size+1:end-half_of_wavelet_size);
        convolution_result = abs(reshape(convolution_result,EEG.pnts,EEG.trials)).^2; % reshape and convert to power

        % "gold standard" is average of all trials
        if dbcorrect
            template = 10*log10( bsxfun(@rdivide,mean(convolution_result,2),mean(mean(convolution_result(baselineidx(1):baselineidx(2),:),1),2)));
            template = template(timeidx(1):timeidx(2));
        else
            template = mean(convolution_result(timeidx(1):timeidx(2),:),2);
        end
        % normalize template for correlation
        template = bsxfun(@rdivide,bsxfun(@minus,template,mean(template)),std(template))';

        for iteri=1:iterations
            for triali=5:EEG.trials % start at 5 trials...

                trials2use = randsample(1:EEG.trials,triali);
                % if you don't have the stats toolbox, use the following:
                % trials2use = randperm(EEG.trials); trials2use = trials2use(1:triali);

                % compute power time series from the random selection of trials, and then normalization
                if dbcorrect
                    tempdat = 10*log10( bsxfun(@rdivide,mean(convolution_result(:,trials2use),2),mean(mean(convolution_result(baselineidx(1):baselineidx(2),trials2use),1),2)));
                    tempdat = tempdat(timeidx(1):timeidx(2));
                else
                    tempdat = mean(convolution_result(timeidx(1):timeidx(2),trials2use),2);
                end
                tempdat = bsxfun(@rdivide,bsxfun(@minus,tempdat,mean(tempdat)),std(tempdat))';

                % compute Pearson correlation. This is a super-fast
                % implementation of a Pearson correlation via least squares
                % fit. 
                powerByTrialFreq(fi,triali) = powerByTrialFreq(fi,triali) + (tempdat*tempdat')\tempdat*template';
            end
        end
    end

    powerByTrialFreq = powerByTrialFreq./iterations;

if input <2
    subplot(2,2,1)
else
    subplot(2,2,3)
end
    plot(5:EEG.trials,squeeze(powerByTrialFreq(:,5:end)))
    xlabel('Number of trials'), ylabel('Power')
    set(gca,'ylim',[-.1 1.1])
    title('Each line is a frequency band')
if input <2
    subplot(2,2,2)
else
    subplot(2,2,4)
end
    contourf(5:EEG.trials,frequencies,squeeze(powerByTrialFreq(:,5:end)),40,'linecolor','none')
    set(gca,'clim',[.6 1])
    xlabel('Number of trials'), ylabel('Frequency (Hz)')
    if dbcorrect, title('DB normalized'), else title('not dB normalized'); 
        colorbar
    end

end

