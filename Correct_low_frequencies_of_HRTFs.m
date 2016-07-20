% CORRECT LOW FREQUENCIES OF HRTFS
%
% This script corrects the frequency response of HRTFs at low frequencies by
% setting the magnitude response to a constant value and extrapolating the phase
% linearly. As the group delay at low frequencies is now low enough, the HRTFs can
% be truncated to 512 samples.
% Confer the following paper for details on this method:
%
%  - Xie, B. (2009): On the low frequency characteristics of head-related
%    transfer functions. Chinese J. Acoust. 28(2), pp. 1-13
%
% The HRTFs are saved as a SOFA file and you need the SOFA API for Matlab/Octave
% to run this script. SOFA is a file format for storing spatially oriented
% acoustic data standardised by the Audio Engineering Society (AES) as AES69-2015.
% See http://www.sofaconventions.org for details. 
%
%
% This script has been written in Matlab version R2015a.

%********************************************************************************
% Copyright (c) 2016 Vera Erbes                                                 *
%                                                                               *
% Permission is hereby granted, free of charge, to any person obtaining a copy  *
% of this software and associated documentation files (the "Software"), to deal *
% in the Software without restriction, including without limitation the rights  *
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     *
% copies of the Software, and to permit persons to whom the Software is         *
% furnished to do so, subject to the following conditions:                      *
%                                                                               *
% The above copyright notice and this permission notice shall be included in    *
% all copies or substantial portions of the Software.                           *
%                                                                               *
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    *
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      *
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   *
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        *
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, *
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN     *
% THE SOFTWARE.                                                                 *
%********************************************************************************

close all; clear; clc

SOFAstart

%% Load original data
filename = 'QU_KEMAR_anechoic_3m.sofa';
% Download data if necessary
if ~exist(filename,'file')
    disp(['Downloading ' filename]);
    url = ['https://github.com/sfstoolbox/data/blob/master/HRTFs/', ...
           'QU_KEMAR_anechoic_3m.sofa?raw=true'];
    websave(filename,url);
end
% load HRTFs
original = SOFAload(filename);

%% Prerequisites
fs = original.Data.SamplingRate; %sampling frequency
N = original.API.N; %original length of HRTFs
f = (0:N-1)/N*fs; %frequency vector
%find indices for frequency range used to correct low frequency data
indl = find(f>=100,1,'first');
indh = find(f<300,1,'last');

new_length = 512; %new length of HRTFs

%% Correct and truncate impulse responses
new_irs = zeros(original.API.M,original.API.R,new_length); %preallocation
for n = 1:original.API.M %all azimuths
    for nn = 1:original.API.R %left and right channel
        ir = squeeze(original.Data.IR(n,nn,:));
        HRTF = fft(ir);
        
        %correct magnitude
        %calculate mean magnitude from 100-300 Hz
        mag_mean = mean(abs(HRTF(indl:indh)));
        %replace low frequency magnitude with mean value
        HRTF_mag_corrected = abs(HRTF);
        HRTF_mag_corrected(1:indl-1) = mag_mean;
        
        %correct phase
        HRTF_phase = unwrap(angle(HRTF));
        HRTF_phase_corrected = HRTF_phase;
        %extrapolate phase from data between 100-300 Hz
        HRTF_phase_corrected(1:indl-1) = interp1(f(indl:indh),...
            HRTF_phase(indl:indh),f(1:indl-1),'linear','extrap');
        
        %compose complex spectrum
        HRTF_corrected = HRTF_mag_corrected.*exp(1i*HRTF_phase_corrected);
        %calculate impulse response
        ir_corrected = ifft(HRTF_corrected,'symmetric');
        
        %truncate corrected impulse response
        ir_corrected_trunc = ir_corrected(1:new_length);
        %apply window at end of impulse response
        Nwin = 32; %window length
        win = blackmanharris(Nwin);
        win = win(Nwin/2+1:end); %use only second half
        ir_corrected_trunc(new_length-Nwin/2+1:end) = ...
            ir_corrected_trunc(new_length-Nwin/2+1:end).*win;
        %apply window at start of impulse response
        Nwin = 64; %window length
        win = blackmanharris(Nwin);
        win = win(1:Nwin/2); %use only first half
        ir_corrected_trunc(1:Nwin/2) = ir_corrected_trunc(1:Nwin/2).*win;

        %save new impulse response for SOFA data array
        new_irs(n,nn,:) = ir_corrected_trunc;
        
        %plot HRTFs only for one example azimuth
        azimuth_example = 0; %in degree
        if original.SourcePosition(n,1) == azimuth_example
            HRTF_corrected_trunc = fft(ir_corrected_trunc,2048);
            gd = grpdelay(ir,1,2048,'whole',fs);
            gd_corrected_trunc = grpdelay(ir_corrected_trunc,1,2048,'whole',fs);
            if nn == 1
                side = 'left';
            elseif nn == 2
                side = 'right';
            end
            figure('name',[num2str(original.SourcePosition(n,1)) '° azimuth, '...
                side])
            subplot(2,2,1)
                plot(db(abs(ir))), hold on
                plot(db(abs(ir_corrected_trunc)))
                grid
                axis([0 2048 -140 0])
                xlabel('time in samples'),
                ylabel('impulse response magnitude in dB')
                legend('original','corr. + trunc.')
            subplot(2,2,3)
                semilogx(f,db(abs(HRTF))), hold on
                semilogx(f,db(abs(HRTF_corrected_trunc)))
                grid
                axis([10 fs/2 -60 10])
                xlabel('frequency in Hz'), ylabel('magnitude response in dB')
            subplot(2,2,2)
                semilogx(f,angle(HRTF)), hold on
                semilogx(f,angle(HRTF_corrected_trunc))
                grid
                axis([10 fs/2 -pi pi])
                xlabel('frequency in Hz'), ylabel('phase in rad')
            subplot(2,2,4)
                semilogx(f,gd), hold on
                semilogx(f,gd_corrected_trunc)
                grid
                axis([10 fs/2 -500 3500])
                xlabel('frequency in Hz'), ylabel('group delay in samples')
        end
    end
end

%% Update SOFA fields
new = original; %copy SOFA struct

%replace data array
new.Data.IR = new_irs;

%round source positions
new.SourcePosition = round(new.SourcePosition);

%extend author contact
new.GLOBAL_AuthorContact = [new.GLOBAL_AuthorContact ...
    ', vera.erbes@uni-rostock.de'];

%update history information
new.GLOBAL_History = [new.GLOBAL_History ...
    'Low frequency correction and truncation to ' num2str(new_length) ' samples'];

% %license
% new.GLOBAL_License = '';

%replace URL
new.GLOBAL_Origin = '';

%update title
new.GLOBAL_Title = 'HRTF low frequency corrected';

%% Save new SOFA file
SOFAsave('KEMAR_HRTFs_lfcorr.sofa',new,1)