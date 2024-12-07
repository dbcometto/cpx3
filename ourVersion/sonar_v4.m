function cpx3_sonar_v4
% by the 3 Engineers Plus Ben team
% sends 10KHz pulses, received by 4 channel array, does sonar processing, and plots
% in polar format
% Also can work on canned test data


% Default Parameters
format LONG        % increases the percision of numbers printed on the command line
continuous = 1.0;  % 0 means only do one image; >0 means continuously create sonar images
live = 0.0;        % 0 means artficial data; >0 means live data from A/D board
show_images = 0; % 0 means do not show intermediate stages images; 1 means show all images
persist = 0.75;     % ranges from 0.0 for no persistence to 1.0 for maximum persistence

c = 1136;           %speed of sound in air in ft/sec
SampleRate = 100000;
frequency = 10000;  % of pulse
samples_per_cycle = SampleRate / frequency;
N = 2000;           % number of samples taken per sensor; range = about 10 feet
NUpsampled = 2*N;   % number of samples after upsampling
FrameSize = N;
upsample = 2;       % later upsampling by 2
num_elements = 4;   % number of elements in array
cycles = 6;         % number of cycles in transmitted pulse
Fs_out = 44100;     % sample frequency of the output DAC
max_range = 100;    % range of the output display in # of pixels
min_range = 2;      % range blanked on the output display in # of pixels
pulse_length = cycles*samples_per_cycle; % pulse length in number of samples
top    = 0.9;     % top threshold for contrast enhancement (window/level) (must be between 0 and 1) 
bottom = 0.6;     % bottom thrseshold for contrast enhancement (window/level) (must be between 0 and 1)
colormap(gray);   % grayscale images, not color
image_rows = 100;
image_col = 201;

% NumBeams = 21; 
NumBeams = (samples_per_cycle*upsample) + 1;  %should be 21 for default case
beamShift = round(NumBeams/2); % should be 11 for default case
USampleRate = 2*SampleRate;

% CREATE DATA ARRAYS NEEDED FOR FUNCTIONS
[ind_bkn, ind_bk1n, ind_bkn1, ind_bk1n1, BMAM, BMA, BAM, BA, pixel_per_foot_row, pixels_per_foot_col] = scan_conversion_precompute(frequency, USampleRate, c, NUpsampled, NumBeams, beamShift, image_rows, image_col);
track_arrays = [];

[ind_bkn, ind_bk1n, ind_bkn1, ind_bk1n1, BMAM, BMA, BAM, BA, feet_per_pixel_row, feet_per_pixel_col] = scan_conversion_precompute(frequency, USampleRate, c, NUpsampled, NumBeams, beamShift, image_rows, image_col);

% PREALLOCATED ARRAYS
data2 = zeros(N*upsample,num_elements); % preallocate and create zeros vector for upsampling (designed for upsample = 2)
beams = zeros(NumBeams, N*upsample);
demod_I = zeros(NumBeams, N*upsample);  % Preallocate for Demod Mix
demod_Q = zeros(NumBeams, N*upsample);
demod_I_LPF = zeros(NumBeams, N*upsample);  % Preallocate for Demod LPF
demod_Q_LPF = zeros(NumBeams, N*upsample);  % Preallocate for Demod LPF
size(demod_I)
Mag_image = zeros(N*upsample, NumBeams);
persist_image = zeros(max_range,(2*max_range + 1));
SF = zeros(4,1);

% PRECOMPUTE LOOK-UP TABLES TO SPEED UP FUNCTIONS
denoisers = load("denoise_fils.mat","-mat");
QuadIncr = pi*2*(frequency/(SampleRate*upsample));
table_len = FrameSize*upsample;
cos_table = zeros(1, table_len);
sin_table = zeros(1, table_len);
for i = 1:table_len
    cos_table(i) = cos(i*QuadIncr);
    sin_table(i) = sin(i*QuadIncr);
end


%time-gain computation 
sampleIndex = [1:1:2000;1:1:2000;1:1:2000;1:1:2000]';
timeGainValue = (1 + (sampleIndex.*c)/SampleRate);

load lpf_upsample.bin -mat % load upsampling filer coefficients



%blackman filter coefficients
% For LPF using blackman rectangular filter
WindowLength = 2*samples_per_cycle*upsample; % 2 cycles
a=1:WindowLength;
blackman(a) = 0.42 - 0.5*cos(2*pi*(a-1)/(WindowLength-1)) + 0.08*cos(4*pi*(a-1)/(WindowLength-1));
filter_coef = blackman;

% initialize A/D board, if necessary
if live > 0
    c6x_daq('Init', 'ths1206.out');     % use TI THS1206 EVM
    c6x_daq('FrameSize', N);
    c6x_daq('QueueSize', 2*N);
    SampleRate = c6x_daq('SampleRate', SampleRate);
    numChannels = c6x_daq('NumChannels', 4);
    c6x_daq('TriggerMode', 'Auto');
    c6x_daq('TriggerSlope', '+');
    c6x_daq('TriggerValue', 0.1);
    c6x_daq('TriggerChannel', 0);
    c6x_daq('GetSettings');
    c6x_daq('Version');
end

if live > 0
    % send out test pulse to warm up speaker
    pulse = create_pulse(cycles, frequency, Fs_out);
    sound(pulse,Fs_out);
    pause(1);
end


iii = 1;
game_on = 1;
while game_on > 0

    % we combined a stage1a & 1b, moved stage3 to before
    Stage1b_cal_array_time = 0;
    
    % STAGE 0: TRANSMIT PULSE
    time1 = tic;
    
    if live > 0
        %create pulse to transmit
        pulse = create_pulse(cycles, frequency, Fs_out);
        %send sound out PC speaker port 
        sound(pulse,Fs_out);
        %receive data
        data = c6x_daq('GetFrame');
        time2 = toc(time1);
        if (show_images == 1.0)  % plot pulse
            figure(1);
            plot(pulse);
            title('Transmitted Pulse');
            xlabel('sample number');
            % add labels
        end
    else
        % use canned data
        if continuous == 0
           dataz = load('test_data');
           data = dataz.data;
        else
           dataz = load('test_data3');
           data(1:N,1:num_elements) = dataz.data_frames(1:N,1:num_elements,iii);
        end
        time2 = toc(time1);
    end
    Stage0_transmit_time = time2
    total_time = Stage0_transmit_time;
    total_time = 0;  % don't include transmit pulse since using canned data
    
    % plot received data
    if (show_images == 1.0)
        figure(2);
        P1=plot(data(:,1),'k');
        hold on
        P2=plot(data(:,2),'r');
        P3=plot(data(:,3),'g');
        P4=plot(data(:,4),'b');
        hold off
        legend('CH 1', 'CH 2', 'CH 3', 'CH 4');
        title('Raw data received');
        xlabel('sample number');
    end

   % STAGE 3: Noise and Bias Removal
   %    This is a good place to add a noise removal LPF,
   %    especially to band limit the signal before upsampling
 
   % First plot spectrum of original data
   % (use simple FFT, rect window)
   % only plot channel 1
    if (show_images == 1.0)
        figure(6);
        a = 1:N;
        x = a*(SampleRate/N);
        plot(x, abs(fft(data(:,1))));
        title('channel 1 spectrum (not shifted) before noise smoothing');
        xlabel('frequency (Hz)');
    end
   
    time1 = tic;
 
    [data] = noise_remove_LPF(data, denoisers); 

    time2 = toc(time1);
    Stage3_noiseLPF_time = time2
    total_time = total_time + Stage3_noiseLPF_time;
    
    % plot modified data
    if (show_images == 1.0)
        figure(7);
        P1=plot(data(:,1),'k');
        hold on
        P2=plot(data(:,2),'r');
        P3=plot(data(:,3),'g');
        P4=plot(data(:,4),'b');
        hold off
        legend('CH 1', 'CH 2', 'CH 3', 'CH 4');
        title('data after noise smoothing');
        xlabel('sample number');
    end
   % Plot spectrum of data after LPF
   % (use simple FFT, rect window)
   % only plot channel 1
    if (show_images == 1.0)
        figure(8);
        a = 1:N;
        x = a*(SampleRate/N);
        plot(x, abs(fft(data(:,1))));
        title('channel 1 spectrum (not shifted) after noise smoothing');
        xlabel('frequency (Hz)');
    end
    
    % STAGE 1: Calibration, remove transimitted "bang" and find start of data
    %          - move data up to start of bang
    %          - blank out transmitted pulse (could impact calibration)
    %          - blank out noise after reflection return
    time1 = tic;
 
    data = cal_start(data, N, pulse_length); 
    
    time2 = toc(time1);
    Stage1a_cal_start_time = time2 
    total_time = total_time + Stage1a_cal_start_time;
    
    % plot modified data
    if (show_images == 1.0)
        figure(3);
        P1=plot(data(:,1),'k');
        hold on
        P2=plot(data(:,2),'r');
        P3=plot(data(:,3),'g');
        P4=plot(data(:,4),'b');
        hold off
        legend('CH 1', 'CH 2', 'CH 3', 'CH 4');
        title('data after blanking and shifting');
        xlabel('sample number');
    end
   
    
  % STAGE 2: Time Gain Compensation
  %  adjust for annenuation of signal over distance
    time1 = tic;
 
    [data] = time_gain_compensation(data, timeGainValue);
   % [data] = time_gain_compensation(data, N); 
    
    time2 = toc(time1);
    Stage2_TGC_time = time2
    total_time = total_time + Stage2_TGC_time;
    
    % plot modified data
    if (show_images == 1.0)
        figure(5);
        P1=plot(data(:,1),'k');
        hold on
        P2=plot(data(:,2),'r');
        P3=plot(data(:,3),'g');
        P4=plot(data(:,4),'b');
        hold off
        legend('CH 1', 'CH 2', 'CH 3', 'CH 4');
        title('data after time gain compensation');
        xlabel('sample number');
    end
    

   % STAGE 4: Upsampling
   %    Upsampling is done to increase the number of beams after
   %    beamforming (assuming simple delay-add beamforming)
   %    Assuming signal is bandlimited (by stage 5), upsampling
   %    can be done my (1) padding with zeroes between samples
   %    (2) then low pass filter the signal "images" created
   %
   %    Note: Dr York's routine is fixed for upsamping by "2".
   %    Call routine multiple times to upsample by 4, 8, 16, etc.
    time1 = tic;
 
    [data2] = upsampling(data,num_upsample,data2,NUpsampled); 

    time2 = toc(time1);
    Stage4_upsample_time = time2
    total_time = total_time + Stage4_upsample_time;
    
    % plot modified data
    if (show_images == 1.0)
        figure(9);
        P1=plot(data2(:,1),'k');
        hold on
        P2=plot(data2(:,2),'r');
        P3=plot(data2(:,3),'g');
        P4=plot(data2(:,4),'b');
        hold off
        legend('CH 1', 'CH 2', 'CH 3', 'CH 4');
        title('data after upsample');
        xlabel('sample number');
    end
   % Plot spectrum of data after upsample
   % (use simple FFT, rect window)
   % only plot channel 1
    if (show_images == 1.0)
        figure(10);
        a = 1:((N*upsample));
        x = a*(SampleRate/N);
        plot(x, abs(fft(data2(:,1))));
        title('channel 1 spectrum (not shifted) after upsample');
        xlabel('frequency (Hz)');
    end
   
    %break


    % STAGE 5:  Beamforming
    % Convert the data from the 4 channels into 21 "beams" using simple
    % delay-sum beamforming
    time1 = tic;
 
    [beams] = beamform(data2, NumBeams, FrameSize*upsample, num_elements);

    time2 = toc(time1);
    Stage5_beamform_time = time2
    total_time = total_time + Stage5_beamform_time;
    
    % plot modified data
    if (show_images == 1.0)
        figure(11);
        P1=plot(beams(11,:),'k');
        hold on
        P2=plot(beams(15,:),'r');
        hold off
        legend('broadside', 'off angle');
        title('Beams k=0 (broadside) and k=4 (off angle)');
        xlabel('sample number');
    end
   % Plot spectrum of broadside beam 
   % (use simple FFT, rect window)
   % only plot channel 1
    if (show_images == 1.0)
        figure(12);

        % a = 1:3984;
        a = 1:(FrameSize*upsample);
        %a = 1:beamlength;
        x = a*(SampleRate/N);
        %length(x)
        %length(abs(fft(beams(11,:))))
        plot(x, abs(fft(beams(11,:))));
        title('Spectrum of Broadside Beam');
        xlabel('frequency (Hz)');
    end
   
    % Plot image of beamformed data
    if (show_images == 1.0)
       figure(13);
       colormap(gray);
       imagesc(transpose(beams));
       title('Image of Beamformed Data');
       xlabel('Beams (each corresponds to a different angle)');
       ylabel('Distance (really sample number)');
   end
    
    
    
   % STAGE 6: Demodulate (get signal envelope or echo image)
   % This is done in 3 substages.  We use a quadrature detector (I + jQ)
   % because we cannot guarrentee to be in phase with the pulse.
   % STAGE 6a: Quad Demod Mixing
   % Demodulate the echo image by finding the envelope of the carrier pulse
   % First, multiply by the carrier (cos + jsin)
   % This moves the echo signal to DC and makes another copy at twice the
   % carrier frequency.
   % STAGE 6b: Quad Demod LPF
   %  LPF to remove harmonic at twice f, keeping DC component
   %  Convolves data with Filter_coef of length WindowLength
   % STAGE 6c: Find magnitude (echo image) from imaginary
   %           Find magnitude of I + jQ
   
   % TEST DATA - UNCOMMENT TO USE
   % Using a sinwave for all channels, we can test the beamform on various
   % frequencies to see waves in our final beamform image.
   %x = 2;    % adjust x for different frequencies
   %t = 1:4000;
   
   %size(beams)
   %NumBeams
   %channel = sin(t/x);
   %beams = [channel', channel', channel', channel', channel', channel', channel', channel', channel', channel', channel', channel', channel', channel', channel', channel', channel', channel', channel', channel', channel'];
   %size(beams.')
   
   
   time1 = tic;

   [demod_I, demod_Q] = quad_demod_mix(beams, NumBeams, cos_table, sin_table);

   %save("quad_test.mat", "beams", "demod_Q", "demod_I", "frequency")

   [demod_I_LPF, demod_Q_LPF] = quad_demod_LPF(demod_I, demod_Q, NumBeams, filter_coef, demod_I_LPF, demod_Q_LPF);
   [Mag_image] = magnitude(demod_I_LPF, demod_Q_LPF);

   time2 = toc(time1);
   Stage6_demod_time = time2
   total_time = total_time + Stage6_demod_time;
   
    % plot modified data
    if (show_images == 1.0)
        figure(14);
        P1=plot(demod_I(11,:),'k');
        hold on
        P2=plot(demod_I(15,:),'r');
        hold off
        legend('broadside', 'off angle');
        title('After Mixing with cosine: component of Beams k=0 (broadside) and k=4 (off angle)');
        xlabel('sample number');
    end
   % Plot spectrum of broadside beam 
   % (use simple FFT, rect window)
   % only plot channel 1
    if (show_images == 1.0)
        figure(15);

        % a = 1:3984;
        a = 1:(FrameSize*upsample);
        %a = 1:beamlength;
        x = a*(SampleRate/N);
        %length(x)
        %length(abs(fft(beams(11,:))))
        plot(x, abs(fft(demod_I(11,:))));
        title('Spectrum of Broadside component of Beam after Mixing: ');
        xlabel('frequency (Hz)');
    end
   
    % Plot image of beamformed data
    if (show_images == 1.0)
       figure(16);
       colormap(gray);
       imagesc(transpose(demod_I));
       title('Image of I component of Beamformed Data After Mixing: ');
       xlabel('Beams (each corresponds to a different angle)');
       ylabel('Distance (really sample number)');
   end
       % plot modified data

    if (show_images == 1.0)
        figure(17);
        P1=plot(demod_I_LPF(11,:),'k');
        hold on
        P2=plot(demod_I_LPF(15,:),'r');
        hold off
        legend('broadside', 'off angle');
        title('After Quad MLPF: I component of Beams k=0 (broadside) and k=4 (off angle)');
        xlabel('sample number');
    end
   % Plot spectrum of broadside beam 
   % (use simple FFT, rect window)
   % only plot channel 1
    if (show_images == 1.0)
        figure(18);

        % a = 1:3984;
        a = 1:(FrameSize*upsample);
        %a = 1:beamlength;
        x = a*(SampleRate/N);
        %length(x)
        %length(abs(fft(beams(11,:))))
        plot(x, abs(fft(demod_I_LPF(11,:))));
        title('Spectrum of Broadside I component of Beam after Quad LPF: ');
        xlabel('frequency (Hz)');
    end
   
    % Plot image of beamformed data
    if (show_images == 1.0)
       figure(19);
       colormap(gray);
       imagesc(transpose(demod_I_LPF));
       title('Image of I component of Beamformed Data After Quad LPF: ');
       xlabel('Beams (each corresponds to a different angle)');
       ylabel('Distance (really sample number)');
   end
   
    if (show_images == 1.0)
        figure(20);
        P1=plot(Mag_image(11,:),'k');
        hold on
        P2=plot(Mag_image(15,:),'r');
        hold off
        legend('broadside', 'off angle');
        title('After demodulation: Beams k=0 (broadside) and k=4 (off angle)');
        xlabel('sample number');
    end
   % Plot spectrum of broadside beam 
   % (use simple FFT, rect window)
   % only plot channel 1
    if (show_images == 1.0)
        figure(21);

        % a = 1:3984;
        a = 1:(FrameSize*upsample);
        %a = 1:beamlength;
        x = a*(SampleRate/N);
        %length(x)
        %length(abs(fft(beams(11,:))))
        plot(x, abs(fft(Mag_image(11,:))));
        title('Spectrum of Broadside Beam after Demodulation: ');
        xlabel('frequency (Hz)');
    end
   
    % Plot image of beamformed data
    if (show_images == 1.0)
       figure(22);
       colormap(gray);
       imagesc(transpose(Mag_image));
       title('Echo Image After Demodulation: ');
       xlabel('Beams (each corresponds to a different angle)');
       ylabel('Distance (really sample number)');
   end
   

   % fprintf('add proper downsampling?');
   
   % STAGE 7: Scan Conversion 
   %   Properly plot the "beams" in their actual spatial geometry
   %   Requires a polar conversion (beam data is stored in a rectangular
   %   array, but must be plotted polar).  Also requires interpolating
   %   pixels that lie in between beam lines.
   %
   %   This can be greatly speeded up by precomputing lookup tables
   
   % scale pixel magnitude to match my scan conversion thresholds
   % the_max = max(max(Mag_image));
   % Mag_image = Mag_image ./ the_max;
   % Mag_image = Mag_image .* 36;

   time1 = tic;
   [sc_image] = scan_conversion(Mag_image, ind_bkn, ind_bk1n, ind_bkn1, ind_bk1n1, BMAM, BMA, BAM, BA);
   time2 = toc(time1);
   Stage7_scan_conversion_time = time2   
   total_time = total_time + Stage7_scan_conversion_time;

   % save("sc_image.mat","sc_image")
   % save("Mag_image.mat","Mag_image")
  
    % Plot scan converted echo image
    if (show_images == 1.0)
       figure(23);
       axis equal
       colormap(gray);
       imagesc(sc_image);
       title('Scan Converted Echo Image');
       xlabel('Distance Along Array');
       ylabel('Distance From Array');
   end
   
   fprintf('need to plot real distances');
   
   % STAGE 8: Image Enhancement
   %    such as speckle reduction or edge enhancement?
   
   
   
   % STAGE 9: Constrast/Brightness 
   %   This function warp's the dynamic range of the image intensity
   %   to emphasis features of interest.  There are many ways to do this
   %   such as:
   %          - logrithmic compression:  I = exp(I)
   %          - "window/level" (known as imadjust in matlab)
   %          - histeq (in matlab)

   % first look at histogram of image before contrast enhancement
   if (show_images == 1.0)
       figure(24);
       % turn image into 1-D vector
       image_vector = sc_image(1,:);
       for y=2:max_range
           image_vector = [image_vector sc_image(y,:)];
       end
       hist(image_vector,256);   
       title('Histogram of Echo Image before contrast enhanced');
       xlabel('Grey levels in image');
       ylabel('# of times each level occurred');
   end

   time1 = tic;
   
   % do "window-level"
   [sc_image] = contrast(sc_image, bottom, top);

   time2 = toc(time1);
   Stage9_contrast_time = time2 
   total_time = total_time + Stage9_contrast_time;
   
   %  look at histogram of image after contrast enhancement
   if (show_images == 1.0)
       figure(25);
       % turn image into 1-D vector
       image_vector = sc_image(1,:);
       for y=2:max_range
           image_vector = [image_vector sc_image(y,:)];
       end
       hist(image_vector,256);   
       title('Histogram of Echo Image after contrast enhanced');
       xlabel('Grey levels in image');
       ylabel('# of times each level occurred');
   end   

   % Plot scan converted echo image
   if (show_images == 1.0)
       figure(26);
       axis equal
       colormap(gray);
       imagesc(sc_image);
       title('Echo Image after contrast enhancement');
       xlabel('Distance Along Array');
       ylabel('Distance From Array');
   end
 
    
   % STAGE 10: Persistence (IIR Filter)
   %    A low pass filter over time can be done known as Persistence.
   %    Persistence is averaging the previous output image with the newly
   %    calculated image.  By changing the weight of the old vs new you
   %    can change the amount of persistence (or how long it takes old
   %    stuff to fade in the image over time).  This filter is great for
   %    remove speckle noise for slowly moving (or none moving) targets.
   %    However, it causes blurring of fast moving targets.
   
   if (persist == 0)
       
       figure(27);
       axis equal
       colormap(gray);
       imagesc(sc_image);
       title('Final Echo Image');
       xlabel('Distance Along Array');
       ylabel('Distance From Array');
       drawnow
   else
       % alternative data for persistence
       %sc_image = rand(max_range, 2*max_range + 1);
       %sc_image(25:30,47:53) = 0.5 * ones(6,7);
       
       time1 = tic;
       
       %persistence filter
       [persist_image] = persistence(sc_image, persist_image, persist);

       time2 = toc(time1);
       Stage10_persist_time = time2 
       total_time = total_time + Stage10_persist_time;
      
       figure(27);
       axis equal
       colormap(gray);
       imagesc(persist_image);
       title('Final Echo Image after persistence');
       xlabel('Distance Along Array');
       ylabel('Distance From Array');
       drawnow
   	% drawnow
    
   end
   
   % STAGE 11: Tracking and/or Velocity Estimation
   % 
   % 2-d correlates the template image
   % the template image is a block 2x2 1's to represent the most filled in
   % part of the image
   %

   % SINGLE object tracking
   
        template = ones(4,4);
        
        %returns the value to plot for tracking
        [track_row,track_col] = tracking(persist_image, template);
       
        %overlays onto figure(27)
        hold on;
        % x over target
        plot(track_col, track_row, 'x', 'MarkerSize',20, 'LineWidth', 2, 'Color', 'red')
        hold off;

        
    % Calculate Velocity
       
        
        if iii == 5
            track_arrays = [track_arrays; [track_row,track_col]];
        elseif iii == 50
            track_arrays = [track_arrays; [track_row,track_col]];
            velocity = calc_velocity(track_arrays(1, 1),track_arrays(1,2), track_arrays(2, 1),track_arrays(2,2), 2*pixel_per_foot_row, 2*pixels_per_foot_col, USampleRate, 45);
        end

       

    %TWO object tracking (does not work with velocity)
%         track_raw_data = persist_image;
%     
%         template = ones(4,4);
%         offset = 15;
%         
%         %trim data to avoid edge case
%         track_data = track_raw_data(1+offset:end-offset,1+offset:end-offset);
%         
%         %performs correlation
%         [track_row, track_col] = tracking(track_data, template);
%         
%         %plot image
%         figure(28)
%         imagesc(track_raw_data)
%         
%         %overlays onto figure(27)
%         hold on;
%         % x over target
%         plot(track_col+offset, track_row+offset, '+', 'MarkerSize',20, 'LineWidth', 2, 'Color', 'red')
%         % circle over target
%         %     plot(track_col, track_row, 'o', 'MarkerSize',50, 'LineWidth', 2, 'Color', 'red')
%         hold off;
%         
%         %remove max values for second target tracking
%         track_data(track_row-offset:track_row+offset, track_col-offset:track_col+offset) = 0;
%         
%         
%         [track_row, track_col] = tracking(track_data, template);
%         
%         
%         %overlays onto figure(27)
%         hold on;
%         % x over target
%         plot(track_col+offset, track_row+offset, 'x', 'MarkerSize',20, 'LineWidth', 2, 'Color', 'red')
%         % circle over target
%         %     plot(track_col, track_row, 'o', 'MarkerSize',50, 'LineWidth', 2, 'Color', 'red')
%         hold off;


   
    
   
   
   
   total_time
   iii = iii + 1
   if (iii == 2)
       stage1_time = Stage1a_cal_start_time + Stage1b_cal_array_time;
       stage2_time = Stage2_TGC_time;
       stage3_time = Stage3_noiseLPF_time;
       stage4_time = Stage4_upsample_time;
       stage5_time = Stage5_beamform_time;
       stage6_time = Stage6_demod_time;
       stage7_time = Stage7_scan_conversion_time;
       stage9_time = Stage9_contrast_time;
       stage10_time = Stage10_persist_time;
   end
   if (continuous == 0)
        game_on = 0;
   else
      if (iii > 2)
          temp_time = Stage1a_cal_start_time + Stage1b_cal_array_time;
          stage1_time = [stage1_time temp_time];
          stage2_time = [stage2_time Stage2_TGC_time];
          stage3_time = [stage3_time Stage3_noiseLPF_time];
          stage4_time = [stage4_time Stage4_upsample_time];
          stage5_time = [stage5_time Stage5_beamform_time];
          stage6_time = [stage6_time Stage6_demod_time];
          stage7_time = [stage7_time Stage7_scan_conversion_time];
          stage9_time = [stage9_time Stage9_contrast_time];
          stage10_time = [stage10_time Stage10_persist_time];
      end
      if iii == 142   % for test_data3
          game_on = 0.0;
      end
               
    end
    

end % while
Stage1_average = mean(stage1_time)
Stage2_average = mean(stage2_time)
Stage3_average = mean(stage3_time)
Stage4_average = mean(stage4_time)
Stage5_average = mean(stage5_time)
Stage6_average = mean(stage6_time)
Stage7_average = mean(stage7_time)
Stage9_average = mean(stage9_time)
Stage10_average = mean(stage10_time)

velocity



