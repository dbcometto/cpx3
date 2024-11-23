function [beams] = beamform_test(data2, NumBeams, FrameSize, num_elements)
%    Convert the data from the 4 channels into 21 "beams" using simple
%    delay-sum beamforming
%
%
% by VICTOR CHEN
%
%   INPUTS
%       data2:          4 channel data from the phased array
%       NumBeams:       number of beams = ratio of Fs/F + 1
%                       we upsampled to increase the number of beams
%       FrameSize:      Number of samples in each channel (after upsample)
%       num_elements:   Number of channels in the phased array
%
%   OUTPUTS
%       beams:           beamformed data
%

    % Initialize n and k terms
    k_max = floor(NumBeams/2);
    k_range = -k_max:k_max;   % -10 to 10

    n_min = 3*k_max + 1;
    n_max = FrameSize - 3*k_max;
    n_range = n_min:n_max; % n is from 31 to 3970

    kk_range = k_range + k_max + 1;

    % Create array to hold beam vals
    beams = zeros(21, FrameSize);   % Initialize beamformed array

    % Precompute offsets for all sensors and k values
    k_offsets = [k_range', 2 * k_range', 3 * k_range']; % Precompute offsets

    % Beamforming computation
    for k = 1:length(k_range)
        kk = kk_range(k);     % Shifted index for kk (1 to 21)

        if kk < 11
            % Apply beamforming equation (vectorized for all valid n)
            temp = -k_offsets(k,3);
            beams(kk, n_range) = data2(n_range + temp, 1) + ...
                                 data2(n_range + temp + k_offsets(k,1), 2) + ...
                                 data2(n_range + temp + k_offsets(k,2), 3) + ...
                                 data2(n_range + temp + k_offsets(k,3), 4);

        else
            % Apply beamforming equation (vectorized for all valid n)
            beams(kk, n_range) = data2(n_range, 1) + ...
                                 data2(n_range + k_offsets(k,1), 2) + ...
                                 data2(n_range + k_offsets(k,2), 3) + ...
                                 data2(n_range + k_offsets(k,3), 4);
        end
    end
end

