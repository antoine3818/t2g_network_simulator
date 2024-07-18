

function outputData = channel(rxInfo,txData)


outputData =  txData;

distance = norm(txData.TransmitterPosition - rxInfo.Position);    % Distance between the transmitter and the receiver
lambda = physconst('LightSpeed')/txData.CenterFrequency;          % Wavelength calculation using the speed of light and center frequency.
pathLoss = fspl(distance, lambda);                                % Free space path loss
outputData.Power = outputData.Power - pathLoss;                   % Adjust the power of the output data based on the

% Retrieve the channel model information
channelModelInfo = getSetChannelModel();

% Obtain the path filter from the channel model information
pathFilter = channelModelInfo.PathFilter;

% Obtain the maximum channel delay matrix from the channel model information
maxChannelDelayMatrix = channelModelInfo.MaxChannelDelayMatrix;

if ~isempty(channelModelInfo.TDLChannels{txData.TransmitterID,rxInfo.ID})
   
    % A channel exists between the transmitter and receiver nodes
    release(channelModelInfo.TDLChannels{txData.TransmitterID,rxInfo.ID})
    channelModelInfo.TDLChannels{txData.TransmitterID,rxInfo.ID}.InitialTime = outputData.StartTime;
    
    if outputData.Abstraction == 0                           % Full PHY
        rxWaveform = [txData.Data; zeros(maxChannelDelayMatrix(txData.TransmitterID,rxInfo.ID), ...
            size(txData.Data,2))];

        [outputData.Data,outputData.Metadata.Channel.PathGains, outputData.Metadata.Channel.SampleTimes] = ...
            channelModelInfo.TDLChannels{txData.TransmitterID,rxInfo.ID}(rxWaveform);

        outputData.Data = outputData.Data.*db2mag(-pathLoss);
        outputData.Duration = outputData.Duration + (1/outputData.SampleRate)*maxChannelDelayMatrix(txData.TransmitterID,rxInfo.ID);


    else                                                     % Abstract PHY
        channelModelInfo.TDLChannels{txData.TransmitterID,rxInfo.ID}.NumTimeSamples =  ...
            txData.Metadata.NumSamples + maxChannelDelayMatrix(txData.TransmitterID,rxInfo.ID);
        [outputData.Metadata.Channel.PathGains, outputData.Metadata.Channel.SampleTimes] = ...
            channelModelInfo.TDLChannels{txData.TransmitterID,rxInfo.ID}();
    end
    outputData.Metadata.Channel.PathFilters = ...
        pathFilter{txData.TransmitterID,rxInfo.ID};
    outputData.Metadata.Channel.PathDelays = [0 3.0000e-07 6.0000e-07] ;

else
    % Set default values for channel parameters
    outputData.Metadata.Channel.PathGains = ...
        permute(ones(outputData.NumTransmitAntennas,rxInfo.NumReceiveAntennas),[3 4 1 2])/ ...
        sqrt(rxInfo.NumReceiveAntennas);
    outputData.Metadata.Channel.PathFilters = 1;
    outputData.Metadata.Channel.SampleTimes = 0;

    if outputData.Abstraction == 0                            % Full PHY
        outputData.Data = outputData.Data.*db2mag(-pathLoss);
        numTxAnts = outputData.NumTransmitAntennas;
        numRxAnts = rxInfo.NumReceiveAntennas;
        H = fft(eye(max([numTxAnts numRxAnts])));
        H = H(1:numTxAnts,1:numRxAnts);
        H = H/norm(H);
        outputData.Data = txData.Data*H;                      % Apply channel to the waveform
    end
end
end

