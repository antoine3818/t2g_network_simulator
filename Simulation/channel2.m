function outputData = channel2(rxInfo, txData)



    waveformInfo = nrOFDMInfo(51,30000/1e3);

    persistent  tdl 
    if isempty(tdl)
        tdl= nrTDLChannel;
        tdl.DelayProfile = 'TDL-C';
        tdl.DelaySpread = 300e-9;
        tdl.NumReceiveAntennas =4;
        tdl.NumTransmitAntennas =4;
        speed = 500/3.6 ;
      
        tdl.MaximumDopplerShift = speed/physconst('LightSpeed')*txData.CenterFrequency ;
        tdl.SampleRate = waveformInfo.SampleRate;
    end

    pathFilter = getPathFilters(tdl).';
    chInfo = info(tdl)   ;

    maxChannelDelayMatrix = ceil(max(chInfo.PathDelays*waveformInfo.SampleRate)) + chInfo.ChannelFilterDelay;

    % Copy input data to output
    outputData = txData;

    % Calculate distance and path loss
    distance = norm(txData.TransmitterPosition - rxInfo.Position);
    lambda = physconst('LightSpeed') / txData.CenterFrequency;
    pathLoss = fspl(distance, lambda);
    outputData.Power = outputData.Power - pathLoss;

    


    release(tdl) 
    tdl.InitialTime = outputData.StartTime;   
    if outputData.Abstraction == 0                           % Full PHY
        rxWaveform = [txData.Data; zeros(maxChannelDelayMatrix, ...
            size(txData.Data,2))];

        [outputData.Data,outputData.Metadata.Channel.PathGains, outputData.Metadata.Channel.SampleTimes] = ...
            tdl(rxWaveform);

        outputData.Data = outputData.Data.*db2mag(-pathLoss);
        outputData.Duration = outputData.Duration + (1/outputData.SampleRate)*maxChannelDelayMatrix;
      

    else                                                     % Abstract PHY
        tdl.NumTimeSamples =  ...
            txData.Metadata.NumSamples + maxChannelDelayMatrix;
        [outputData.Metadata.Channel.PathGains, outputData.Metadata.Channel.SampleTimes] = ...
            tdl();
    end
    outputData.Metadata.Channel.PathFilters = pathFilter;
    
end

