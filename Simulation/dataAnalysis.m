function analysedData = dataAnalysis(data)
    % Initialization of latency array
    latency = [];
    
    % Loop through each User Equipment (UE)
    for i = 1:data.numUEs 
        % Calculate latency for each packet for the current UE
        latency = [latency, (data.packet(2:end-1, i+1) - data.packet(2:end-1, 1))]; 
        % Find the maximum latency for the current UE
        latencyMax(i) = max(latency(:, i));
        % Find the minimum latency for the current UE
        latencyMin(i) = min(latency(:, i));
        % Calculate the average latency for the current UE
        latencyAverage(i) = mean(latency(:, i));
        % Display latency statistics for the current UE
        fprintf("Node(%d): max latency %d, min latency %d, average latency %d\r\n", i, latencyMax(i), latencyMin(i), latencyAverage(i))
    end
    
    % Get and display gNB statistics
    gnbStats = statistics(data.gnb);
    fprintf("nb packet transmitted: %d\r\n", gnbStats.PHY.TransmittedPackets)
    
    % Loop through each UE to get and display packet reception and error statistics
    for i = 1:data.numUEs
        % Get statistics for the current UE
        uesStats(i) = statistics(data.ues(i));
        % Calculate the bit rate for the current UE
        bitRate(i) = (uesStats(i).App.ReceivedBytes * 8) / data.simulationTime;
        % Calculate the Packet Error Rate (PER) for the current UE
        PER(i) = uesStats(i).PHY.DecodeFailures / gnbStats.PHY.TransmittedPackets;
        % Display the number of packets received by the current UE
        fprintf("Node(%d): nb packet received: %d\n", i, uesStats(i).PHY.ReceivedPackets)
    end
    
    % Loop through each UE to display PER statistics
    for i = 1:data.numUEs
        fprintf("Node(%d): PER = %.2f%%\n", i, PER(i) * 100)
    end
    
    % Loop through each UE to display bit rate statistics
    for i = 1:data.numUEs
        fprintf("Node(%d): bitRate = %.2f Mb/s\n", i, bitRate(i) / 1e6)
    end
    
    % Store the analyzed data for each UE in the output structure
    for i = 1:data.numUEs 
        analysedData(i).Latency = latencyAverage(i);
        analysedData(i).PER = PER(i);
        analysedData(i).bitRate = bitRate(i);
    end
end