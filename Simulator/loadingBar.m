function loadingBar(~,~)
    global data 
    currentTime = data.Simulator.CurrentTime;
    totalTime  = data.simulationTime ;
    
    % Calculate the percentage of completion
    percentComplete = (currentTime / totalTime) * 100;
    
    % Determine the number of characters for the progress bar
    barLength = 50; % Length of the progress bar in characters
    numBars = floor((currentTime / totalTime) * barLength);
    
    % Create the progress bar string
    progressBarString = ['[' repmat('=', 1, numBars) repmat(' ', 1, barLength - numBars) ']'];
    
    % Display the progress bar with the percentage
    if currentTime == 0
        %fprintf('loading bar: %s %6.2f%%', progressBarString, percentComplete);
        fprintf("\n")
        fprintf('loading bar: %s %6.2f%%', progressBarString, percentComplete);
    else
        fprintf(repmat('\b', 1, barLength + 23)); % Clear the previous progress bar
        fprintf('loading bar: %s %6.2f%%', progressBarString, percentComplete);
        
    end
    
    % Newline if progress is complete
    if currentTime >= totalTime
        fprintf('\n');
    end

end
