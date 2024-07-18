function channelModelInfo = getSetChannelModel(varargin)
    persistent channelModels;
    if nargin == 1 % Set
        channelModels = varargin{1};
    end
    channelModelInfo = channelModels;
end