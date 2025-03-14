function ColourPalette = CommonColourPalette(RewardProbCategories)
% colour palette for events (suitable for most colourblind people)
ColourPalette.Incorrect = [254, 60, 60]/255; % scarlet
ColourPalette.General = [31, 54, 104]/255; % denim
ColourPalette.Rewarded = [0, 162, 254]/255; % azure

ColourPalette.NotBaited = [26, 255, 26]/255; % neon_green
ColourPalette.SkippedBaited = [168, 12, 180]/255; % neon_purple

ColourPalette.Left = [225, 190 106]/255; % sand
ColourPalette.Right = [64, 176, 166]/255; % turquoise
ColourPalette.LeftRight = [ColourPalette.Left; ColourPalette.Right]; % LRPalette

ColourPalette.Explore = [230, 97, 0]/255; % carrot
ColourPalette.Exploit = [93, 58, 155]/255; % violet

ColourPalette.Session = [1, 1, 1] * 0.85;
ColourPalette.Pooled = [0, 0, 0];

ColourPalette.RL = [0, 0, 255]/255; % blue
ColourPalette.GLM = [255, 0, 0]/255; % red

% colour palette for cues: (1- P(r)) * 128 + 127
% P(0) = white; P(1) = smoky gray
if nargin >= 1
    ColourPalette.RewardProbLight = ((1 - RewardProbCategories) * [128 128 128] + 127)/255; % CuedPalette 
    ColourPalette.RewardProbDark = ((1 - RewardProbCategories) * [128 128 128] + 96)/255;
end

end

