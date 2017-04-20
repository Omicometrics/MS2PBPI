%% Assign regions to the input X and predict intesities...
function t = ms2predict(X,pathname,c,modtype)
%
% This function is to assign regions to the input X and predict intensity
% for it. However, this output "t" is only the prediction of input X that
% is derived from one fragmentation pathway.
% The input "pathname" denotes the name of fragmentation pathways listed
% below: 
% 
% "c" is the charge state of precursors which generate the input
% fragmentation variable X.
% "modtype" denotes the type of modifications in peptide sequence with
% values 'n', 'm' and 'p' that represent no modification, oxidation of
% metheonine and other modifications, and phosphorylation respectively.
%

% You must put the model folders in the same directory of the functions.
dirName = cd;

% Specifying the directory of models
switch modtype
    case 'n'
        modelDir = [dirName '\models\mods\' pathname '_c' c '\'];
    case 'm'
        modelDir = [dirName '\models\mods_m\' pathname '_c' c '\'];
    case 'p'
        modelDir = [dirName '\models\mods_p\' pathname '_c' c '\'];
end

if ~isdir(modelDir)
    warning('MATLAB:fileAbsent',...
        'No such directory existed, predicted value of 0 is output...');
    t = 0;
    return;
end

load([modelDir 'treemod.mat']); % Tree models
load([modelDir 'modinfo.mat']); % model information

% Assign regions to the input variable set X
yRegion = regionassign(X,evStr,nID,ib,pv);
numSample = size(X,1);
t = zeros(numSample,1);

% Predicting intensity for each region
regionIdx = unique(yRegion);
for ii = 1:numel(regionIdx)
    currSamidx = yRegion==regionIdx(ii);
    if m(regionIdx(ii)) == 0
        t(currSamidx) = v{regionIdx(ii)};
    else
        load([modelDir v{regionIdx(ii)} '.mat']);
        ypre = sgbtreepredict(X(currSamidx,:),params,evStr,nID,ib,pv);
        t(currSamidx) = ypre;
    end
end