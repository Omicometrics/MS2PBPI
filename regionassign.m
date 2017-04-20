%% Predict y from stochastic gradient boosting tree (SGBTree)
function Y = regionassign(X,evalStr,nodeID,isBranch,pathVal)
%
% Predict the intensity of mass spectra according to the models constructed
% by stochastic gradient boosting tree (SGBTee)
%
% Nai-ping Dong, PolyU in HK
% Email: np.dong572@gmail.com
% 1/2/2014

m = size(X,1);
Y = zeros(m,1);

nodeIdx = ones(m,1);
sampleIdx = (1:m)';
remainNum = true;

while remainNum
    
    ib = isBranch(nodeIdx(sampleIdx(1)));
    currNodeidx = nodeIdx(sampleIdx(1));
    
    while ib
        
        currNodeidx = nodeIdx(sampleIdx(1));
        t = eval(evalStr{currNodeidx});
        
        if any(t)
            nodeIdx(sampleIdx(t)) = nodeID(currNodeidx,1);
            nodeIdx(sampleIdx(~t))= nodeID(currNodeidx,2);
            sampleIdx = sampleIdx(t);
        else
            nodeIdx(sampleIdx) = nodeID(currNodeidx,[t(1) ~t(1)]);
        end
        
        ib = isBranch(currNodeidx); % If the child of current node is leaf, jump out this while loop
    end
    
    Y(sampleIdx) = pathVal(currNodeidx);
    nodeIdx(sampleIdx) = 0;
    remainCheck = nodeIdx~=0;
    sampleIdx = find(nodeIdx==sum(min(nodeIdx(remainCheck)))); % To avoid empty min output when no non-zero value exists in "nodeIdx"
    remainNum = any(remainCheck);
    
end