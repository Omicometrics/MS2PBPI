%% Predict y from stochastic gradient boosting tree (SGBTree)
function Y = sgbtreepredict(X,params,evalStr,nodeID,isBranch,pathVal)
%
% Predict the intensity of mass spectra according to the models constructed
% by stochastic gradient boosting tree (SGBTee)
%
% Nai-ping Dong, PolyU in HK
% Email: np.dong572@gmail.com
% 1/2/2014

% if nargin < 2
%     error('The Variable Matrix Must Be Input with Models!');
% elseif nargin > 2
%     error('Too Many Inputs!');
% end

X(:,params.delidx) = [];
m = size(X,1);

Y = zeros(m,1);
Y = Y+0;

for jj = 1:100
%     [evalStr,nodeID,isBranch,pathVal] = matrix2tree(models.fun(jj).m);
    
    Yp = zeros(m,1);
    nodeIdx = ones(m,1);
    sampleIdx = (1:m)';
    remainNum = true;
    
    while remainNum
        
        ib = isBranch(jj).v(nodeIdx(sampleIdx(1)));
        currNodeidx = nodeIdx(sampleIdx(1));
        
        while ib
            
            currNodeidx = nodeIdx(sampleIdx(1));
            t = eval(evalStr(jj).s{currNodeidx});
            
            if any(t)
                nodeIdx(sampleIdx(t)) = nodeID(jj).v(currNodeidx,1);
                nodeIdx(sampleIdx(~t))= nodeID(jj).v(currNodeidx,2);
                sampleIdx = sampleIdx(t);
            else
                nodeIdx(sampleIdx) = nodeID(jj).v(currNodeidx,[t(1) ~t(1)]);
            end
            
            ib = isBranch(jj).v(currNodeidx); % If the child of current node is leaf, jump out this while loop
        end
        
        Yp(sampleIdx) = pathVal(jj).v(currNodeidx);
        nodeIdx(sampleIdx) = 0;
        remainCheck = nodeIdx~=0;
        sampleIdx = find(nodeIdx==sum(min(nodeIdx(remainCheck)))); % To avoid empty min output when no non-zero value exists in "nodeIdx"
        remainNum = any(remainCheck); 
        
    end
    
    Y = Y+params.shr*params.m(jj)*Yp;
    
end