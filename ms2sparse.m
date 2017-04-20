%% Predicting mass spectra
function ms = ms2sparse(pepinfo)
%
% Predicting mass spectra for peptides input in variable "pepinfo", which
% is a structure format data containing fields 'mod_info', 'charge' and
% 'pepseq'. 
% The other variable "writetag" is an indicator to tell the function to
% write the predicted mass spectra into ..mgf file... The value 'w'
% indicates writing to the ..mgf file, and 'o' or no value indicates
% outputting them directly.
%

n = numel(pepinfo);
ms = cell(n,1);

pathname = {'yion';'bion';'proly';'prolb';'asply';'asplb';'gluly';'glulb';...
    'yb2y';'yb2b';'pion';'aion';'bw';'yw';'bn';'yn';...
    'pw';'pn';'pwe';'qn';'ypho';'bpho';'ppho';'ymo';...
    'bmo';'pmo';};
numPathname = numel(pathname);

% Generating variables for peptide mass spectra prediction according to the
% modifications...
phosphoPepsidx = zeros(n,1);
numMod = zeros(n,1);
for ii = 1:n
    numMod(ii) = numel(pepinfo(ii).mod_infor);
    for jj = 1:numMod(ii)
        if numel(strfind(lower(pepinfo(ii).mod_infor(jj).name),'phospho'))
            phosphoPepsidx(ii) = 1;
            break;
        end
    end
end

% Divid the peptide information into three parts, no modification,
% phosphorylation and other PTMs, which will be predicted by distinct
% models.
modPepsidx = find(phosphoPepsidx==0&numMod~=0);
phosphoPepsidx = find(phosphoPepsidx==1);
nonmodIdx = find(numMod==0);
npho = numel(phosphoPepsidx);
nmod = numel(modPepsidx);
nnmod = numel(nonmodIdx);

cStr = '123';
uc = [1;2;3];

% Phosphorylated peptides ....
if npho > 0
    % Generating variables...
    nVar = [81;81;81;81;81;81;81;81;81;81;72;81;81;81;81;81;72;...
        72;72;72;81;81;72;81;81;72];
    c = cat(1,pepinfo(phosphoPepsidx).charge);
    
    [theoFraginfo,ns] = vargenerator(pepinfo(phosphoPepsidx),...
        nVar,pathname,numPathname,'p');
    fprintf('Variables of total %d phosphorylated peptides are generated...\n',npho);
    nsnzIdx = find(ns~=0);
    theoFraginfo = rmfield(theoFraginfo,pathname(ns==0));
    numPathname = numel(nsnzIdx);
    
    %%%
    % Predicting peptide intensity
    disp('Predicting fragment ion intensity for phosphorylated peptides...') % Predicting peptides
    for ii = 1:numPathname
        % Predicting intensities for each fragmentation pathway in "pathname"
        allCharge = c(theoFraginfo.(pathname{nsnzIdx(ii)}).info(:,2));
        for cc = 1:3
            cidx = find(allCharge==uc(cc));
            xx = theoFraginfo.(pathname{nsnzIdx(ii)}).v(cidx,:);
            if size(xx,1) > 0
                yypre = ms2predict(xx,pathname{nsnzIdx(ii)},cStr(cc),'p');
                theoFraginfo.(pathname{nsnzIdx(ii)}).info(cidx,3) = yypre;
            end
        end
        fprintf('Intensity prediction for fragmentation pathway of %s is finished.\n',...
            pathname{nsnzIdx(ii)});
    end
    
    % Orgnizing the predicted mass spectra 
    ms2pred = ms2org(theoFraginfo,pathname(nsnzIdx),npho,numPathname);
    ms(phosphoPepsidx) = ms2pred;
end

% Other modified peptides ....
if nmod > 0
    % Generating variables...
    nVar = [59;59;54;59;55;59;55;59;59;59;49;59;59;59;59;59;...
        49;49;49;49;58;58;49;];
    modpathname = pathname([1:20 24:end]);
    
    c = cat(1,pepinfo(modPepsidx).charge);
    
    [theoFraginfo,ns] = vargenerator(pepinfo(modPepsidx),...
        nVar,modpathname,23,'m');
    fprintf('Variables of total %d modified peptides are generated...\n',nmod);
    nsnzIdx = find(ns~=0);
    theoFraginfo = rmfield(theoFraginfo,modpathname(ns==0));
    numPathname = numel(nsnzIdx);
    
    %%%
    % Predicting peptide intensity
    disp('Predicting fragment ion intensity for modified peptides...') % Predicting peptides
    for ii = 1:numPathname
        % Predicting intensities for each fragmentation pathway in "pathname"
        allCharge = c(theoFraginfo.(modpathname{nsnzIdx(ii)}).info(:,2));
        for cc = 1:3
            cidx = find(allCharge==uc(cc));
            xx = theoFraginfo.(modpathname{nsnzIdx(ii)}).v(cidx,:);
            if size(xx,1) > 0
                yypre = ms2predict(xx,modpathname{nsnzIdx(ii)},cStr(cc),'m');
                theoFraginfo.(modpathname{nsnzIdx(ii)}).info(cidx,3) = yypre;
            end
        end
        fprintf('Intensity prediction for fragmentation pathway of %s is finished.\n',...
            modpathname{nsnzIdx(ii)});
    end
    
    % Orgnizing the predicted mass spectra 
    ms2pred = ms2org(theoFraginfo,modpathname(nsnzIdx),nmod,numPathname);
    ms(modPepsidx) = ms2pred;
end

% Non modified peptides ....
if nnmod > 0
    % Generating variables...
    nVar = [49;49;48;48;48;48;48;48;51;51;44;49;49;49;49;49;44;44;44;39;];
    
    c = cat(1,pepinfo(nonmodIdx).charge);
    
    [theoFraginfo,ns] = vargenerator(pepinfo(nonmodIdx),...
        nVar,pathname,20,'n');
    fprintf('Variables of total %d non modified peptides are generated...\n',nnmod);
    nsnzIdx = find(ns~=0);
    theoFraginfo = rmfield(theoFraginfo,pathname(ns==0));
    numPathname = numel(nsnzIdx);
    
    disp('Predicting fragment ion intensity for non modified peptides...') % Predicting peptides
    for ii = 1:numPathname
        % Predicting intensities for each fragmentation pathway in "pathname"
        allCharge = c(theoFraginfo.(pathname{nsnzIdx(ii)}).info(:,2));
        for cc = 1:3
            cidx = find(allCharge==uc(cc));
            xx = theoFraginfo.(pathname{nsnzIdx(ii)}).v(cidx,:);
            if size(xx,1) > 0
                yypre = ms2predict(xx,pathname{nsnzIdx(ii)},cStr(cc),'n');
                theoFraginfo.(pathname{nsnzIdx(ii)}).info(cidx,3) = yypre;
            end
        end
        fprintf('Intensity prediction for fragmentation pathway of %s is finished.\n',...
            pathname{nsnzIdx(ii)});
    end
    
    % Orgnizing the predicted mass spectra
    ms2pred = ms2org(theoFraginfo,pathname(nsnzIdx),nnmod,numPathname);
    ms(nonmodIdx) = ms2pred;
end


%%%%%
%% Variable generation
function [varInfo,nSample] = vargenerator(pepinfo,nVar,pathname,npname,modtype)
%
% Generating variable for mass spectra prediction....

% Initialization
varInfo = cell2struct(cell(size(pathname)),pathname,1);
n = numel(pepinfo);
for ii = 1:npname
    varInfo.(pathname{ii}).v = zeros(n*60,nVar(ii)); % variables
    varInfo.(pathname{ii}).info = zeros(n*60,3); % m/z, indices and predicted intensities
end

disp('Extracting variables and fragment ion information...');
nSample = zeros(size(nVar));

switch modtype
    case 'p'
        
        for ii = 1:n
            res_mass = residumasscal(pepinfo(ii).pepseq,...
                pepinfo(ii).mod_infor,0);
            ionInfo = phosphoions(res_mass,...
                pepinfo(ii).pepseq,...
                pepinfo(ii).charge,...
                pepinfo(ii).mod_infor);
            
            for jj = 1:npname
                if numel(ionInfo.(pathname{jj})) > 0
                    ns = size(ionInfo.(pathname{jj}).vars,1);
                    varInfo.(pathname{jj}).v(nSample(jj)+1:nSample(jj)+ns,:) = ionInfo.(pathname{jj}).vars;
                    varInfo.(pathname{jj}).info(nSample(jj)+1:nSample(jj)+ns,1:2) = ...
                        [ionInfo.(pathname{jj}).mass(:,1) ...
                        repmat(ii,ns,1)];
                    nSample(jj) = nSample(jj)+ns;
                end
            end
            
            if mod(ii,100) == 0
                fprintf('%d peptides have been processed...\n',ii);
            end
        end
        
    otherwise
        
        for ii = 1:n
            res_mass = residumasscal(pepinfo(ii).pepseq,...
                pepinfo(ii).mod_infor,0);
            ionInfo = ions(res_mass,pepinfo(ii).pepseq,...
                pepinfo(ii).charge,pepinfo(ii).mod_infor);
            
            for jj = 1:npname
                if numel(ionInfo.(pathname{jj})) > 0
                    ns = size(ionInfo.(pathname{jj}).vars,1);
                    varInfo.(pathname{jj}).v(nSample(jj)+1:nSample(jj)+ns,:) = ionInfo.(pathname{jj}).vars;
                    varInfo.(pathname{jj}).info(nSample(jj)+1:nSample(jj)+ns,1:2) = ...
                        [ionInfo.(pathname{jj}).mass(:,1) ...
                        repmat(ii,ns,1)];
                    nSample(jj) = nSample(jj)+ns;
                end
            end
            
            if mod(ii,100) == 0
                fprintf('%d peptides have been processed...\n',ii);
            end
        end
        
end

for jj = 1:npname
    varInfo.(pathname{jj}).v = varInfo.(pathname{jj}).v(1:nSample(jj),:);
    varInfo.(pathname{jj}).info = varInfo.(pathname{jj}).info(1:nSample(jj),:);
end

%%%
%% Mass spectra re-orgnization 
function ms2pred = ms2org(varInfo,pathname,n,npname)
%
% Re-orgnizing and processing the predicted intensity for each peptide
% sequence.

ms2pred = cell(n,1);

for ii = 1:n
    predms = zeros(1500,2);
    nMS = 0;
    for jj = 1:npname
        msIdx = find(varInfo.(pathname{jj}).info(:,2)==ii);
        predms(nMS+1:nMS+numel(msIdx),:) = varInfo.(pathname{jj}).info(msIdx,[1 3]);
        nMS = nMS+numel(msIdx);
    end
    
    % %% Processing the predicted intensity
    predms(predms(:,2)<=0,:) = [];
    
    % cluster peaks with m/z within 0.1
    tempMS = predms;
    [~,sortIdx] = sort(tempMS(:,2),'descend');
    tempMS = tempMS(sortIdx,:);
    numMS = size(tempMS,1);
    numMSidx = 0;
    delidx = false(numMS,1);
    while numMSidx < numMS && ...
            ~all(delidx(numMSidx+1:end))
        
        numMSidx = numMSidx+find(delidx(numMSidx+1:end)==0,1,'first');
        selidx = find(tempMS(:,1)>=tempMS(numMSidx,1)-0.1&...
            tempMS(:,1)<=tempMS(numMSidx,1)+0.1);
        delidx(selidx(selidx~=numMSidx)) = true;
        
    end
    predms = predms(sortIdx(~delidx),:);
    predms(:,2) = predms(:,2)/max(predms(:,2));
    
    % Sort predicted mass spectra according to the ascending order of m/z.
    [~,sortIdx] = sort(predms(:,1));
    ms2pred{ii,1} = predms(sortIdx,:);
end

disp('Prediction of all peptides is finished...');