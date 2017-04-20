%% Ions and variables calculation
function ionInfo = phosphoions(res_mass,pepseq,charge,mod_infor)
%
% This function is to calculate m/z of ions and generate variables for
% phosphorylated peptides.
%

pathname = {'yion';'bion';'proly';'prolb';'asply';'asplb';'gluly';'glulb';...
    'yb2y';'yb2b';'pion';'aion';'bw';'yw';'bn';'yn';...
    'pw';'pn';'pwe';'qn';'ypho';'bpho';'ppho';'ymo';...
    'bmo';'pmo';};

ionInfo = cell2struct(cell(size(pathname)),pathname,1);

numPath = numel(pathname);
pm_avg = sum(residumasscal(pepseq,mod_infor,1))+18.015;
pm = sum(res_mass)+18.015;

comass = 28.0102;
watermass = 18.0153;
hydrogenmass = 1.0079;
ammoniamass = 17.0306;
resnum = numel(pepseq);

% AA identity.
C = 'ARNDCQEGHILKMFPSTWYV';
identAA = zeros(numel(pepseq),1);
for ii = 1:resnum
    identAA(ii) = strfind(C,pepseq(ii));
end

proIdx = strfind(pepseq,'P'); proIdx(proIdx==1) = [];
aspIdx = strfind(pepseq,'D'); aspIdx(aspIdx==resnum) = [];
gluIdx = strfind(pepseq,'E'); gluIdx(gluIdx==resnum) = [];


%%%%%%%%% Ion mass calculating and feature assigning...
% All series...
bseries = zeros(resnum-1,2);
yseries = bseries;
bseries_name = cell(resnum-1,1);
yseries_name = bseries_name;
aseries_name = bseries_name;
bwloss_name = bseries_name;
ywloss_name = bseries_name;
bnloss_name = bseries_name;
ynloss_name = bseries_name;
for ii = 1:resnum-1
    bseries(ii,1) = sum(res_mass(1:ii))+hydrogenmass;
    bseries(ii,2) = ii;
    
    yseries(ii,1) = sum(res_mass(ii+1:end))+hydrogenmass+watermass;
    yseries(ii,2) = resnum-ii;
    
    bseries_name{ii} = ['b' num2str(ii)];
    yseries_name{ii} = ['y' num2str(resnum-ii)];
    aseries_name{ii} = ['a' num2str(ii)];
    bwloss_name{ii} = ['b' num2str(ii) '-H2O'];
    ywloss_name{ii} = ['y' num2str(resnum-ii) '-H2O'];
    bnloss_name{ii} = ['b' num2str(ii) '-NH3'];
    ynloss_name{ii} = ['y' num2str(resnum-ii) '-NH3'];
end

bVars = ionfeatures(1:resnum-1,charge,'b',identAA,1,bseries(:,1),pm_avg,mod_infor);
yVars = ionfeatures(2:resnum,charge,'y',identAA,1,yseries(:,1),pm_avg,mod_infor);

if numel(proIdx) > 0 % proline effect
    proIonb = bseries(proIdx-1,:);
    ionInfo.('prolb').name = 'prolb';
    ionInfo.('prolb').ionname = bseries_name(proIdx-1);
    ionInfo.('prolb').mass = proIonb;
    ionInfo.('prolb').vars = bVars(proIdx-1,:);
    
    proIony = yseries(proIdx-1,:);
    ionInfo.('proly').name = 'proly';
    ionInfo.('proly').ionname = yseries_name(proIdx-1);
    ionInfo.('proly').mass = proIony;
    ionInfo.('proly').vars = yVars(proIdx-1,:);
end

if numel(aspIdx) > 0 % aspartic acid effect
    aspIonb = bseries(aspIdx,:);
    ionInfo.('asplb').name = 'asplb';
    ionInfo.('asplb').ionname = bseries_name(aspIdx);
    ionInfo.('asplb').mass = aspIonb;
    ionInfo.('asplb').vars = bVars(aspIdx,:);
    
    aspIony = yseries(aspIdx,:);
    ionInfo.('asply').name = 'asply';
    ionInfo.('asply').ionname = yseries_name(aspIdx);
    ionInfo.('asply').mass = aspIony;
    ionInfo.('asply').vars = yVars(aspIdx,:);
end

if numel(gluIdx) > 0 % glutamic acid effect
    gluIonb = bseries(gluIdx,:);
    ionInfo.('glulb').name = 'glulb';
    ionInfo.('glulb').ionname = bseries_name(gluIdx);
    ionInfo.('glulb').mass = gluIonb;
    ionInfo.('glulb').vars = bVars(gluIdx,:);
    
    gluIony = yseries(gluIdx,:);
    ionInfo.('gluly').name = 'gluly';
    ionInfo.('gluly').ionname = yseries_name(gluIdx);
    ionInfo.('gluly').mass = gluIony;
    ionInfo.('gluly').vars = yVars(gluIdx,:);
end

%%% b2-yn-2
% y
ionInfo.('yb2y').name = 'yb2y';
ionInfo.('yb2y').ionname = yseries_name(2);
ionInfo.('yb2y').mass = yseries(2,:);
ionInfo.('yb2y').vars = yVars(2,:);

ionInfo.('yb2b').name = 'yb2b';
ionInfo.('yb2b').ionname = bseries_name(2);
ionInfo.('yb2b').mass = bseries(2,:);
ionInfo.('yb2b').vars = bVars(2,:);

% bx-yz
delidx = [proIdx-1 aspIdx gluIdx 2];
bidx = (1:resnum-1)'; bidx(delidx) = [];
yidx = (1:resnum-1)'; yidx(delidx) = [];
if size(bidx,1) > 0
    ionInfo.('bion').name = 'bion';
    ionInfo.('bion').ionname = bseries_name(bidx);
    ionInfo.('bion').mass = bseries(bidx,:);
    ionInfo.('bion').vars = bVars(bidx,:);
    
    ionInfo.('yion').name = 'yion';
    ionInfo.('yion').ionname = yseries_name(yidx);
    ionInfo.('yion').mass = yseries(yidx,:);
    ionInfo.('yion').vars = yVars(yidx,:);
end

%%% Precursor ions
pionmass = [pm_avg+hydrogenmass resnum];
ionInfo.('pion').name = 'pion';
ionInfo.('pion').ionname = {'p'};
ionInfo.('pion').mass = pionmass;
tempVar = ionfeatures(0,1,'p',identAA,1,pionmass(1),pm,mod_infor);
ionInfo.('pion').vars = tempVar;

%%% a ions and neutral loss of y and b ions
bwloss = [bseries(:,1)-watermass bseries(:,2)];
ywloss = [yseries(:,1)-watermass yseries(:,2)];
pwloss = [pionmass(1)-watermass resnum];
bnloss = [bseries(:,1)-ammoniamass bseries(:,2)]; 
ynloss = [yseries(:,1)-ammoniamass yseries(:,2)];
pnloss = [pionmass(1)-ammoniamass resnum];

% a ions
ionInfo.('aion').name = 'aion';
ionInfo.('aion').ionname = aseries_name;
atempVar = bVars; % The variables can be derived from b ions directly
atempVar(:,12) = atempVar(:,12)-comass;
atempVar(:,13) = atempVar(:,12)./((atempVar(:,11)+(charge-1)*hydrogenmass)/charge);
atempVar(:,14) = 1-(atempVar(:,12)-hydrogenmass)./(atempVar(:,11)-hydrogenmass);
ionInfo.('aion').mass = [bseries(:,1)-comass bseries(:,2)];
ionInfo.('aion').vars = atempVar;

% water loss of y and b ions
% // This section is not correct, now is corrected
ionInfo.('bw').name = 'bw';
ionInfo.('bw').ionname = bwloss_name;
ionInfo.('bw').mass = bwloss;
bwtempVar = bVars; % The variables can be derived from b ions directly
bwtempVar(:,12) = bwloss(:,1);
bwtempVar(:,13) = bwtempVar(:,12)./((bwtempVar(:,11)+(charge-1)*hydrogenmass)/charge);
bwtempVar(:,14) = 1-(bwloss(:,1)-hydrogenmass)./(bwtempVar(:,11)-hydrogenmass);
ionInfo.('bw').vars = bwtempVar;

ionInfo.('yw').name = 'yw';
ionInfo.('yw').ionname = ywloss_name;
ionInfo.('yw').mass = ywloss;
ywtempVar = yVars; % The variables can be derived from y ions directly
ywtempVar(:,12) = ywloss(:,1);
ywtempVar(:,13) = ywtempVar(:,12)./((ywtempVar(:,11)+(charge-1)*hydrogenmass)/charge);
ywtempVar(:,14) = 1-(ywloss(:,1)--hydrogenmass)./(ywtempVar(:,11)-hydrogenmass);
ionInfo.('yw').vars = ywtempVar;

%%% water loss of precursor ions
ionInfo.('pw').name = 'pw';
ionInfo.('pw').ionname = {'p-H2O'};
ionInfo.('pw').mass = pwloss;
ionInfo.('pw').vars = ionInfo.('pion').vars;

%%% ammonia loss of y, b and precursor ions
% // This section is not correct, now is corrected
ionInfo.('bn').name = 'bn';
ionInfo.('bn').ionname = bnloss_name;
ionInfo.('bn').mass = bnloss;
bntempVar = bVars; % The variables can be derived from b ions directly
bntempVar(:,12) = bnloss(:,1);
bntempVar(:,13) = bntempVar(:,12)./((bntempVar(:,11)+(charge-1)*hydrogenmass)/charge);
bntempVar(:,14) = 1-(bnloss(:,1)--hydrogenmass)./(bntempVar(:,11)-hydrogenmass);
ionInfo.('bn').vars = bntempVar;

ionInfo.('yn').name = 'yn';
ionInfo.('yn').ionname = ynloss_name;
ionInfo.('yn').mass = ynloss;
yntempVar = yVars; % The variables can be derived from y ions directly
yntempVar(:,12) = ynloss(:,1);
yntempVar(:,13) = yntempVar(:,12)./((yntempVar(:,11)+(charge-1)*hydrogenmass)/charge);
yntempVar(:,14) = 1-(ynloss(:,1)--hydrogenmass)./(yntempVar(:,11)-hydrogenmass);
ionInfo.('yn').vars = yntempVar;

ionInfo.('pn').name = 'pn';
ionInfo.('pn').vars = ionInfo.('pion').vars;
ionInfo.('pn').ionname = {'p-NH3'};
ionInfo.('pn').mass = pnloss;

%%%% Modification infor
% Get the information of modified residues.
metoxyidx = zeros(size(mod_infor));
stphosidx = zeros(size(mod_infor));
% //
if numel(mod_infor) > 0
    modaapos = cat(1,mod_infor.pos);
    modname = cat(1,{mod_infor.name})';
    for ii = 1:numel(modname)
        currmodname = modname{ii};
        if numel(strfind(lower(currmodname),'oxidation')) > 0
            metoxyidx(ii) = modaapos(ii);
        elseif numel(strfind(lower(currmodname),'phospho')) > 0
            stphosidx(ii) = modaapos(ii);
        end
    end
end
% //
metoxyidx(metoxyidx==0) = [];
stphosidx(stphosidx==0) = [];

% Because the deamidation of Asn is derived from the PNGase F-catalyzed
% deglycosylation of Asn in glycopeptides, the residue is conversed into
% aspartic acid.
% nonasnseq = pepseq;
% if numel(asndeamidx) > 0
%     nonasnseq(asndeamidx) = 'D';
% end

%%%% ====================================================================
% Add the information of modification of residues into the ion series.
% oxidation of Metheonine
if numel(metoxyidx) > 0
    
    bmetoxyidx = metoxyidx(metoxyidx~=resnum);
    ymetoxyidx = metoxyidx(metoxyidx~=1);
    if numel(bmetoxyidx) > 0
        bsloss = [bseries(min(bmetoxyidx):end,1)-63.9982 bseries(min(bmetoxyidx):end,2)];
        bsloss_name = cell(resnum-min(bmetoxyidx),1);
        for ii = min(bmetoxyidx):resnum-1
            bsloss_name{ii-min(bmetoxyidx)+1} = [bseries_name{ii} '-CH3SOH'];
        end
        
        ionInfo.('bmo').name = 'bmo';
        ionInfo.('bmo').ionname = bsloss_name;
        ionInfo.('bmo').mass = bsloss;
        bstempVar = bVars(min(bmetoxyidx):end,:); % The variables can be derived from b ions directly
        bstempVar(:,12) = bsloss(:,1);
        bstempVar(:,13) = bstempVar(:,12)./((bstempVar(:,11)+(charge-1)*hydrogenmass)/charge);
        bstempVar(:,14) = 1-(bsloss(:,1)-hydrogenmass)./(bstempVar(:,11)-hydrogenmass);
        ionInfo.('bmo').vars = bstempVar;
    end
    
    if numel(ymetoxyidx) > 0
        ysloss = [yseries(1:max(ymetoxyidx)-1,1)-63.9982 yseries(1:max(ymetoxyidx)-1,2)];
        ysloss_name = cell(max(ymetoxyidx)-1,1);
        for ii = 1:max(resnum-ymetoxyidx+1)
            ysloss_name{ii,1} = [yseries_name{ii} '-CH3SOH'];
        end
        
        ionInfo.('ymo').name = 'ymo';
        ionInfo.('ymo').ionname = ysloss_name;
        ionInfo.('ymo').mass = ysloss;
        ystempVar = yVars(1:max(ymetoxyidx)-1,:); % The variables can be derived from b ions directly
        ystempVar(:,12) = ysloss(:,1);
        ystempVar(:,13) = ystempVar(:,12)./((ystempVar(:,11)+(charge-1)*hydrogenmass)/charge);
        ystempVar(:,14) = 1-(ysloss(:,1)-hydrogenmass)./(ystempVar(:,11)-hydrogenmass);
        ionInfo.('ymo').vars = ystempVar;
    end
    
%     psloss = [pionmass(1)-63.9982 resnum];
%     ionInfo.('pmo').name = 'pmo';
%     ionInfo.('pmo').ionname = {'p-CH3SOH'};
%     ionInfo.('pmo').mass = psloss;
%     ionInfo.('pmo').vars = ionInfo.('pion').vars;
end

% phosphorylations of Serine and Threonine
if numel(stphosidx) > 0
    
    bstphosidx = stphosidx(stphosidx~=resnum);
    ystphosidx = stphosidx(stphosidx~=1);
    if numel(bstphosidx) > 0
        bploss = [bseries(min(bstphosidx):end,1)-97.9768 bseries(min(bstphosidx):end,2)];
        bploss_name = cell(resnum-min(bstphosidx),1);
        for ii = min(bstphosidx):resnum-1
            bploss_name{ii-min(bstphosidx)+1} = [bseries_name{ii} '-H3PO4'];
        end
        
        ionInfo.('bpho').name = 'bpho';
        ionInfo.('bpho').ionname = bploss_name;
        ionInfo.('bpho').mass = bploss;
        bptempVar = bVars(min(bstphosidx):end,:); % The variables can be derived from b ions directly
        bptempVar(:,12) = bploss(:,1);
        bptempVar(:,13) = bptempVar(:,12)./((bptempVar(:,11)+(charge-1)*hydrogenmass)/charge);
        bptempVar(:,14) = 1-(bploss(:,1)-hydrogenmass)./(bptempVar(:,11)-hydrogenmass);
        ionInfo.('bpho').vars = bptempVar;
    end
    
    if numel(ystphosidx) > 0
        yploss = [yseries(1:max(ystphosidx)-1,1)-97.9768 yseries(1:max(ystphosidx)-1,2)];
        yploss_name = cell(max(ystphosidx)-1,1);
        for ii = 1:max(resnum-ystphosidx+1)
            yploss_name{ii,1} = [yseries_name{ii} '-H3PO4'];
        end
        
        ionInfo.('ypho').name = 'ypho';
        ionInfo.('ypho').ionname = yploss_name;
        ionInfo.('ypho').mass = yploss;
        yptempVar = yVars(1:max(ystphosidx)-1,:); % The variables can be derived from b ions directly
        yptempVar(:,12) = yploss(:,1);
        yptempVar(:,13) = yptempVar(:,12)./((yptempVar(:,11)+(charge-1)*hydrogenmass)/charge);
        yptempVar(:,14) = 1-(yploss(:,1)-hydrogenmass)./(yptempVar(:,11)-hydrogenmass);
        ionInfo.('ypho').vars = yptempVar;
    end
    
    pploss = [pionmass(1)-97.9768 resnum];
    ionInfo.('ppho').name = 'ppho';
    ionInfo.('ppho').ionname = {'p-H3PO4'};
    ionInfo.('ppho').mass = pploss;
    ionInfo.('ppho').vars = ionInfo.('pion').vars;
end

%%%% ====================================================================
% doubly and triply charged, with doubly and triply charged fragments.
if charge >= 2;
    for ii = 1:numPath
        if numel(ionInfo.(pathname{ii})) > 0
            c2idx = find(ionInfo.(pathname{ii}).mass(:,2) >= 3);
            c3idx = find(ionInfo.(pathname{ii}).mass(:,2) >= 5);
            if numel(c2idx) > 0
                c2Ionmass = zeros(numel(c2idx),2); c2Ionname = cell(size(c2idx));
                for jj = 1:numel(c2idx)
                    c2Ionmass(jj,1) = (ionInfo.(pathname{ii}).mass(c2idx(jj),1)+hydrogenmass)/2;
                    c2Ionmass(jj,2) = ionInfo.(pathname{ii}).mass(c2idx(jj),2);
                    c2Ionname{jj} = [ionInfo.(pathname{ii}).ionname{c2idx(jj)} '^2'];
                end
                
                tempVar = ionInfo.(pathname{ii}).vars(c2idx,:);
                tempVar(:,50) = 2;
                
                if ~strcmp(ionInfo.(pathname{ii}).ionname{1}(1),'p') && ~strcmp(ionInfo.(pathname{ii}).ionname{1}(1),'q')
                    tempVar(:,12) = (tempVar(:,12)+hydrogenmass)/2;
                    tempVar(:,13) = tempVar(:,12)./((tempVar(:,11)+(charge-1)*hydrogenmass)/charge);
                end
                
                ionInfo.(pathname{ii}).vars = [ionInfo.(pathname{ii}).vars;tempVar];
                ionInfo.(pathname{ii}).mass = [ionInfo.(pathname{ii}).mass;c2Ionmass];
                ionInfo.(pathname{ii}).ionname = [ionInfo.(pathname{ii}).ionname;c2Ionname];
            end
            
            if charge == 3
                if numel(c3idx) > 0
                    c3Ionmass = zeros(numel(c3idx),2); c3Ionname = cell(size(c3idx));
                    for jj = 1:numel(c3idx)
                        c3Ionmass(jj,1) = (ionInfo.(pathname{ii}).mass(c3idx(jj),1)+2*hydrogenmass)/3;
                        c3Ionmass(jj,2) = ionInfo.(pathname{ii}).mass(c3idx(jj),2);
                        c3Ionname{jj} = [ionInfo.(pathname{ii}).ionname{c3idx(jj)} '^3'];
                    end
                    
                    tempVar = ionInfo.(pathname{ii}).vars(c3idx,:);
                    tempVar(:,50) = 3;
                    
                    if ~strcmp(ionInfo.(pathname{ii}).ionname{1}(1),'p') && ~strcmp(ionInfo.(pathname{ii}).ionname{1}(1),'q')
                        tempVar(:,12) = (tempVar(:,12)+2*hydrogenmass)/3;
                        tempVar(:,13) = tempVar(:,12)./((tempVar(:,11)+(charge-1)*hydrogenmass)/charge);
                    end
                    
                    ionInfo.(pathname{ii}).vars = [ionInfo.(pathname{ii}).vars;tempVar];
                    ionInfo.(pathname{ii}).mass = [ionInfo.(pathname{ii}).mass;c3Ionmass];
                    ionInfo.(pathname{ii}).ionname = [ionInfo.(pathname{ii}).ionname;c3Ionname];
                end
            end
        end
    end
end

%%%%%%%%%%%
%% Feature extraction for the prediction of tandem mass spectra.
function X = ionfeatures(frag_pos,charge,iontype,identAA,ioncharge,ionmass,precursormass,mod_infor)
%
% Features: Peptide length; identity of residue positions N-terminal to
% fragmentation site; identity of residue positions C-terminal to
% fragmentation site; fragmentation position (distance to C-termini and
% N-termini); proton mobility, charge state, difference between precursor
% m/z and fragment m/z, position of residue 'Arg' in the fragment ion to
% N-terminal and C-terminal side, if no, 0, if multiple, retain the
% minimum.

load AAdescriptor;

if ~strcmp(iontype,'p')
    X = zeros(numel(frag_pos),81);
    for ii = 1:numel(frag_pos)
        X(ii,1) = identAA(1);
        X(ii,2) = identAA(end);
        X(ii,3) = numel(identAA);
        switch iontype
            case 'y'
                X(ii,4) = identAA(frag_pos(ii)-1);
                X(ii,5) = identAA(frag_pos(ii));
                X(ii,6) = numel(find(identAA(frag_pos(ii):end)==15));
                X(ii,7) = numel(find(identAA(frag_pos(ii):end)==4));
                X(ii,8) = numel(find(identAA(frag_pos(ii):end)==7));
                X(ii,9) = numel(find(identAA(frag_pos(ii):end)==9))+...
                    numel(find(identAA(frag_pos(ii):end)==12))+...
                    numel(find(identAA(frag_pos(ii):end)==2));
                
                if frag_pos(ii) <= X(ii,3)-1 % C-terminal
                    X(ii,19) = identAA(frag_pos(ii)+1);
                end
                if frag_pos(ii) >= 3 % N-terminal
                    X(ii,20) = identAA(frag_pos(ii)-2);
                end
                X(ii,21) = mean(AAdesc.d(identAA(frag_pos(ii):end),3));
                
                if frag_pos(ii) >= 4
                    X(ii,51) = identAA(frag_pos(ii)-3);
                end
                if frag_pos(ii) <= X(ii,3)-2
                    X(ii,52) = identAA(frag_pos(ii)+2);
                end
                X(ii,53) = (sum(AAdesc.numatom(identAA(frag_pos(ii):end)))+3)*3-6;
                
            case {'b';'a'}
                X(ii,4) = identAA(frag_pos(ii)); % N-terminal
                X(ii,5) = identAA(frag_pos(ii)+1); % C-terminal
                X(ii,6) = numel(find(identAA(1:frag_pos(ii))==15));
                X(ii,7) = numel(find(identAA(1:frag_pos(ii))==4));
                X(ii,8) = numel(find(identAA(1:frag_pos(ii))==7));
                X(ii,9) = numel(find(identAA(1:frag_pos(ii))==9))+...
                    numel(find(identAA(1:frag_pos(ii))==12))+...
                    numel(find(identAA(1:frag_pos(ii))==2));
                
                % ...
                if frag_pos(ii) <= X(ii,3)-2 % C-terminal
                    X(ii,19) = identAA(frag_pos(ii)+2);
                end
                if frag_pos(ii) >= 2 % N-terminal
                    X(ii,20) = identAA(frag_pos(ii)-1);
                end
                X(ii,21) = mean(AAdesc.d(identAA(1:frag_pos(ii)),3));
                
                % ...
                if frag_pos(ii) >= 3
                    X(ii,51) = identAA(frag_pos(ii)-2);
                end
                if frag_pos(ii) <= X(ii,3)-3
                    X(ii,52) = identAA(frag_pos(ii)+3);
                end
                X(ii,53) = (sum(AAdesc.numatom(identAA(1:frag_pos(ii))))+3)*3-6;
        end
                
        if numel(find(identAA==2)) >= charge
            X(ii,10) = 0;
        elseif numel(find(identAA==9))+numel(find(identAA==12)) >= charge
            X(ii,10) = 1;
        else
            X(ii,10) = 2;
        end
        
        X(ii,11) = precursormass+1.0079;
        X(ii,12) = (ionmass(ii)+(ioncharge-1)*1.0079)/ioncharge;
        X(ii,13) = X(ii,12)/((X(ii,11)+(charge-1)*1.0079)/charge);
        X(ii,14) = 1-(ionmass(ii)-1.0079)/(X(ii,11)-1.0079);
        X(ii,15) = frag_pos(ii);
        X(ii,16) = X(ii,3)-frag_pos(ii);
        X(ii,17) = frag_pos(ii)/X(ii,3);
        X(ii,18) = 1-X(ii,17);
        
        X(ii,22) = AAdesc.d(X(ii,4),1);
        X(ii,23) = AAdesc.d(X(ii,4),2);
        X(ii,24) = AAdesc.d(X(ii,4),3);
        X(ii,25) = AAdesc.d(X(ii,4),4);
        X(ii,26) = AAdesc.d(X(ii,5),1);
        X(ii,27) = AAdesc.d(X(ii,5),2);
        X(ii,28) = AAdesc.d(X(ii,5),3);
        X(ii,29) = AAdesc.d(X(ii,5),4);
        
        aanum = zeros(1,20);
        for jj = 1:20
            aanum(jj) = numel(find(identAA==jj));
        end
        X(ii,30:49) = aanum;
        
        X(ii,50) = ioncharge;
        for jj = 1:numel(mod_infor)
            if numel(strfind(mod_infor(jj).name,'Carbamidomethyl')) > 0
                X(ii,54) = X(ii,54)+1;
            end
            if numel(strfind(mod_infor(jj).name,'Acetyl')) > 0
                X(ii,55) = X(ii,55)+1;
            end
            if numel(strfind(mod_infor(jj).name,'Deamidat')) > 0
                X(ii,56) = X(ii,56)+1;
            end
            if numel(strfind(lower(mod_infor(jj).name),'pyro-glu')) > 0
                X(ii,57) = X(ii,57)+1;
            end
            if numel(strfind(mod_infor(jj).name,'Oxidation')) > 0
                X(ii,58) = X(ii,58)+1;
            end
        end
        
        modpos = cat(1,mod_infor.pos);
        X(ii,59) = numel(find(modpos>=frag_pos(ii)-2&modpos>=frag_pos(ii)+1));
        
        phosphoPos = zeros(size(modpos));
        for jj = 1:numel(modpos)
            if identAA(modpos(jj)) == 16
                X(ii,60) = X(ii,60)+1;
                phosphoPos(jj) = modpos(jj);
            elseif identAA(modpos(jj)) == 19
                X(ii,61) = X(ii,61)+1;
                phosphoPos(jj) = modpos(jj);
            elseif identAA(modpos(jj)) == 17
                X(ii,62) = X(ii,62)+1;
                phosphoPos(jj) = modpos(jj);
            end
        end
        X(ii,63) = sum(X(ii,60:62));
        
        % The special amino acids bonded to the modified residues
        % The special residues considered are: Pro, Arg, Lys, His, Tyr,
        % Ser, Thr, oxidated Met.
        phosphoPos(phosphoPos==0) = [];
        for jj = 1:numel(phosphoPos)
            if phosphoPos(jj) ~= X(ii,3)
                if identAA(phosphoPos(jj)+1) == 15
                    X(ii,64) = X(ii,64)+1;
                elseif identAA(phosphoPos(jj)+1) == 2
                    X(ii,65) = X(ii,65)+1;
                elseif identAA(phosphoPos(jj)+1) == 12
                    X(ii,66) = X(ii,66)+1;
                elseif identAA(phosphoPos(jj)+1) == 9
                    X(ii,67) = X(ii,67)+1;
                elseif identAA(phosphoPos(jj)+1) == 19
                    X(ii,68) = X(ii,68)+1;
                elseif identAA(phosphoPos(jj)+1) == 16
                    X(ii,69) = X(ii,69)+1;
                elseif identAA(phosphoPos(jj)+1) == 17
                    X(ii,70) = X(ii,70)+1;
                elseif identAA(phosphoPos(jj)+1) == 13 && ...
                        numel(find(modpos==phosphoPos(jj)+1)) > 0
                    if numel(strfind(mod_infor(modpos==phosphoPos(jj)+1).name,...
                            'Oxidation')) > 0
                        X(ii,71) = X(ii,71)+1;
                    end
                end
            end
            if phosphoPos(jj) ~= 1
                if identAA(phosphoPos(jj)-1) == 15
                    X(ii,72) = X(ii,72)+1;
                elseif identAA(phosphoPos(jj)-1) == 2
                    X(ii,73) = X(ii,73)+1;
                elseif identAA(phosphoPos(jj)-1) == 12
                    X(ii,74) = X(ii,74)+1;
                elseif identAA(phosphoPos(jj)-1) == 9
                    X(ii,75) = X(ii,75)+1;
                elseif identAA(phosphoPos(jj)-1) == 19
                    X(ii,76) = X(ii,76)+1;
                elseif identAA(phosphoPos(jj)-1) == 16
                    X(ii,77) = X(ii,77)+1;
                elseif identAA(phosphoPos(jj)-1) == 17
                    X(ii,78) = X(ii,78)+1;
                elseif identAA(phosphoPos(jj)-1) == 13 && ...
                        numel(find(modpos==phosphoPos(jj)-1)) > 0
                    if numel(strfind(mod_infor(modpos==phosphoPos(jj)-1).name,...
                            'Oxidation')) > 0
                        X(ii,79) = X(ii,79)+1;
                    end
                end
            end
        end
        X(ii,80) = min(phosphoPos);
        X(ii,81) = max(phosphoPos);
    end
else
    
    X = zeros(1,72);
    
    X(1) = identAA(1);
    X(2) = identAA(end);
    X(3) = numel(identAA);
    X(4) = numel(find(identAA==15));
    X(5) = numel(find(identAA==4));
    X(6) = numel(find(identAA==7));
    X(7) = numel(find(identAA==9))+...
        numel(find(identAA==12))+...
        numel(find(identAA==2));
    if numel(find(identAA==2)) >= charge
        X(8) = 0;
    elseif numel(find(identAA==9))+numel(find(identAA==12)) >= charge
        X(8) = 1;
    else
        X(8) = 2;
    end
    
    for jj = 1:numel(identAA)-1
        if (identAA(jj) == 4 && identAA(jj+1) == 15) || ...
                (identAA(jj) == 7 && identAA(jj+1) == 15)
            X(9) = 1;
            break;
        end
    end
    
    if identAA(3) == 15; X(10) = 1; end
    if identAA(3) == 4; X(11) = 1; end
    
    X(12) = precursormass+1.0079;
    X(13) = identAA(2);
    X(14) = identAA(end-1);
    X(15) = AAdesc.d(X(1),1);
    X(16) = AAdesc.d(X(1),2);
    X(17) = AAdesc.d(X(1),3);
    X(18) = AAdesc.d(X(1),4);
    X(19) = mean(AAdesc.d(identAA,3));
    
    aanum = zeros(1,20);
    for jj = 1:20
        aanum(jj) = numel(find(identAA==jj));
    end
    X(20:39) = aanum;
    X(40) = identAA(end-2);
    
    X(41) = identAA(3);
    X(42) = identAA(end-2);
    X(43) = (sum(AAdesc.numatom(identAA))+3)*3-6;
    
    for jj = 1:numel(mod_infor)
        if numel(strfind(mod_infor(jj).name,'Carbamidomethyl')) > 0
            X(44) = X(44)+1;
        end
        if numel(strfind(mod_infor(jj).name,'Acetyl')) > 0
            X(45) = X(45)+1;
        end
        if numel(strfind(mod_infor(jj).name,'Deamidat')) > 0
            X(46) = X(46)+1;
        end
        if numel(strfind(lower(mod_infor(jj).name),'pyro-glu')) > 0
            X(47) = X(47)+1;
        end
        if numel(strfind(mod_infor(jj).name,'Oxidation')) > 0
            X(48) = X(48)+1;
        end
    end
    
    modpos = cat(1,mod_infor.pos);
    X(49) = numel(find(modpos==1));
    X(50) = ioncharge;
    
    phosphoPos = zeros(size(modpos));
    for jj = 1:numel(modpos)
        if identAA(modpos(jj)) == 16
            X(51) = X(51)+1;
            phosphoPos(jj) = modpos(jj);
        elseif identAA(modpos(jj)) == 19
            X(52) = X(52)+1;
            phosphoPos(jj) = modpos(jj);
        elseif identAA(modpos(jj)) == 17
            X(53) = X(53)+1;
            phosphoPos(jj) = modpos(jj);
        end
    end
    X(54) = sum(X(51:53));
    
    % The special amino acids bonded to the modified residues
    % The special residues considered are: Pro, Arg, Lys, His, Tyr,
    % Ser, Thr, oxidated Met.
    phosphoPos(phosphoPos==0) = [];
    for jj = 1:numel(phosphoPos)
        if phosphoPos(jj) ~= X(3)
            if identAA(phosphoPos(jj)+1) == 15
                X(55) = X(55)+1;
            elseif identAA(phosphoPos(jj)+1) == 2
                X(56) = X(56)+1;
            elseif identAA(phosphoPos(jj)+1) == 12
                X(57) = X(57)+1;
            elseif identAA(phosphoPos(jj)+1) == 9
                X(58) = X(58)+1;
            elseif identAA(phosphoPos(jj)+1) == 19
                X(59) = X(59)+1;
            elseif identAA(phosphoPos(jj)+1) == 16
                X(60) = X(60)+1;
            elseif identAA(phosphoPos(jj)+1) == 17
                X(61) = X(61)+1;
            elseif identAA(phosphoPos(jj)+1) == 13 && ...
                    numel(find(modpos==phosphoPos(jj)+1)) > 0
                if numel(strfind(mod_infor(modpos==phosphoPos(jj)+1).name,...
                        'Oxidation')) > 0
                    X(62) = X(62)+1;
                end
            end
        end
        if phosphoPos(jj) ~= 1
            if identAA(phosphoPos(jj)-1) == 15
                X(63) = X(63)+1;
            elseif identAA(phosphoPos(jj)-1) == 2
                X(64) = X(64)+1;
            elseif identAA(phosphoPos(jj)-1) == 12
                X(65) = X(65)+1;
            elseif identAA(phosphoPos(jj)-1) == 9
                X(66) = X(66)+1;
            elseif identAA(phosphoPos(jj)-1) == 19
                X(67) = X(67)+1;
            elseif identAA(phosphoPos(jj)-1) == 16
                X(68) = X(68)+1;
            elseif identAA(phosphoPos(jj)-1) == 17
                X(69) = X(69)+1;
            elseif identAA(phosphoPos(jj)-1) == 13 && ...
                    numel(find(modpos==phosphoPos(jj)-1)) > 0
                if numel(strfind(mod_infor(modpos==phosphoPos(jj)-1).name,...
                        'Oxidation')) > 0
                    X(70) = X(70)+1;
                end
            end
        end
    end
    X(71) = min(phosphoPos);
    X(72) = max(phosphoPos);
end