%% Ions and variables calculation
function ionInfo = ions(res_mass,pepseq,charge,mod_infor)
% ion types: y+, b+, a+, (y-18)+, (b-18)+, (y-17)+, (b-17)+, (b-17-17)+,
%                   (y-17-17)+, (b-18-17)+, (y-18-17)+, p+, (p-18)+,
%                   (p-17)+, (p-17-17)+, (p-18-17)+, (p-18-18)+.
% For water loss, it is considered in each type ion and if Asp (D), Glu
% (E), Ser (S), Thr (T) existed, further water loss is added. For ammonia
% loss, it is considered only when residues Arg (R), Lys (K), Asn (N), Gln
% (Q) existed in the peptide sequence, and an additional ammonia loss is
% considered for Arg (R).
% Note: At most 2 neutral losses are allowed for each ion.
% ionInfo: ion mass, ion features, pathway
%

pathname = {'yion';'bion';'proly';'prolb';'asply';'asplb';'gluly';'glulb';'yb2y';'yb2b';...
    'pion';'aion';'bw';'yw';'bn';'yn';'pw';'pn';'pwe';'qn';'bmo';'ymo';'pmo'};

ionInfo = cell2struct(cell(size(pathname)),pathname,1);

numMod = numel(mod_infor);
varNum = 53;
if numMod > 0
    pm = sum(res_mass)+18.015;
else
    pm = molweight(pepseq);
end
pm_avg = molweight(pepseq);
pmz = (pm+charge*1.0079)/charge;

comass = 28.0102;
watermass = 18.0153;
hydrogenmass = 1.0079; ammoniamass = 17.0306;
resnum = numel(pepseq);

% AA identity.
C = 'ARNDCQEGHILKMFPSTWYV';
identAA = zeros(numel(pepseq),1);
for ii = 1:resnum
    identAA(ii) = strfind(C,pepseq(ii));
end
proIdx = strfind(pepseq,'P'); proIdx = proIdx(proIdx~=1);
aspIdx = strfind(pepseq,'D'); aspIdx = aspIdx(aspIdx~=resnum);
gluIdx = strfind(pepseq,'E'); gluIdx = gluIdx(gluIdx~=resnum);

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

bVars = ionfeatures(1:resnum-1,charge,'b',identAA,1,bseries(:,1),pm);
yVars = ionfeatures(2:resnum,charge,'y',identAA,1,yseries(:,1),pm);

if numMod > 0
    bmodVars = ionfeatures_mod(identAA,1:resnum-1,mod_infor,'b');
    ymodVars = ionfeatures_mod(identAA,2:resnum,mod_infor,'y');
end

if numel(proIdx) > 0 % proline effect
    proIonb = bseries(proIdx-1,:);
    ionInfo.('prolb').name = 'prolb';
    ionInfo.('prolb').ionname = bseries_name(proIdx-1);
    ionInfo.('prolb').mass = proIonb;
    ionInfo.('prolb').vars = bVars(proIdx-1,[1:4 6:25 30:varNum]);
    ionInfo.('prolb').adjIdx = [11 12 13 45];
    
    proIony = yseries(proIdx-1,:);
    ionInfo.('proly').name = 'proly';
    ionInfo.('proly').ionname = yseries_name(proIdx-1);
    ionInfo.('proly').mass = proIony;
    ionInfo.('proly').vars = yVars(proIdx-1,[1:4 6:25 30:varNum]);
    ionInfo.('proly').adjIdx = [11 12 13 45];
    if numMod > 0
        ionInfo.('prolb').vars = [bVars(proIdx-1,:) bmodVars(proIdx-1,1:6)];
        ionInfo.('proly').vars = [yVars(proIdx-1,[1:4 6:25 30:varNum]) ...
            ymodVars(proIdx-1,1:6)];
        ionInfo.('prolb').adjIdx = [12 13 14 50];
        % Though this may be not from my expectation, but it can help to
        % construct the model and perform well, thus for this pathway, I
        % used this variable instead.
        sidx = proIdx+3 <= resnum;
        ionInfo.('proly').vars(sidx,46:47) = [ionInfo.('proly').vars(sidx,19) identAA(proIdx(sidx)+3)];
        ionInfo.('proly').vars(~sidx,46:47) = [ionInfo.('proly').vars(~sidx,19) ...
            zeros(numel(find(sidx==0)),1)];
    end
end

if numel(aspIdx) > 0 % aspartic acid effect
    aspIonb = bseries(aspIdx,:);
    ionInfo.('asplb').name = 'asplb';
    ionInfo.('asplb').ionname = bseries_name(aspIdx);
    ionInfo.('asplb').mass = aspIonb;
    ionInfo.('asplb').vars = bVars(aspIdx,[1:3 5:21 26:varNum]);
    ionInfo.('asplb').adjIdx = [11 12 13 45];
    
    aspIony = yseries(aspIdx,:);
    ionInfo.('asply').name = 'asply';
    ionInfo.('asply').ionname = yseries_name(aspIdx);
    ionInfo.('asply').mass = aspIony;
    ionInfo.('asply').vars = yVars(aspIdx,[1:3 5:21 26:varNum]);
    ionInfo.('asply').vars(:,14:15) = ionInfo.('asply').vars(:,14:15)+...
        [ones(numel(aspIdx),1)*-1 ones(numel(aspIdx),1)];
    ionInfo.('asply').adjIdx = [11 12 13 45];
    if numMod > 0
        ionInfo.('asplb').vars = [bVars(aspIdx,:) bmodVars(aspIdx,1:6)];
        ionInfo.('asply').vars = [yVars(aspIdx,[1:3 5:21 26:varNum]) ymodVars(aspIdx,1:7)];
        ionInfo.('asplb').adjIdx = [12 13 14 50];
    end
end

if numel(gluIdx) > 0 % glutamic acid effect
    gluIonb = bseries(gluIdx,:);
    ionInfo.('glulb').name = 'glulb';
    ionInfo.('glulb').ionname = bseries_name(gluIdx);
    ionInfo.('glulb').mass = gluIonb;
    ionInfo.('glulb').vars = bVars(gluIdx,[1:3 5:21 26:varNum]);
    ionInfo.('glulb').adjIdx = [11 12 13 45];
    
    gluIony = yseries(gluIdx,:);
    ionInfo.('gluly').name = 'gluly';
    ionInfo.('gluly').ionname = yseries_name(gluIdx);
    ionInfo.('gluly').mass = gluIony;
    ionInfo.('gluly').vars = yVars(gluIdx,[1:3 5:21 26:varNum]);
    ionInfo.('gluly').adjIdx = [11 12 13 45];
    if numMod > 0
        ionInfo.('glulb').adjIdx = [12 13 14 50];
        ionInfo.('glulb').vars = [bVars(gluIdx,:) bmodVars(gluIdx,1:6)];
        ionInfo.('gluly').vars = [yVars(gluIdx,[1:3 5:21 26:varNum]) ymodVars(gluIdx,:)];
    end
end

%%% b2-yn-2
% y
ionInfo.('yb2y').name = 'yb2y';
ionInfo.('yb2y').ionname = yseries_name(2);
ionInfo.('yb2y').mass = yseries(2,:);
ionInfo.('yb2y').vars = yVars(2,[1:3 5 4 6:19 21 26:29 22:25 30:50 52:varNum]);
ionInfo.('yb2y').adjIdx = [12:14 49];

ionInfo.('yb2b').name = 'yb2b';
ionInfo.('yb2b').ionname = bseries_name(2);
ionInfo.('yb2b').mass = bseries(2,:);
ionInfo.('yb2b').vars = bVars(2,[1:3 5 4 6:19 21 26:29 22:25 30:50 52:varNum]);
ionInfo.('yb2b').adjIdx = [12:14 49];

if numMod > 0
    ionInfo.('yb2y').adjIdx = [12:14 50];
    ionInfo.('yb2b').adjIdx = [12:14 50];
    ionInfo.('yb2b').vars = [bVars(2,:) bmodVars(2,1:6)];
    ionInfo.('yb2y').vars = [yVars(2,:) ymodVars(2,1:6)];
    ionInfo.('yb2y').vars(:,26:29) = ionInfo.('yb2y').vars(:,22:25);
end

% bx-yz
delidx = [proIdx-1 aspIdx gluIdx 2];
bidx = (1:resnum-1)'; bidx(delidx) = [];
yidx = (1:resnum-1)'; yidx(delidx) = [];
if size(bidx,1) > 0
    ionInfo.('bion').name = 'bion';
    ionInfo.('bion').ionname = bseries_name(bidx);
    ionInfo.('bion').mass = bseries(bidx,:);
    ionInfo.('bion').vars = bVars(bidx,[1:25 30:varNum]);
    ionInfo.('bion').adjIdx = [12:14 46];
    
    ionInfo.('yion').name = 'yion';
    ionInfo.('yion').ionname = yseries_name(yidx);
    ionInfo.('yion').mass = yseries(yidx,:);
    ionInfo.('yion').vars = [yVars(yidx,[1:25 30:49 51:varNum]) zeros(size(yidx))];
    ionInfo.('yion').adjIdx = [12:14 49];
    
    if numMod > 0
        ionInfo.('yion').adjIdx = [12:14 50];
        ionInfo.('bion').adjIdx = [12:14 50];
        ionInfo.('bion').vars = [bVars(bidx,:) bmodVars(bidx,1:6)];
        ionInfo.('yion').vars = [yVars(yidx,:) ymodVars(yidx,1:6)];
    end
end

%%% Precursor ions
pionmass = [sum(res_mass)+hydrogenmass+watermass resnum];
ionInfo.('pion').name = 'pion';
ionInfo.('pion').ionname = {'p'};
ionInfo.('pion').mass = pionmass;
ionInfo.('pion').adjIdx = 43;
tempVar = ionfeatures(0,1,'p',identAA,1,pionmass(1),pm);
if numMod > 0
    pmodVars = ionfeatures_mod(identAA,1,mod_infor,'p');
    ionInfo.('pion').vars = [tempVar([1:3 6:9 13 10:12 14 4:5 15:39 41 40 42]) ...
        pmodVars(1:6) tempVar(43)];
    ionInfo.('pion').vars(12) = pionmass(1);
    ionInfo.('pion').adjIdx = 49;
else
    ionInfo.('pion').vars = tempVar;
    ionInfo.('pion').vars(:,end) = charge;
end

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
ionInfo.('aion').mass = [bseries(:,1)-comass bseries(:,2)];
ionInfo.('aion').adjIdx = [12:14 46];
atempVar = bVars;
atempVar(:,12) = atempVar(:,12)-comass;
atempVar(:,13) = atempVar(:,12)/pmz;
atempVar(:,14) = 1-(atempVar(:,12)-1.0079)/pm_avg;
ionInfo.('aion').vars = atempVar(:,[1:25 30:varNum]);

% water loss of y and b ions
% // This section is not correct, now is corrected
ionInfo.('bw').name = 'bw';
ionInfo.('bw').ionname = bwloss_name;
ionInfo.('bw').mass = bwloss;
bwtempVar = bVars;
bwtempVar(:,12) = bwloss(:,1);
bwtempVar(:,13) = bwloss(:,1)/pmz;
bwtempVar(:,14) = 1-(bwloss(:,1)-1.0079)/pm_avg;
ionInfo.('bw').vars = bwtempVar(:,[1:25 30:49 51:varNum 50]);
ionInfo.('bw').adjIdx = [12:14 49];

ionInfo.('yw').name = 'yw';
ionInfo.('yw').ionname = ywloss_name;
ionInfo.('yw').mass = ywloss;
ywtempVar = yVars;
ywtempVar(:,12) = ywloss(:,1);
ywtempVar(:,13) = ywloss(:,1)/pmz;
ywtempVar(:,14) = 1-(ywloss(:,1)-1.0079)/pm_avg;
ionInfo.('yw').vars = ywtempVar(:,[1:25 30:49 51:varNum 50]);
ionInfo.('yw').adjIdx = [12:14 49];

if numMod > 0
    ionInfo.('aion').vars = [atempVar bmodVars(:,1:6)]; % for a ions
    ionInfo.('aion').adjIdx = [12:14 50];
    ionInfo.('bw').vars = [bwtempVar bmodVars(:,1:6)]; % for bw fragments
    ionInfo.('bw').adjIdx = [12:14 50];
    ionInfo.('yw').vars = [ywtempVar ymodVars(:,1:6)]; % for yw fragments
    ionInfo.('yw').adjIdx = [12:14 50];
end

%%% water loss of precursor ions
if strcmp(pepseq(1),'E') % N-terminal glutamine or glutamic acid effect
    thisIonname = 'pwe';
else
    thisIonname = 'pw';
end
ionInfo.(thisIonname).name = thisIonname;
ionInfo.(thisIonname).ionname = {'p-H2O'};
ionInfo.(thisIonname).mass = pwloss;
ionInfo.(thisIonname).vars = ionInfo.('pion').vars;
ionInfo.(thisIonname).adjIdx = ionInfo.('pion').adjIdx;

%%% ammonia loss of y, b and precursor ions
% // This section is not correct, now is corrected
ionInfo.('bn').name = 'bn';
ionInfo.('bn').ionname = bnloss_name;
ionInfo.('bn').mass = bnloss;
bntempVar = bVars;
bntempVar(:,12) = bnloss(:,1);
bntempVar(:,13) = bnloss(:,1)/pmz;
bntempVar(:,14) = 1-(bnloss(:,1)-1.0079)/pm_avg;
ionInfo.('bn').vars = bntempVar(:,[1:25 30:49 51:varNum 50]);
ionInfo.('bn').adjIdx = [12:14 49];

ionInfo.('yn').name = 'yn';
ionInfo.('yn').ionname = ynloss_name;
ionInfo.('yn').mass = ynloss;
yntempVar = yVars;
yntempVar(:,12) = ynloss(:,1);
yntempVar(:,13) = ynloss(:,1)/pmz;
yntempVar(:,14) = 1-(ynloss(:,1)-1.0079)/pm_avg;
ionInfo.('yn').vars = yntempVar(:,[1:25 30:49 51:varNum 50]);
ionInfo.('yn').adjIdx = [12:14 49];

if numMod > 0
    ionInfo.('bn').vars = [bntempVar bmodVars(:,1:6)]; % bn
    ionInfo.('bn').adjIdx = [12:14 50];
    ionInfo.('yn').vars = [yntempVar ymodVars(:,1:6)]; % yn
    ionInfo.('yn').adjIdx = [12:14 50];
end

if strcmp(pepseq(1),'N') || strcmp(pepseq(1),'Q') % N-terminal glutamine or asparagine effect
    ionInfo.('qn').ionname = {'p-NH3'};
    ionInfo.('qn').name = 'qn';
    ionInfo.('qn').mass = pnloss;
    if numMod > 0
        ionInfo.('qn').vars = ionInfo.('pion').vars;
        ionInfo.('qn').adjIdx = 49;
    else
        ionInfo.('qn').vars = ionInfo.('pion').vars(1:39);
        ionInfo.('qn').adjIdx = [];
    end
else
    ionInfo.('pn').name = 'pn';
    ionInfo.('pn').vars = ionInfo.('pion').vars;
    ionInfo.('pn').adjIdx =  ionInfo.('pion').adjIdx;
    ionInfo.('pn').ionname = {'p-NH3'};
    ionInfo.('pn').mass = pnloss;
end

%%%% Modification infor
% Get the information of modified residues.
metoxyidx = zeros(size(mod_infor));
% //
if numMod > 0
    modaapos = cat(1,mod_infor.pos);
    modname = cat(1,{mod_infor.name})';
    for ii = 1:numel(modname)
        currmodname = modname{ii};
        if numel(strfind(lower(currmodname),'oxidation')) > 0
            metoxyidx(ii) = modaapos(ii);
        end
    end
end
% //
metoxyidx(metoxyidx==0) = [];

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
        bstempVar = bVars(min(bmetoxyidx):end,:);
        bstempVar(:,12) = bsloss(:,1);
        bstempVar(:,13) = bsloss(:,1)/pmz;
        bstempVar(:,14) = 1-(bsloss(:,1)-1.0079)/(pm+1.0079);
        bstempVar(:,53) = bstempVar(:,53)-18;
        ionInfo.('bmo').vars = [bstempVar bmodVars(min(bmetoxyidx):end,1:5)];
        ionInfo.('bmo').vars = [ionInfo.('bmo').vars(:,[1:49 51:58]) ...
            ionInfo.('bmo').vars(:,50)];
        ionInfo.('bmo').adjIdx = [12:14 58];
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
        ystempVar = yVars(1:max(ymetoxyidx)-1,:);
        ystempVar(:,12) = ysloss(:,1);
        ystempVar(:,13) = ysloss(:,1)/pmz;
        ystempVar(:,14) = 1-(ysloss(:,1)-1.0079)/(pm+1.0079);
        ystempVar(:,53) = ystempVar(:,53)-18;
        ionInfo.('ymo').vars = [ystempVar ymodVars(1:max(ymetoxyidx)-1,1:5)];
        ionInfo.('ymo').vars = [ionInfo.('ymo').vars(:,[1:49 51:58]) ...
            ionInfo.('ymo').vars(:,50)];
        ionInfo.('ymo').adjIdx = [12:14 58];
    end
    
    psloss = [pionmass(1)-63.9982 resnum];
    ionInfo.('pmo').name = 'pmo';
    ionInfo.('pmo').ionname = {'p-CH3SOH'};
    ionInfo.('pmo').mass = psloss;
    ionInfo.('pmo').vars = ionInfo.('pion').vars;
    ionInfo.('pmo').adjIdx = 49;
end

%%%% ====================================================================
% doubly and triply charged, with doubly and triply charged fragments.
if charge >= 2 && resnum >= 3
    
    for ii = 1:23
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
                
%                 disp(pathname{ii});
                tempVar = ionInfo.(pathname{ii}).vars(c2idx,:);
                adjIdx = ionInfo.(pathname{ii}).adjIdx;
                if ~isempty(adjIdx) % not N-terminal Q or N effect ions
                    
                    tempVar(:,adjIdx(end)) = 2;
                    if numel(adjIdx) > 1 % non precursor ions....
                        tempVar(:,adjIdx(1)) = c2Ionmass(:,1);
                        tempVar(:,adjIdx(2)) = tempVar(:,adjIdx(1))/pmz;
                    end
                    
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
                    adjIdx = ionInfo.(pathname{ii}).adjIdx;
                    if ~isempty(adjIdx) % not N-terminal Q or N effect ions
                        
                        tempVar(:,adjIdx(end)) = 3;
                        if numel(adjIdx) > 1 % non precursor ions....
                            tempVar(:,adjIdx(1)) = c3Ionmass(:,1);
                            tempVar(:,adjIdx(2)) = tempVar(:,adjIdx(1))/pmz;
                        end
                        
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
% Feature extraction for the prediction of tandem mass spectra.
function X = ionfeatures(frag_pos,charge,iontype,identAA,ioncharge,ionmass,precursormass)
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
    X = zeros(numel(frag_pos),53);
    for ii = 1:numel(frag_pos)
        X(ii,1) = identAA(1);
        X(ii,2) = identAA(end);
        X(ii,3) = numel(identAA);
        switch iontype
            case 'y'
                X(ii,4) = identAA(frag_pos(ii)-1); % N-terminal
                X(ii,5) = identAA(frag_pos(ii)); % C-terminal
                X(ii,6) = numel(find(identAA(frag_pos(ii):end)==15));
                X(ii,7) = numel(find(identAA(frag_pos(ii):end)==4));
                X(ii,8) = numel(find(identAA(frag_pos(ii):end)==7));
                X(ii,9) = numel(find(identAA(frag_pos(ii):end)==9))+...
                    numel(find(identAA(frag_pos(ii):end)==12))+...
                    numel(find(identAA(frag_pos(ii):end)==2));
                
                % ...
                if frag_pos(ii) <= X(ii,3)-1 % C-terminal
                    X(ii,19) = identAA(frag_pos(ii)+1);
                end
                if frag_pos(ii) >= 3 % N-terminal
                    X(ii,20) = identAA(frag_pos(ii)-2);
                end
                X(ii,21) = mean(AAdesc.d(identAA(frag_pos(ii):end),3));
                
                % ...
                if frag_pos(ii) >= 4
                    X(ii,51) = identAA(frag_pos(ii)-3);
                end
                if frag_pos(ii) <= X(ii,3)-2
                    X(ii,52) = identAA(frag_pos(ii)+2);
                end
                X(ii,53) = (sum(AAdesc.numatom(identAA(frag_pos(ii):end)))+3)*3-6;
            case {'b','a'}
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
        X(ii,14) = 1-(ionmass(ii)-1.0079)/molweight(AAdesc.aa(identAA));
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
    end
else
    X(1) = identAA(1);
    X(2) = identAA(end);
    X(3) = numel(identAA);
    X(4) = identAA(2);
    X(5) = identAA(end-1);
    X(6) = numel(find(identAA==15));
    X(7) = numel(find(identAA==4));
    X(8) = numel(find(identAA==7));
    X(9) = numel(find(identAA==9))+...
        numel(find(identAA==12))+...
        numel(find(identAA==2));
    for jj = 1:numel(identAA)-1
        if identAA(jj) == 4 && identAA(jj+1) == 15
            X(10) = 1;
            break;
        end
    end
    if identAA(3) == 15; X(11) = 1; end
    if identAA(2) == 4; X(12) = 1; end
    
    if numel(find(identAA==2)) >= charge
        X(13) = 0;
    elseif numel(find(identAA==9))+numel(find(identAA==12)) >= charge
        X(13) = 1;
    else
        X(13) = 2;
    end
    
    X(14) = precursormass+1.0079;
    X(15) = AAdesc.d(X(1),1);
    X(16) = AAdesc.d(X(1),2);
    X(17) = AAdesc.d(X(1),3);
    X(18) = AAdesc.d(X(1),4);
    zaanum = zeros(1,20);
    for jj = 1:20
        zaanum(jj) = numel(find(identAA==jj));
    end
    X(19) = mean(AAdesc.d(identAA,3));
    X(20:39) = zaanum;
    
    X(40) = identAA(3);
    X(41) = identAA(end-2);
    X(42) = (sum(AAdesc.numatom(identAA))+3)*3-6;
    X(43) = ioncharge;
    X(44) = charge;
end

%%%%%%%%%%%
% Feature extraction for the prediction of tandem mass spectra generated by modified peptides.
function X = ionfeatures_mod(~,frag_pos,mod_infor,iontype)

%
npos = numel(frag_pos);
numMod = numel(mod_infor);
modvars = zeros(1,5);
modionvars = zeros(npos,2);
% X = zeros(numel(frag_pos),29);
for ii = 1:numMod
    if numel(strfind(mod_infor(ii).name,'Carbamidomethyl')) > 0
        modvars(1) = modvars(1)+1;
    end
    if numel(strfind(mod_infor(ii).name,'Acetyl')) > 0
        modvars(2) = modvars(2)+1;
    end
    if numel(strfind(mod_infor(ii).name,'Deamidat')) > 0
        modvars(3) = modvars(3)+1;
    end
    if numel(strfind(lower(mod_infor(ii).name),'pyro-glu')) > 0
        modvars(4) = modvars(4)+1;
    end
    if numel(strfind(mod_infor(ii).name,'Oxidation')) > 0
        modvars(5) = modvars(5)+1;
    end
end

mpos = cat(1,mod_infor.pos);
for ii = 1:npos
    switch iontype
        case 'y'
            modionvars(ii,1) = numel(find(mpos>=frag_pos(ii)-2&mpos>=frag_pos(ii)+1));
        case {'b';'a'}
            modionvars(ii,1) = numel(find(mpos>=frag_pos(ii)-1&mpos>=frag_pos(ii)+2));
        case 'p'
            modionvars(ii,1) = numel(find(mpos==1));
    end
    % Identify whether Asp is modified...
    if numel(find(mpos==frag_pos(ii)-1)) > 0
        modionvars(ii,2) = 1;
    end
end

X = [repmat(modvars,npos,1) modionvars];