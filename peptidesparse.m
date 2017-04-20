%% Sparse peptides in text files for MS/MS spectra prediction...
function pepinfo = peptidesparse(dirname,filename)
%
% This function is used to sparse peptides stored in text files for peptide
% fragment ion mass spectra prediction. All peptides must be stored as one
% sequence per line with charge state (currently 1-3 charge states are
% allowd) of it in the last column. If you want MS2PBPI to predict MS/MS
% spectra for all charge states (currently 1-3) of a peptide, you must
% assign 0 for it in the last column.
% To predict MS/MS spectra for modified peptides, you can specify the
% modifications for a peptide in two ways:
% 1. Three-column-way: a column specified by modification information
%     must be inserted between peptide sequences and charge states in
%     following format: Mods=#/n,aa,tag/n,aa,tag..., where "#" is the total
%     number of modifications in peptide, "n" is the position of each
%     modification (i.e. the position of modified residue in peptide
%     calculated from N-terminus, starting from 1), "aa" is the one letter
%     abbreviation of modified residue and "tag" is the name of
%     modifications. Thus if no modification, "#" is 0. Multiple
%     modifications are seperated by "/". It should be noted here that no
%     blank or tab space exists in the format line, otherwise error or
%     incorrect MS/MS spectrum is occurred. The name of modification must
%     be specified as follows: 
%           Oxidation: oxidition of metheonine, +16;
%           Carbamidomethyl: carbamidomethylation of cysteine, +57;
%           Acetyl: acetylation, +42;
%           Glu->pyro-Glu: Glu to pyro-Glu, -18;
%           Gln->pyro-Glu: Gln to pyro-Glu, -17;
%           Deamidation: deamidation, +0.98;
%           Phospho: Phosphorylation of serine and threonine, +80;
%     ref: NIST Library of Peptide Ion Fragmentation Spectra, version 2006,
%            at http://chemdata.nist.gov/mass-spc/ftp/mass-spc/PepLib.pdf.
% 2. In-sequence: the second way is specifying modifications in peptide
%     sequences directly. In this way, you can either calculate the residue mass
%     of modified residues or  specify a character to indicate the
%     modification and put them in square brackets adjacent to the modified
%     residues. For example, oxidation of metheonine can be specified as
%     M[147]. MS2PBPI uses this information to calculate m/z values of
%     fragment ions and convert it to structure format which stores the
%     information of modification, i.e. oxidation of metheonine in this
%     example. If N-terminus of the peptide is modified, you must put
%     square brackets to the right site of N-terminal residue.
% Since MS2PBPI only considered five kinds of modifications in MS/MS
% spectra prection, specifying other kinds of modifications or wrong
% modifications (e.g. Carbamidomethyl to metheonine) will obtain a warning
% message and prediction of this line will be skipped. 
% Also, we strongly suggest that you should set your peptide sequence file
% in strict accordance with the above ways we provided. Otherwise warnings
% will be generated.
% However, MATLAB errors may also occur if you do not follow our formats
% due to the unexpected output or input that could not match the require of
% MATLAB's function. If MATLAB's error occurs, this package will be
% terminated without any indication of the peptide sequence line in your
% input text file that causes the error.
% If you provide the peptide information in both ways, we still generate
% MATLAB warning for this because we could not sure what the really
% modifications are, even the information are consistent in both way.
% Though we have tried our best to exhaust all  possible kinds of format
% you may set in the input peptide sequence file for peptide MS/MS spectra
% prediction to run the package without any barrier or sudden termination,
% the unexpected error would still occur, so we strongly suggest that you
% should set your input file in strict accordance with our instruction.
% 
% Nai-ping Dong. PolyU in HK.
% Email: np.dong572@gmail.com
% 10/1/2014

if ischar(filename)
    filename = {filename};
end

totalPeps = 0;
warnPeps = 0;

Allpepinfo = struct('info',cell(size(filename)));

for ii = 1:numel(filename)
    
    % Read the total lines of current file for data initialization
    fid = fopen([dirname filename{ii}]);
    totalLines = numel(cell2mat(textscan(fid,'%1c%*[^\n]')));
    fclose(fid);
    
    h = waitbar(0,'Reading peptide sequences....');
    
    % Peptide sequence information structure initialization
    pepinfo = struct('pepseq',cell(totalLines*2,1),...
        'mod_infor',cell(totalLines*2,1),...
        'charge',cell(totalLines*2,1));
    
    nt = 0;
    d = 1;
    nlines = 0;
    fid = fopen([dirname filename{ii}]);
    n = 0;
    while ~feof(fid)
        nt = nt+1;
        tline = fgetl(fid);
        spaceIdx = find(isspace(tline)==1);
        
        if numel(spaceIdx) == 0 % Bad peptide sequence information
            warning('MATLAB:PepSeqInfoInsufficient',...
                'The infomation for peptide MS/MS spectra seems insufficient in this line: %s',...
                tline);
            continue;
        end
        
        if numel(spaceIdx) >= 1
            deltaSpace = [spaceIdx numel(tline)] - [0 spaceIdx];
            infoIdx = find(deltaSpace>1);
            if deltaSpace(infoIdx(1)) < 5 % The length of peptide sequence must larger than 5
                warning('MATLAB:PepSeqInsufficient',...
                    'The infomation for peptide MS/MS spectra seems insufficient in this line: %s',...
                    tline);
                continue;
            end
        end
        
        % This is probably caused by inserting blanks in the front of the
        % line, whereas no other information such as modifications and
        % charge states are missed...
        if infoIdx(1) > numel(spaceIdx)
            warning('MATLAB:PepSeqInsufficient',...
                'The length of peptide must be larger than 5: %s',...
                tline);
            continue;
        end
        
        % Get peptide sequence (and modifications if are set inside sequence) ....
        pepseq = tline(1:spaceIdx(infoIdx(1)));
        pepseq = pepseq(~isspace(pepseq));
        
        % Charge state
        if numel(tline) - max(spaceIdx) >= 2
            warning('MATLAB:BadChargeInfo',...
                'Current charge state assignment is invalid: %s',...
                tline);
            continue;
        end
        
        if max(spaceIdx) == numel(tline)
            charge = str2double(tline(spaceIdx(infoIdx(end)-1):end));
        else
            charge = str2double(tline(spaceIdx(infoIdx(end)):end));
        end
        if isnan(charge) % Bad peptide sequence line: can not read charge state of current peptide sequence
            warning('MATLAB:BadChargeInfo',...
                'Cannot extract charge state information in this line: %s',...
                tline);
            continue;
        end
        
        if charge > 3 || charge < 0 % Other charge states except 0 and 1-3 are not accepted...
            warning('MATLAB:BadChargeInfo',...
                'Charge of precursor ions must be 0, 1, 2 or 3: %s',...
                tline);
            continue;
        end
            
        if numel(infoIdx) == 1 || (numel(infoIdx)==2 && deltaSpace(infoIdx(2))==2)
            bidx1 = strfind(pepseq,'[');
            bidx2 = strfind(pepseq,']');
            if numel(bidx1) ~= numel(bidx2) % The square brackets must be in pairs...
                warning('MATLAB:BadPeptideSequence',...
                    ['The peptide sequence of current line is bad: %s. \n' ...
                    'This is possibly caused by unrecognized modification formats...'],...
                    tline);
                continue;
            end
            
            if numel(bidx1) == 0 % No modification
                
                if numel(intersect(pepseq,'ARNDCQEGHILKMFPSTWYV')) ~= ...
                        numel(unique(pepseq)) % Unrecognized residue beyond 20 common residues exists...
                    warning('MATLAB:BadPeptideSequence',...
                        'Unrecognized residues existed in current peptide sequence: %s. ', ...
                        tline);
                    continue;
                end
                
                % If charge state is set to 0, the charge is ranging from 1
                % to 3...
                if charge == 0
                    for cc = 1:3
                        n = n+1;
                        pepinfo(n,1).pepseq = pepseq;
                        pepinfo(n,1).charge = cc;
                    end
                else
                    n = n+1;
                    pepinfo(n,1).pepseq = pepseq;
                    pepinfo(n,1).charge = charge;
                end
                nlines = nlines+1;
                
            else
                
                if min(bidx1) == 1 % The modification for N-terminus should be assigned at the right site of this residue
                    warning('MATLAB:BadPeptideSequence',...
                        ['The peptide sequence of current line is bad: %s.\n ' ...
                        'You must assign the modification information of N-terminus residue to the right site of it...'],...
                        tline);
                    continue;
                end
                
                delidx = false(size(pepseq));
                modcheck = false(size(bidx1));
                for kk = 1:numel(bidx1)
                    delidx(bidx1(kk):bidx2(kk)) = true;
                    modcheck(kk) = ~isnan(str2double(pepseq(bidx1(kk)+1:bidx2(kk)-1)));
                end
                
                if ~any(delidx) % I think this is redundant...
                    warning('MATLAB:BadPeptideSequence',...
                        ['The peptide sequence of current line is bad: %s. \n' ...
                        'This is possibly caused by unrecognized modifications...'],...
                        tline);
                    continue;
                end
                
                if ~any(modcheck) % Could not get the mass of modified residue...
                    warning('MATLAB:BadPeptideSequence',...
                        'The modification information inserted in peptide sequence could not be recognized: %s. ',...
                        tline);
                    continue;
                end
                
                mod_infor = struct('aa',cell(size(bidx1')),...
                    'pos',cell(size(bidx1')),...
                    'name',cell(size(bidx1')));
                modpos = getposition(bidx1,bidx2);
                [modname,isValidmod] = getmodname(pepseq,bidx1,bidx2);
                
                if ~isValidmod % Unrecognized modifications
                    warning('MATLAB:BadModification',...
                        'Unrecognized modifications existed: %s. ',...
                        tline);
                    continue;
                end
                
                for kk = 1:numel(bidx1)
                    mod_infor(kk).aa = pepseq(bidx1(kk)-1);
                    mod_infor(kk).pos = modpos(kk);
                    mod_infor(kk).name = modname{kk};
                end
                
                if charge == 0
                    for cc = 1:3
                        n = n+1;
                        pepinfo(n,1).pepseq = pepseq(~delidx);
                        pepinfo(n,1).charge = cc;
                        pepinfo(n,1).mod_infor = mod_infor;
                    end
                else
                    n = n+1;
                    pepinfo(n,1).pepseq = pepseq(~delidx);
                    pepinfo(n,1).charge = charge;
                    pepinfo(n,1).mod_infor = mod_infor;
                end
                nlines = nlines+1;
            end
            
        elseif numel(infoIdx) == 2 || (numel(infoIdx) == 3 && deltaSpace(infoIdx(3))==2)
            
            if numel(intersect(pepseq,'ARNDCQEGHILKMFPSTWYV')) ~= ...
                    numel(unique(pepseq)) % Unrecognized residue beyond 20 common residues exists...
                warning('MATLAB:BadPeptideSequence',...
                    'Unrecognized residues existed in current peptide sequence: %s. ', ...
                    tline);
                continue;
            end
            
            modStr = tline(spaceIdx(infoIdx(1)):spaceIdx(infoIdx(2)));
            modStr = modStr(~isspace(modStr));
            
            [mod_infor,isValidmod] = getmodinfo(modStr);
            if ~isValidmod % Unrecognized modifications
                warning('MATLAB:BadModification',...
                    'Unrecognized modifications existed: %s. ',...
                    tline);
                continue;
            end
            
            if charge == 0
                for cc = 1:3
                    n = n+1;
                    pepinfo(n,1).pepseq = pepseq;
                    pepinfo(n,1).charge = cc;
                    pepinfo(n,1).mod_infor = mod_infor;
                end
            else
                n = n+1;
                pepinfo(n,1).pepseq = pepseq;
                pepinfo(n,1).charge = charge;
                pepinfo(n,1).mod_infor = mod_infor;
            end
            nlines = nlines+1;
            
        else
            warning('MATLAB:BadPeptideSequenceInfo',...
                'Unrecognized information in this line: %s. ', ...
                tline);
        end
        
        % Show the progress in percentage...
        nPercent = nt/totalLines*100;
        if mod(round(nPercent),10) == 0 && nPercent-floor(nPercent/10)*10 < d ...
                && nPercent > 10
            currPercent = round(nPercent)/10;
%             disp([num2str(round(nPercent)) '% of peptide sequences of ' filename{ii} ...
%                 ' have been processed...']);
            waitbar(currPercent/10,h,...
                sprintf('%d%% of peptide sequences of "%s" \nhave been processed...',...
                currPercent*10,filename{ii}));
        end
        d = nPercent-floor(nPercent/10)*10;
        
    end
    fclose(fid);
    delete(h);
    
    Allpepinfo(ii,1).info = pepinfo(1:n);
    totalPeps = totalPeps+n;
    warnPeps = warnPeps+totalLines-nlines;
end

pepinfo = cat(1,Allpepinfo.info);
clear Allpepinfo;
fprintf('Total %d peptides have been input with %d warnings...\n',...
    totalPeps,warnPeps);

%% Get position of modification in a peptide sequence
function pos = getposition(bracketIdx1,bracketIdx2)

if numel(bracketIdx1) == 1
    pos = bracketIdx1-1;
    return
end

pos = zeros(size(bracketIdx1));
pos(1) = bracketIdx1(1)-1;

for ii = 1:numel(bracketIdx2)-1
    pos(ii+1) = pos(ii)+bracketIdx1(ii+1)-bracketIdx2(ii)-1;
end


%% Get names of modification in a peptide sequence
function [modname,isValidmod] = getmodname(pepseq,bracketIdx1,bracketIdx2)
%
% Modifications and corresponding residues modified considered currently
% are:
% N-terminal glutamic acid and glutamine ------ pyro-Glu (-18, -17)
% metheonine ------ oxidation (+16)
% serine and threonine ------- phosphorylation (+80)
% N-terminal amino acid ------ acetylation (+42)
% asparagine and glutamine ------- deamidation (+0.98)
% cysteine ------- carbamidomethylation (+57)

mAA.C = 103.1388;
mAA.E = 129.1155;
mAA.Q = 128.1307;
mAA.S = 87.0782;
mAA.T = 101.1051;
mAA.M = 131.1926;
mAA.A = 71.0788;
mAA.D = 115.0886;
mAA.F = 147.1766;
mAA.G = 57.0519;
mAA.H = 137.1411;
mAA.I = 113.1594;
mAA.L = 113.1594;
mAA.K = 128.1741;
mAA.N = 114.1038;
mAA.P = 97.1167;
mAA.R = 156.1875;
mAA.V = 99.1326;
mAA.W = 186.2132;
mAA.Y = 163.176;

modname = cell(size(bracketIdx1));
isValidmod = true;

for ii = 1:numel(bracketIdx1)
    
    modM = str2double(pepseq(bracketIdx1(ii)+1:bracketIdx2(ii)-1));
    
    if ~isnan(modM)
        deltaMmod = modM-mAA.(pepseq(bracketIdx1(ii)-1));
        
        if strcmp(pepseq(bracketIdx1(ii)-1),'M') && abs(deltaMmod-16) < 5
            modname{ii} = 'Oxidation';
        elseif strcmp(pepseq(bracketIdx1(ii)-1),'C') && abs(deltaMmod-57) < 5
            modname{ii} = 'Carbamidomethyl';
        elseif (strcmp(pepseq(bracketIdx1(ii)-1),'N') || strcmp(pepseq(bracketIdx1(ii)-1),'Q')) ...
                && abs(deltaMmod-0.98) < 5
            modname{ii} = 'Deamidation';
        elseif (strcmp(pepseq(bracketIdx1(ii)-1),'T') || strcmp(pepseq(bracketIdx1(ii)-1),'S') || ...
                strcmp(pepseq(bracketIdx1(ii)-1),'Y')) && abs(deltaMmod-80) < 5
            modname{ii} = 'Phospho';
        elseif bracketIdx1(ii)-1 == 1
            if strcmp(pepseq(1),'E') && abs(deltaMmod+18) < 5
                modname{ii} = 'Glu->pyro-Glu';
            elseif strcmp(pepseq(1),'Q') && abs(deltaMmod+17) < 5
                modname{ii} = 'Gln->pyro-Glu';
            elseif abs(deltaMmod-42) < 5
                modname{ii} = 'Acetyl';
            else
                isValidmod = false;
            end
        else
            isValidmod = false;
        end
        
    else
        
        isValidmod = false;
        
    end
end


%% Get modification information from the string...
function [mod_infor,isValidmod] = getmodinfo(modtext)
%
%

mod_name = ['oxidation;' ...
    'carbamidomethyl;' ...
    'acetyl;' ...
    'deamidation;' ...
    'gln->pyro-glu;' ...
    'glu->pyro-glu;' ...
    'phospho;'];
modAA = 'M;C;ARNDCQEGHILKMFPSTWYV;NQ;Q;E;STY;';
semicolonidxName = [0 strfind(mod_name,';')];
semicolonidx = [0 strfind(modAA,';')];

slashidx = strfind(modtext,'/');

if numel(slashidx) == 0 
    isValidmod = ~isnan(str2double(modtext(end)));
    mod_infor = [];
    return;
end

if numel(strfind(modtext,'Mods=')) == 0 % In this way, modification must be started with "Mods="
    isValidmod = false;
    mod_infor = [];
    return;
end

nummod = str2double(modtext(slashidx(1)-1));
if isnan(nummod) || ... % Could not recognize the number
        numel(slashidx) ~= nummod
    isValidmod = false;
    mod_infor = [];
    return;
end

if nummod == 0
    isValidmod = true;
    mod_infor = [];
    return;
end

mod_infor = struct('aa',cell(nummod,1),...
    'pos',cell(nummod,1),...
    'name',cell(nummod,1));

isValidmod = true;
commaidx = [strfind(modtext,',') numel(modtext)];
for ii = 1:nummod-1
    mod_infor(ii,1).aa = modtext(commaidx((ii-1)*2+1)+1);
    mod_infor(ii,1).name = modtext(commaidx(ii*2)+1:slashidx(ii+1)-1);
    mod_infor(ii,1).pos = str2double(modtext(slashidx(ii)+1:commaidx((ii-1)*2+1)-1));
end

mod_infor(nummod,1).aa = modtext(commaidx((nummod-1)*2+1)+1);
mod_infor(nummod,1).name = modtext(commaidx(nummod*2)+1:end);
mod_infor(nummod,1).pos = str2double(modtext(slashidx(nummod)+1:commaidx((nummod-1)*2+1)-1));

% Checking...
for ii = 1:nummod
    if numel(strfind(mod_name,lower(mod_infor(ii).name))) == 0 % unrecognized modification
        isValidmod = false;
    elseif numel(strfind(mod_name,lower(mod_infor(ii).name))) == 2 % Possibly occurs when setting modification name as pyro-glu
        if strcmp(mod_infor(ii).aa,'E')
            mod_infor(ii).name = 'Glu->pyro-Glu';
        elseif strcmp(mod_infor(ii).aa,'Q')
            mod_infor(ii).name = 'Gln->pyro-Glu';
        else
            isValidmod = false;
        end
    else
        modidx = strfind(mod_name,lower(mod_infor(ii).name));
        scidx = find(semicolonidxName<modidx,1,'last');
        modAAref = modAA(semicolonidx(scidx)+1:semicolonidx(scidx+1));
        if numel(strfind(modAAref,mod_infor(ii).aa)) == 0 % incorrect modification assignments...
            isValidmod = false;
        end
    end
end