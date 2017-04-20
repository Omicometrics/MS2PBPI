function residuemass = residumasscal(peptide_seq,mod_infor,mwtype)
%
% Calculate each residue mass in given peptide sequence (input variable
% peptide_seq) and modification information with molecular weight type.
% mwtype: molecular weight type, 1 = average molecular weight, 0 =
%                 monoisotope molecular weight.
%
% Naiping Dong, College of Chemistry and Chemical Engineering, Central
% South University, Hunan, Changsha, R. P. China, 410083.
% 2011/11/27

residuemass = zeros(numel(peptide_seq),1);

% delete oxidation modification of methionine
oidx = strfind(peptide_seq,'O');
if numel(oidx) ~= 0
    peptide_seq([oidx oidx-1 oidx+1]) = [];
end

pepmw.monoisoweight = zeros(numel(peptide_seq),1); % monoisotope molecular weight.
pepmw.averageweight = zeros(numel(peptide_seq),1); % average molecular weight.

for ii = 1:numel(peptide_seq)
    
    % one letter code of amino acid.
    resstrc = lower(peptide_seq(ii));
    switch resstrc
        case 'a'
            resmw_avg = 71.0788;
            resmw_mono = 71.03712;
        case 'c'
            resmw_avg = 103.1388;
            resmw_mono = 103.00919;
        case 'd'
            resmw_avg = 115.0886;
            resmw_mono = 115.02695;
        case 'e'
            resmw_avg = 129.1155;
            resmw_mono = 129.04260;
        case 'f'
            resmw_avg = 147.1766;
            resmw_mono = 147.06842;
        case 'g'
            resmw_avg = 57.0519;
            resmw_mono = 57.02147;
        case 'h'
            resmw_avg = 137.1411;
            resmw_mono = 137.05891;
        case 'i'
            resmw_avg = 113.1594;
            resmw_mono = 113.08407;
        case 'k'
            resmw_avg = 128.1741;
            resmw_mono = 128.09497;
        case 'l'
            resmw_avg = 113.1594;
            resmw_mono = 113.08407;
        case 'm'
            resmw_avg = 131.1926;
            resmw_mono = 131.04049;
        case 'n'
            resmw_avg = 114.1038;
            resmw_mono= 114.04293;
        case 'p'
            resmw_avg = 97.1167;
            resmw_mono = 97.05277;
        case 'q'
            resmw_avg = 128.1307;
            resmw_mono = 128.05858;
        case 'r'
            resmw_avg = 156.1875;
            resmw_mono= 156.10112;
        case 's'
            resmw_avg = 87.0782;
            resmw_mono = 87.03203;
        case 't'
            resmw_avg = 101.1051;
            resmw_mono = 101.04768;
        case 'v'
            resmw_avg = 99.1326;
            resmw_mono = 99.06842;
        case 'w'
            resmw_avg = 186.2132;
            resmw_mono = 186.07932;
        case 'y'
            resmw_avg = 163.176;
            resmw_mono = 163.06333;
    end
    pepmw.monoisoweight(ii) = resmw_mono;
    pepmw.averageweight(ii) = resmw_avg;
end

% residue mass
switch mwtype
    case 1
        residuemass = pepmw.averageweight;
    case 0
        residuemass = pepmw.monoisoweight;
end

mod_name = ['oxidation;' ...
    'carbamidomethyl;' ...
    'icat_light;' ...
    'icat_heavy;' ...
    'ab_old_icatd0;' ...
    'ab_old_icatd8;' ...
    'acetyl;' ...
    'deamidation;' ...
    'pyro-cmc;' ...
    'gln->pyro-glu;' ...
    'glu->pyro-glu;' ...
    'amide;' ...
    'phospho;' ...
    'methyl;' ...
    'carbamyl;' ...
    'pyro-carbamidomethyl;' ...
    'propionamide'];
mod_mass = [15.9949;57.0215;227.12;236.12;442.20;450.20;42.0106;0.9840;-17.0265;-17.027;-18.015;-0.9840;79.9663;14.0157;43.0058;40.0125;73.09];

semicolonidx = [0 strfind(mod_name,';')];

% modification information
if numel(mod_infor) ~= 0
    for ii = 1:numel(mod_infor)
        idx = isspace(mod_infor(ii).name)==1;
        mod_infor(ii).name(idx) = [];
        pos = mod_infor(ii).pos;
        mod_name_pos = strfind(mod_name,lower(mod_infor(ii).name));
        mod_name_idx = find(mod_name_pos(1)-semicolonidx>0,1,'last');
        residuemass(pos) = residuemass(pos)+mod_mass(mod_name_idx);
    end
end
