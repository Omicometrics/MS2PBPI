%% Main function of MS2PBPI
function ms2pbpi
% This is the main function of mass spectral prediction framework MS2PBPI
% (tandem mass spectra prediction boosting peptide identification).
%
% Naiping Dong. PolyU HK
% Email: np.dong572@gmail.com
% 1/23/2014


[filename,dirname] = uigetfile({'*.txt','Text File(*.txt)'},...
    'Select Peptide Sequence InfoFile',...
    'MultiSelect', 'on');
if isnumeric(dirname)
    return;
end

pepinfo = peptidesparse(dirname,filename);
n = numel(pepinfo);
if n == 0
    disp('No peptide input! Program is terminated...');
    return;
end

ms = ms2sparse(pepinfo);
t = 10000;
nf = ceil(n/t);

disp('Write predicted mass spectra into ..mgf files');
if nf == 1
    mgfpredwriter(ms,pepinfo,[dirname 'ms2pbpi_MS.mgf']);
else
    for ii = 1:nf
        mgfpredwriter(ms(t*(ii-1)+1:min(ii*t,n)),...
            pepinfo(t*(ii-1)+1:min(ii*t,n)),...
            [dirname 'ms2pbpi_MS_' int2str(ii) '.mgf']);
    end
end