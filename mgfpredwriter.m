function mgfpredwriter(ms2data,pepinfo,filename)
%
% Writing predicted mass spectra into Mascot Generic Format (.mgf) file.
% The introduction of this format file can be found at:
% "http://www.thegpm.org/GPM/faq.html#faq4" and
% "http://www.matrixscience.com/help/data_file_help.html#GEN"
% "filename" must contain the full path.
%
% Note:
% The default set of charge states for each spectrum is 1 to 3, if not
% specified by the spectrum such as introdued by dta files.
% Naiping Dong. PolyU HK
% Email: np.dong572@gmail.com
% 23/1/2014

n = numel(pepinfo);
fid = fopen(filename,'a+');
fprintf(fid,'%s\n','GENERATOR=MS2PBPI MASS SPECTRAL Predictor');
fprintf(fid,'%s\n','REPTYPE=Peptide');

for ii = 1:n
    
    fprintf(fid,'%s\n','BEGIN IONS');
    
    peaks = ms2data{ii};
    
    res_mass = residumasscal(pepinfo(ii).pepseq,pepinfo(ii).mod_infor,1);
    precursormz = (sum(res_mass)+18.015+pepinfo(ii).charge*1.0079)/pepinfo(ii).charge; % Precursor m/z
    fprintf(fid,'%s%f\n','PEPMASS=',precursormz);
    
    fprintf(fid,'%s%d%s\n','CHARGE=',pepinfo(ii).charge,'+');
    fprintf(fid,'%s\n','TITLE=MS2PBPI Predicted Mass Spectra');
    
    peaks(:,1) = round(peaks(:,1)*10000)/10000;
    dlmwrite(filename,peaks,'delimiter', '\t','precision','%.4f','-append');
    
    fprintf(fid,'%s\n','END IONS');
    fprintf(fid,'%s\n',''); % Leave one for blank line to seperate different MS2 query
    
end
fclose(fid);