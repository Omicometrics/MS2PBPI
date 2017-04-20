#MS2PBPI#
---------------------
**MS2PBPI** predicts peptide tandem mass spectra generated from collision induced dissociation (CID) for facilitation of peptide identificaition in shotgun proteomics. It is coded by MATLAB. The prediction algorithm is based on stochastic gradient boosting trees (SGBTree) ([Friedman J H, 1999](https://statweb.stanford.edu/~jhf/ftp/stobst.pdf)), and trained by large number of tandem mass spectra, thus can provide accurate prediction of mass spectra.  
Currently, _MS2PBPI_ is only implemented by MATLAB. In order to provide more easier way to reimplement by other programming language if is required, text versions of each SGBTree model as well as binary tree is also generated and shared in the same compressed file .

#Usage#
---------------
Simply type `ms2pbpi` in MATLAB command line, then a dialog will be raised to select file that stores the peptide sequences for prediction. To avoid `Undefined function or variable 'ms2pbpi'` error, put all `.m` files and folder `models` under working directory or set the path of MS2PBPI to the search path of MATLAB by `path` command or `Set Path` dialog.    
Once the prediction is finished, an `Mascot generic format (.mgf)` is generated to store all the predicted tandem mass spectra.

**Note** MS2PBPI firstly generates fragments and corresponding variables for each peptide sequence input, and then analyzes the variables to assign regions to these fragments. Finally the package loads SGBTree model to predict the intensities. This procedure will take lots of time, especially for loading regression models. Thus we generate fragments and variables for all input peptide sequences and attempt to load each regression model only once. Considering this, we suggest that all peptide sequences should be input via  a single file. However, much more memory is required for predicting tandem mass spectra for all sequences. So, consider the available memory of the computer, input as more peptide sequences as possible, but should not exceed the memory, or `out of memory error` will be raised.

#Input File Format  
-------------------------
The input file can be any format but should be readable as a text (i.e., ASCII). All peptides must be stored as one sequence per line with charge state (currently 1-3 charge states are allowd) at the last column. If one wants _MS2PBPI_ to predict tandem mass spectra for all charge states (currently 1-3) for a peptide sequence, **0** should be assigned at last column instead.
To predict tandem mass spectra for modified peptides, one can specify the modifications for a peptide in two ways:  

##Three-column
A column specified by modification information must be inserted between peptide sequences and charge states in following format: `Mods=#/n,aa,tag/n,aa,tag...`, where `#` is the total number of modifications in peptide, `n` is the position of each modification (i.e. the position of modified residue in peptide calculated from N-terminus, started by 1), `aa` is the one letter abbreviation of modified residue and `tag` is the name of modification. Thus if there is no modification, `#` is 0. Multiple modifications are separated by `/`. It should be noted here that no blank or tab space exists in the format line, otherwise error or tandem mass spectrum can be incorrectly predicted. The name of modification must be specified as follows:  
>Oxidation: oxidation of methionine, +16;  
Carbamidomethyl: carbamidomethylation of cysteine, +57;  
Acetyl: acetylation, +42;  
Glu->pyro-Glu: Glu to pyro-Glu, -18;  
Gln->pyro-Glu: Gln to pyro-Glu, -17;  
Deamidation: deamidation, +0.98;  
Phospho: Phosphorylation of serine and threonine, +80;  

([NIST Library of Peptide Ion Fragmentation Spectra, version 2006](http://chemdata.nist.gov/mass-spc/ftp/mass-spc/PepLib.pdf))

##In-sequence
The second way is specifying modifications in peptide sequences directly. In this way, one can either calculate the residue mass of modified residues or specify a character to indicate the modifications and put them in square brackets adjacent to the modified residues. For example, oxidation of methionine can be specified as `M[147]`. _MS2PBPI_ uses this information to calculate m/z values of fragment ions and convert it to structure format which stores the information of modifications, e.g. oxidation of methionine in this example. If N-terminus of the peptide is modified, one must put square bracket to the right site of N-terminal residue. The recognized modifications by MS2PBPI can be expressed as follows:
>Oxidation: M[147];  
Carbamidomethyl: C[160];  
Acetyl: AA[mAA+42];  
Glu->pyro-Glu: E[111];  
Gln->pyro-Glu: Q[111];  
Deamidation: AA[mAA+0.98];  
Phospho: S[167], Y[243];  

Where `AA` represents any of 20 amino acid residues, `mAA` stands for the mass of corresponding residue.

**Note** Since _MS2PBPI_ currently only considers five kinds of modifications in tandem mass spectrum prediction, specifying other kinds of modifications or wrong modifications (e.g. Carbamidomethyl to methionine) will obtain a warning message and the prediction of this sequence will be skipped. MATLAB errors may also occur if the sequences do not follow our formats due to the unexpected input that could not match the requirement of MATLAB's function. If MATLAB's error occurs, this package will be terminated without any indication of the peptide sequence in the input text file that causes the error.
If peptide modifications are provided in both ways, we still generate MATLAB warning for this because we could not sure what the really modifications are, even the information are consistent between both ways. Though we have tried our best to exhaust all  possible kinds of format one may set in the input peptide sequence file for peptide tandem mass spectra prediction to run the package without any barrier or sudden termination, the unexpected error would still occur.
**So, we strongly suggest that the peptide sequence file should be input in strict accordance with the above ways we provided or example file in the package.**

#Citation
-----------------
`Dong NP, Liang YZ, Xu QS, et al. Prediction of Peptide Fragment Ion Mass Spectra by Data Mining Techniques. Anal Chem. 2014, 86(15), 7446-7454.`  

#Copy Right
-----------------
&copy; Naiping Dong  
Email: np.dong572@gmail.com

#License
----------
_MS2PBPI_ is distributed under _Apache License 2.0_.
