
# TMCrys

TMCrys was developed to help the target selection of structural genomics projects by providing prediction
for the propensity of the solubilization, purification and crystallization steps of the crystallization 
process, as well as a prediction for the whole process.

## Citation
If you find TMCrys useful, please cite:
Julia Varga and Gábor E. Tusnády  
TMCrys:...

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites
Packages and modules, please copy and paste code below to install.  
TMCrys was developed with R v3.4.1 and Perl v5.18.2. Lower versions may not work properly.

R packages - from R shell

```	
install.packages("xgboost", repos='http://cran.rstudio.com/')
install.packages("caret", repos='http://cran.rstudio.com/')	
install.packages("docopt", repos='http://cran.rstudio.com/')	
install.packages("protr", repos='http://cran.rstudio.com/')
```


Perl Modules - you may need to add `sudo` before the commands.

```
cpan install XML::LibXML
cpan install Bio::Tools::Protparam
cpan install Getopt::Std
```

You will also need a modified version of the OB module (used for OB-score calculation), it is downloaded together with TMCrys to tools directory. Please do not remove it or data (data/zmat.dat) belonging to it.

TMCrys requires an installed copy of BioPerl (http://bioperl.org/INSTALL.html) for running properly. BioPerl could also be installed during the installation of Bio::Tools::Protparam when installer ask about whether to install all modules. 


### Installing

Download or clone git folder from https://github.com/tmcrys/tmcrys/.
```
git clone https://github.com/tmcrys/tmcrys/
```
If downloaded as a compressed file, please uncompress it to a folder.

Add $TMCRYS to the environmental variables with

```
export TMCRYS=/path/to/tmcrys/folder  
```
Alternatively, you may copy it to ~/.bashrc or ~/.profile or ~/.bash_profile according to your system settings.
If you want to make it permanent, write `TMCRYS=/path/to/tmcrys/folder` to /etc/environment.


### Test
Please run
```
cd $TMCRYS
./tmcrys --test
```

### Usage

For running TMCrys you will need:
1. Sequence and topology of transmembrane protein(s).  There are multiple options for input.
	- A single CCTOP result file, containing one CCTOP entry. Use `-i <CCTOPFILE>` option with tmcrys.
	- A directory of CCTOP files. Use `-d <CCTOPDIR>` option. If you use this possibility, please provide a name for the project with `--name NAME` option.
	- Alternatively, you may also use a space delimited file where lines look as follow: 'proteinID sequence topology). Here, a string represents topology as in `test/test.txt` file. Use `-s <DELIMITEDFILE>` option.
	You may predict the topology of your protein with CCTOP at http://cctop.enzim.ttk.mta.hu. For multiple proteins, a python script is available at http://cctop.enzim.ttk.mta.hu/?_=/documents/direct_interface.html.
2. NetSurfP result .rsa files. Please provide them with `-n <NETSURFPFILE>` option. It may contain results for or multiple proteins. NetSurfP may be run or downloaded from http://www.cbs.dtu.dk/services/NetSurfP/.
3. A working directory, specified with `--wd <DIR>` option.

For test purposes, all these are included in the ./test folder.

To run please type
```
cd $TMCRYS
./tmcrys (-i <CCTOPFILE> | -d <CCTOPDIR> | -s <DELIMITEDFILE>) -n <NETSURFPFILE> --wd <DIR>
```

Help for every script is available by typing `-h` or `--help` or no arguments when running commands.

## Authors
Julia Varga  
Gábor E. Tusnády

## Troubleshooting
If you encounter any problems, please feel free to open an issue or contact: (email@ttk.mta.hu)

## License
This project is licensed under the GNU License - see the LICENSE.md file for details.

## References


1. Dobson,L., Reményi,I. and Tusnády,G.E. (2015) CCTOP: a Consensus Constrained TOPology prediction web server. Nucleic Acids Res., 43, W408–W412.
2. Xiao,N., Cao,D.-S., Zhu,M.-F. and Xu,Q.-S. (2015) protr/ProtrWeb: R package and web server for generating various numerical representation schemes of protein sequences. Bioinformatics, 31, 1857–1859.
3. Walker,J.M. ed. (2005) The Proteomics Protocols Handbook Humana Press, Totowa, NJ.
4. Overton,I.M. and Barton,G.J. (2006) A normalised scale for structural genomics target ranking: The OB-Score. FEBS Lett., 580, 4005–4009.
5. Petersen,B., Petersen,T.N., Andersen,P., Nielsen,M. and Lundegaard,C. (2009) A generic method for assignment of reliability scores applied to solvent accessibility predictions. BMC Struct. Biol., 9, 51.
6. Chen,T., He,T., Benesty,M., Khotilovich,V. and Tang,Y. (2017) xgboost: Extreme Gradient Boosting.
7. Kuhn,M. et al. (2017) caret: Classification and Regression Training.
8. Kawashima,S., Ogata,H. and Kanehisa,M. (1999) AAindex: Amino Acid Index Database. Nucleic Acids Res., 27, 368–369.
9. Yan,Y. (2016) rBayesianOptimization: Bayesian Optimization of Hyperparameters.

