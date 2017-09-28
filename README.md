
# TMCrys

TMCrys was developed to help the target selection of structural genomics projects by providing prediction
for the propensity of the solubilization, purification and crystallization steps of the crystallization 
process, as well as a prediction for the whole process.


## Citation
If you find it useful, please cite:
Julia Varga and Gábor E. Tusnády  
TMCrys:...

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

### Prerequisites
Packages and modules, copy and paste code below to install.

R packages - from R shell

	```	
	install.packages("xgboost", repos='http://cran.rstudio.com/')
	install.packages("caret", repos='http://cran.rstudio.com/')	
	install.packages("docopt", repos='http://cran.rstudio.com/')	
	install.packages("protr", repos='http://cran.rstudio.com/')
	```


Perl Modules - you may need to add `sudo` before the commands

	```	
	cpan install XML::LibXML	
	cpan install Bio::Tools::Protparam
	cpan install Getopt::Std	
	```

You will also need a modified version of the OB module (used for OB-score calculation), it is downloaded together with TMCrys to tools directory. Please do not remove it.

TMCrys requires an installed copy of BioPerl (http://bioperl.org/INSTALL.html) for running properly. BioPerl could also be installed during the installation of Bio::Tools::Protparam when installer ask about whether to install all modules. 

TMCrys was developed with R v3.4.1 and Perl v5.18.2. Lower versions may not work properly.

### Installing

Download or clone git folder from https://github.com/tmcrys/tmcrys/  
If downloaded as a compressed file, please uncompress it to a folder.

Add $TMCRYS to the environmental variables with

	`export TMCRYS=/path/to/tmcrys/folder`
Or you may copy it to ~/.bashrc or ~/.profile or ~/.bash_profile according to your system settings.
If you want to make it permanent, write `TMCRYS=/path/to/tmcrys/folder to /etc/environment`


### Test
Please run
```
cd $TMCRYS
./tmcrys --test
```

### Usage
For running TMCrys you will need:
1. Sequence and topology of transmembrane protein(s).  There are multiple options for input.
	- A single CCTOP result file, containing one CCTOP entry. Use -i <CCTOPFILE> option with tmcrys.
	- A directory of CCTOP files. Use -d <CCTOPDIR> option.
	- Alternatively, you may also use a space delimited file where lines look as follow: 'proteinID sequence topology). Here, a string represents topology as in 
	You may predict the topology of your protein with CCTOP at http://cctop.enzim.ttk.mta.hu. For multiple proteins, a python script is available at http://cctop.enzim.ttk.mta.hu/?_=/documents/direct_interface.html.
2. netsurfp result .rsa files. Please provide them with -n <NETSURFPFILE> option. It may contain results for or multiple proteins.

For test, these are included in the ./test folder.

## Authors
Julia Varga  
Gábor E. Tusnády

## Troubleshooting
If you encounter any problems, please feel free to open an issue or contact: ....

## License
This project is licensed under the GNU License - see the LICENSE.md file for details.

## References

