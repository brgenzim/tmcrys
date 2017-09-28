
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

R packages - from R shell

	```	
	install.packages("xgboost", dependencies=TRUE, repos='http://cran.rstudio.com/')
	install.packages("caret", dependencies=TRUE, repos='http://cran.rstudio.com/')	
	install.packages("docopt", dependencies=TRUE, repos='http://cran.rstudio.com/')	
	install.packages("protr", dependencies=TRUE, repos='http://cran.rstudio.com/')
	```


Perl Modules - you may need to add `sudo` before the commands

	```	
	cpan install XML::LibXML	
	cpan install Bio::Tools::Protparam
	cpan install Getopt::Std	
	```

and a modified version of the OB module (used for OB score calculation) that can be downloaded together with TMCrys. It requires the installation of BioPerl (http://bioperl.org/INSTALL.html) for running properly. BioPerl could also be installed during the installation of Bio::Tools::Protparam.



### Installing

Download or clone git folder from https://github.com/tmcrys/tmcrys/  
If zipped, please unzip it to a folder.

Add $TMCRYS to the environmental variables with

	`export TMCRYS=/path/to/tmcrys/folder`
Or you may copy it to ~/.bashrc or ~/.profile or ~/.bash_profile according to your system settings.
If you want to make it permanent, write `TMCRYS=/path/to/tmcrys/folder to /etc/environment


### Test
Please run
```
cd $TMCRYS
./tmcrys --test
```
## Authors
Julia Varga  
Gábor E. Tusnády

## Troubleshooting
If you encounter any problems, please feel free to open an issue or contact: ....

## License
This project is licensed under the GNU License - see the LICENSE.md file for details.
