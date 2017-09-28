#!/usr/bin/env Rscript

####################
#
#    TMCrys version 0.1
#    Copyright 2017 Julia Varga and Gábor E. Tusnády
#
#    If you use TMCrys, please cite: 
#    Julia Varga and Gábor E. Tusnády TMCRys...
#
#
#    This file is part of TMCrys.
#
#    TMCrys is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    TMCrys is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with TMCrys.  If not, see <http://www.gnu.org/licenses/>.
#
####################

options(verbose = F, warn = -1, width = 200, show.error.locations=TRUE)
library(caret, verbose = F, quietly = T, warn.conflicts = F)
library(xgboost, verbose = F, quietly = T, warn.conflicts = F)
library(docopt, quietly = T, warn.conflicts = F, verbose = F)
'Usage:
   prediction.R (-i FILE) [-o FILE] [--verbose] [--Rdata RFILE] [--version]

Options:
   -h --help  	print this help
   -i FILE  	file with calculated features
   -o FILE  	output file, if not specified, results will be print to STDOUT
   --verbose  	print result to STDOUT besides writing to file
   --Rdata  	path to tmcrys.Rdata if not in the same directory as script or was renamed [default: ./tmcrys.Rdata]
   --version   print version and exit

TMCrys version 0.1
If you use it, please cite Julia Varga and Gábor E. Tusnády TMCrys... 
' -> doc

opts <- docopt(doc, help = T, version = "TMCrys version 0.1\n")
if (opts$Rdata == TRUE ){
	if (file.exists(opts$RFILE)) {
		load(opts$RFILE)
	}
	else{
		cat(paste("Could not find ", opts$RFILE, " please specifiy path with --Rdata RFILE option\n", sep = ""))
		quit(save="no", status = 1 )
	}
} else{
	if (file.exists("tmcrys.Rdata")) {
		load("tmcrys.Rdata")
	}
	else{
		cat(paste("Could not find ", "tmcrys.Rdata", " please specifiy path with --Rdata RFILE option\n", sep = ""))
		quit(save="no", status = 1 )
	}
}
# opts <- list()
# opts$i <- "/bigdisk/users/jvarga/crys/tools/vegso/test/testxml.features"
if (file.exists(opts$i)){
	data <- read.table(opts$i, header = T, row.names = 1, sep = "\t")
} else {
	cat(paste("Could not find input file ", opts$i, " please specifiy correct path with -i FILE option\n", sep = ""))
	quit(save="no", status = 1 )
}

row.has.na <- apply(data, 1, function(x){any(is.na(x))})
omit <- data[row.has.na,]
data <- data[!row.has.na,]

data1 <-  data[,which(colnames(data) %in% features1)]
data2 <-  data[,which(colnames(data) %in% features2)]
data3 <-  data[,which(colnames(data) %in% features3)]


#1st step - solubilization
dataOut1 <- data.frame(rownames(data1))
for (i in 2:ncol(data1)){
	col <- data.frame(data1[,c(1,i)])
	colname <- colnames(col)[2]
	if ( colname == "pI"){
		col$dipI <- NA
		col$dipI[col$pI > 7 ] <- 1
		col$dipI[col$pI <= 7 ] <- 0
		dataOut1 <- cbind(dataOut1 , col$dipI)
	} 	else if ( grepl("lag", colname)  || 
			   colname == "Xc2.lambda.1" || 
			   colname == "Xc2.lambda.25" || 
			   colname == "Xc2.lambda.30" ||
			   colname == "Xc2.lambda.10" ||
			   colname == "Xc2.lambda.15" ||
			   colname == "Xc2.lambda.2" ||
			   colname == "Xc2.lambda.20" ||
			   colname == "Xc2.lambda.3" 
	) {
		
		lowerBound <- dev1[colname,"lowerBound"]
		upperBound <- dev1[colname,"upperBound"]
		
		if(is.na(lowerBound) || is.na(upperBound) ){
			next
		}
		
		# print(paste(colname,lowerBound, upperBound, colname, collapse = "\t"))
		
		col[col[,2] <= lowerBound | col[,2] >= upperBound, 2] <- 1
		col[col[,2] > lowerBound & col[,2] < upperBound, 2] <- 0
		
		colsum <- sum(col[,2])
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	} 	else if (	colname %in% vector1) {
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (grepl("Xc1", colname) || 
			 colname == "TM.C" || 
			 colname == "TM.D" || 
			 colname == "TM.H" ||
			 colname == "TM.K" ||
			 colname == "buriedratio" ||
			 colname == "length" ||
			 colname == "lengthNonTM" ||
			 colname == "lengthTM"||
			 colname == "numTM" ||
			 colname == "polar" ||
			 colname == "pos" ||
			 colname == "prop1.Tr1221" ||
			 colname == "prop1.Tr1331" ||
			 colname == "prop3.Tr2332" ||
			 colname == "prop5.Tr1331" ||
			 colname == "prop5.Tr2332" ||
			 colname == "prop6.Tr1331" ||
			 colname == "NonTM.W" ||
			 colname == "E" ||
			 colname == "NonTM.pos" ||
			 colname == "Q" ||
			 colname == "TM.E" ||
			 colname == "TM.P" ||
			 colname == "TM.Q" ||
			 colname == "TM.R" ||
			 colname == "TMratio" ||
			 colname == "lengthTM" ||
			 colname == "longestNonTM" ||
			 colname == "neg" ||
			 colname == "TM.W" ||
			 colname == "Xc1.P"
			 
	) {
		newcol <- data.frame(log(col[,colname]))
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "C") {
		col[col[,2] <= 0.005, 2] <- 1
		col[col[,2] > 0.005, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "K") {
		col[col[,2] <= 0.045, 2] <- 1
		col[col[,2] > 0.045, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	} 	else if (colname == "NonTM.E") {
		col[col[,2] <= 0.05, 2] <- 1
		col[col[,2] > 0.05, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "half_life") {
		col[col[,2] < 20, 2] <- 1
		col[col[,2] >= 20, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "NonTM.H" || 
			 colname == "NonTM.I" || 
			 colname == "NonTM.N" ||  
			 colname == "NonTM.Q" ||  
			 colname == "NonTM.S" ||  
			 colname == "neg" ||
			 colname == "TM.T" ||
			 colname == "TM.charge"||
			 colname == "TM.pos"
	) {
		col[col[,2] <= 0.01, 2] <- 1
		col[col[,2] > 0.01, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "NonTM.K"	) {
		col[col[,2] <= 0.08, 2] <- 1
		col[col[,2] > 0.08, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "TM.noCharge"	) {
		col[col[,2] <= 0.98, 2] <- 1
		col[col[,2] > 0.98, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	} 	else if (colname == "TM.sulfur"	) {
		col[col[,2] <= 0.015, 2] <- 1
		col[col[,2] > 0.015, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "avgTM"	) {
		col[col[,2] <= 18, 2] <- 1
		col[col[,2] > 18, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "charge"	) {
		col[col[,2] <= 0.2, 2] <- 1
		col[col[,2] > 0.2, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "fractionTM"	) {
		col[col[,2] <= 0.3, 2] <- 1
		col[col[,2] > 0.3, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "noCharge"	) {
		col[col[,2] <= 0.8, 2] <- 1
		col[col[,2] > 0.8, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "prop7.Tr1221"	) {
		col[col[,2] <= 0.22, 2] <- 1
		col[col[,2] > 0.22, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "NonTM.P") {
		col[col[,2] <= 0.01 || col[,2] >= 0.04, 2] <- 1
		col[col[,2] > 0.01 && col[,2] < 0.04, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "NonTM.T") {
		col[col[,2] <= 0.01 || col[,2] >= 0.05, 2] <- 1
		col[col[,2] > 0.01 && col[,2] < 0.05, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "OB") {
		col[col[,2] <= 2.5 || col[,2] >= 5, 2] <- 1
		col[col[,2] > 2.5 && col[,2] < 5, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "R") {
		col[col[,2] <= 0.025 || col[,2] >= 0.06, 2] <- 1
		col[col[,2] > 0.025 && col[,2] < 0.06, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "TM.N" || colname == "TM.neg") {
		col[col[,2] <= 0.01 || col[,2] >= 0.045, 2] <- 1
		col[col[,2] > 0.01 && col[,2] < 0.045, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "TM.M") {
		col[col[,2] <= 0.01 || col[,2] >= 0.075, 2] <- 1
		col[col[,2] > 0.01 && col[,2] < 0.075, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}	else if (colname == "TM.Y") {
		col[col[,2] <= 0.01 || col[,2] >= 0.08, 2] <- 1
		col[col[,2] > 0.01 && col[,2] < 0.08, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut1 <- cbind(dataOut1, newcol)
	}
}	
rownames(dataOut1) <- rownames(data1)

#2nd step - purification
dataOut2 <- data.frame(rownames(data2))
for (i in 2:ncol(data2)){
	col <- data.frame(data2[,c(1,i)])
	colname <- colnames(col)[2]
	if ( colname == "pI"){
		col$dipI <- NA
		col$dipI[col$pI > 7 ] <- 1
		col$dipI[col$pI <= 7 ] <- 0
		dataOut2 <- cbind(dataOut2 , col$dipI)
	} 	else if ( grepl("lag", colname)  ||
			   grepl("lambda", colname) ||
			   colname == "NonTM.alifatic" ||
			   colname == "NonTM.aromatic" ||
			   colname == "NonTM.charge" ||
			   colname == "NonTM.hydroxil" ||
			   colname == "NonTM.neg" ||
			   colname == "NonTM.noCharge" ||
			   colname == "NonTM.nonpolar" ||
			   colname == "NonTM.polar"
	) {
		
		lowerBound <- dev2[colname,"lowerBound"]
		upperBound <- dev2[colname,"upperBound"]
		
		if(is.na(lowerBound) || is.na(upperBound) ){
			next
		}
		
		col[col[,2] <= lowerBound | col[,2] >= upperBound, 2] <- 1
		col[col[,2] > lowerBound & col[,2] < upperBound, 2] <- 0
		
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	} 	else if (	colname %in% vector2) {
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	}	else if ( colname == "D" || 
			  colname == "E" || 
			  colname == "I" ||
			  colname == "K" ||
			  colname == "L" ||
			  colname == "M" ||
			  colname == "TM.C" ||
			  colname == "TM.D" ||
			  colname == "TM.E" ||
			  colname == "TM.H" ||
			  colname == "TM.I" ||
			  colname == "TM.K" ||
			  colname == "TM.L" ||
			  colname == "TM.N" ||
			  colname == "TM.P" ||
			  colname == "TM.W" ||
			  colname == "TM.aromatic" ||
			  colname == "TM.hydroxil" ||
			  colname == "TM.neg" ||
			  colname == "TM.polar" ||
			  colname == "Xc1.C" ||
			  colname == "Xc1.D" ||
			  colname == "Xc1.E" ||
			  colname == "Xc1.F" ||
			  colname == "Xc1.G" ||
			  colname == "Xc1.H" ||
			  colname == "Xc1.K" ||
			  colname == "Xc1.L" ||
			  colname == "Xc1.R" ||
			  colname == "Xc1.Q" ||
			  colname == "Xc1.S" ||
			  colname == "Xc1.T" ||
			  colname == "Xc1.V" ||
			  colname == "Xc1.W" ||
			  colname == "alifatic" ||
			  colname == "aromatic" ||
			  colname == "buriedratio" ||
			  colname == "length" ||
			  colname == "lengthNonTM" ||
			  colname == "lengthTM" ||
			  colname == "longestNonTM" ||
			  colname == "neg" ||
			  colname == "noCharge" ||
			  colname == "numTM" ||
			  colname == "TMratio"
	) {
		newcol <- data.frame(log(col[,colname]))
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	}	else if (colname == "OB") {
		col[col[,2] <= 2.5 || col[,2] >= 5, 2] <- 1
		col[col[,2] > 2.5 && col[,2] < 5, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	}	else if (colname == "TM.M") {
		col[col[,2] <= 0.02 , 2] <- 1
		col[col[,2] > 0.02 , 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	}	else if (colname == "TM.N" || colname == "TM.R") {
		col[col[,2] <= 0.01 , 2] <- 1
		col[col[,2] > 0.01 , 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	}	else if (colname == "TM.Q" ) {
		col[col[,2] <= 0.005 , 2] <- 1
		col[col[,2] > 0.005 , 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	} 	else if (colname == "TM.noCharge") {
		col[col[,2] <= 0.98 , 2] <- 1
		col[col[,2] > 0.98 , 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	} 	else if (colname == "TM.sulfur") {
		col[col[,2] <= 0.02 , 2] <- 1
		col[col[,2] > 0.02 , 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	}	else if (colname == "Xc1.M") {
		col[col[,2] <= 3 , 2] <- 1
		col[col[,2] > 3 , 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	}	else if (colname == "charge" || colname == "fractionTM") {
		col[col[,2] <= 0.2 , 2] <- 1
		col[col[,2] > 0.2 , 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	}	else if (colname == "gravy") {
		col[col[,2] <= 1.2 , 2] <- 1
		col[col[,2] > 1.2 , 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	} 	else if (colname == "TM.S") {
		col[col[,2] <= 0.02 || col[,2] >= 0.09, 2] <- 1
		col[col[,2] > 0.02 && col[,2] < 0.09, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	}	else if (colname == "TM.T") {
		col[col[,2] <= 0.025 || col[,2] >= 0.08, 2] <- 1
		col[col[,2] > 0.025 && col[,2] < 0.08, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	}	else if (colname == "TM.Y") {
		col[col[,2] <= 0.01 || col[,2] >= 0.08, 2] <- 1
		col[col[,2] > 0.01 && col[,2] < 0.08, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	}	else if (colname == "TM.charge") {
		col[col[,2] <= 0.01 || col[,2] >= 0.09, 2] <- 1
		col[col[,2] > 0.01 && col[,2] < 0.09, 2] <- 0
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut2 <- cbind(dataOut2, newcol)
	}
}	
rownames(dataOut2) <- rownames(data2)

#3rd step - crystallization
dataOut3 <- data.frame(rownames(data3))
for (i in 2:ncol(data3)){
	col <- data.frame(data3[,c(1,i)])
	colname <- colnames(col)[2]
	if (  grepl("lag", colname)  || colname == "Xc2.lambda.20" || colname == "Xc2.lambda.25" || colname == "Xc2.lambda.30") {
		
		lowerBound <- dev3[colname,"lowerBound"]
		upperBound <- dev3[colname,"upperBound"]
		
		if(is.na(lowerBound) || is.na(upperBound) ){
			next
		}
		
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut3 <- cbind(dataOut3, newcol)
	} 	else if (	colname %in% vector3) {
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut3 <- cbind(dataOut3, newcol)
	} 	else if (grepl("prop", colname)) {
		newcol <- data.frame(col[,colname])
		colnames(newcol) <- c(colname)
		dataOut3 <- cbind(dataOut3, newcol)
	}	else if (grepl("Xc1", colname) || 
			 colname == "Y" || 
			 colname == "aromatic" || 
			 colname == "length" ||
			 colname == "lengthNonTM" ||
			 colname == "lengthTM" ||
			 colname == "longestNonTM" ||
			 colname == "TM.neg" ||
			 colname == "numTM") {
		newcol <- data.frame(log(col[,colname]))
		colnames(newcol) <- c(colname)
		dataOut3 <- cbind(dataOut3, newcol)
	}
}
rownames(dataOut3) <- rownames(data3)

pr1 <- predict(model1, dataOut1[,-1], type="prob")
pr2 <- predict(model2, dataOut2[,-1], type="prob")
pr3 <- predict(model3, dataOut3[,-1], type="prob")

prediction <- data.frame(cbind(pr1[1], pr2[1], pr3[1]))
rownames(prediction) <- rownames(data)
colnames(prediction) <- c("pr1", "pr2", "pr3")
prediction$prw <- NA
prediction$prw <- as.numeric(prediction$pr1) + as.numeric(prediction$pr2) + as.numeric(prediction$pr3)

prediction$pred1 <- NA
prediction$pred2 <- NA
prediction$pred3 <- NA
prediction$predw <- NA

prediction$pred1[prediction$pr1>=thr1] <- "success"
prediction$pred1[prediction$pr1<thr1] <- "failure"

prediction$pred2[prediction$pr2>=thr2] <- "success"
prediction$pred2[prediction$pr2<thr2] <- "failure"

prediction$pred3[prediction$pr3>=thr3] <- "success"
prediction$pred3[prediction$pr3<thr3] <- "failure"

prediction$predw[prediction$prw>=thrw] <- "success"
prediction$predw[prediction$prw<thrw] <- "failure"


if ( !is.null(opts$o) ) {
	write.table(format(prediction, digits=4), opts$o, quote = F, sep = "\t", col.names=NA)
}

if (opts$verbose == T || is.null(opts$o)){
	print(format(prediction, digits=4), sep = "\t")
}

