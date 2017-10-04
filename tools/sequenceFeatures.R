################################################################################
#										#
#    TMCrys version 1.0								#
#    Copyright 2017 Julia Varga and Gábor E. Tusnády				#
#										#
#    If you use TMCrys, please cite: 						#
#    Julia Varga and Gábor E. Tusnády TMCRys...					#
#										#
#										#
#    This file is part of TMCrys.						#
#										#
#    TMCrys is free software: you can redistribute it and/or modify		#
#    it under the terms of the GNU General Public License as published by	#
#    the Free Software Foundation, either version 3 of the License, or		#
#    (at your option) any later version.					#
#										#
#    TMCrys is distributed in the hope that it will be useful,			#
#    but WITHOUT ANY WARRANTY; without even the implied warranty of		#
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the		#
#    GNU General Public License for more details.				#
#										#
#    You should have received a copy of the GNU General Public License		#
#    along with TMCrys.  If not, see <http://www.gnu.org/licenses/>.		#
#										#
#################################################################################

library("protr", verbose = F, quietly = T, warn.conflicts = F)

args = commandArgs(trailingOnly=TRUE)
seq=args[1]

moreau <- extractMoreauBroto(seq)
moran <- extractMoran(seq)
paac <- extractPAAC(seq)
trans <- extractCTDT(seq)

out <- c(as.list(moreau), as.list(moran), as.list(paac), as.list(trans))
write.table(t(as.data.frame(out)), file='', sep="\t", quote=FALSE)
