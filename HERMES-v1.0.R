#######################################################################################################################
#                                                                                                                     #
# HERMES is a straightforward index which tries to summarize the mitochondrial evolution pace using a single number.  #
# Several mitogenomic features are evaluated in a factor analysis framework; namely, in the current version, they are #
# %URs, AMIGA, SU skew, root-to-tip distance and pairwise ML distance.                                                #
#                                                                                                                     #
# Copyright (C) 2016 Guglielmo Puccio, Federico Plazzi                                                                #
#                                                                                                                     #
# This program is free software: you can redistribute it and/or modify                                                #
# it under the terms of the GNU General Public License as published by                                                #
# the Free Software Foundation, either version 3 of the License, or                                                   #
# (at your option) any later version.                                                                                 #
#                                                                                                                     #
# This program is distributed in the hope that it will be useful,                                                     #
# but WITHOUT ANY WARRANTY; without even the implied warranty of                                                      #
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                                       #
# GNU General Public License for more details.                                                                        #
#                                                                                                                     #
# You should have received a copy of the GNU General Public License                                                   #
# along with this program.  If not, see <http://www.gnu.org/licenses/>.                                               #
#                                                                                                                     #
#######################################################################################################################

#HERMES version: 1.0

#Load library psych for factor analysis.

library(psych)

#Reading the command line searching for the alpha argument.

HERMES.args <- commandArgs(trailingOnly=TRUE)
if (length(HERMES.args) == 0) {
	alpha <- 0.05
	} else alpha <- as.numeric(HERMES.args[1])

#Read input file and parsing shades (note that the maximum allowed number of groups/shades is 657 - see the colors/colours function).

HERMES.file <- read.table(file="./Results/HERMES_variables.txt",header=TRUE,sep="\t",row.names=1)
if ("Shades" %in% colnames(HERMES.file)) {
	HERMES.shades <- HERMES.file[,"Shades"]
	HERMES.file <- HERMES.file[colnames(HERMES.file)[colnames(HERMES.file) != "Shades"]]
	} else HERMES.shades <- rep("black",length(rownames(HERMES.file)))
HERMES.communalities <- paste(colnames(HERMES.file)," communality",sep="")

#Factor analysis.

HERMES.fa <- fa(HERMES.file,rotate="varimax",scores="tenBerge",fm="ml",cor="cor",normalize=TRUE,missing=TRUE,impute="median",alpha=alpha)
HERMES.fa.scores <- HERMES.fa$scores
plot(HERMES.fa.scores,main="Hyper-Empirical Relative Mitochondrial Evolutionary Speed",xlab="Species",ylab="HERMES",col=colors()[HERMES.shades],bg=colors()[HERMES.shades],pch=21,xaxt="n")
text(x=c(1:length(rownames(HERMES.file))),y=HERMES.fa.scores,labels=c(1:length(rownames(HERMES.file))),pos=2,col=colors()[HERMES.shades])

#Print parameter table.

HERMES.KMO <- KMO(HERMES.file)$MSA
HERMES.TLI <- HERMES.fa$TLI
HERMES.SRMR <- HERMES.fa$rms
HERMES.RMSEA <- as.numeric(HERMES.fa$RMSEA[1])
HERMES.RMSEA.lower <- as.numeric(HERMES.fa$RMSEA[2])
HERMES.RMSEA.upper <- as.numeric(HERMES.fa$RMSEA[3])
HERMES.total.communality <- mean(as.numeric(HERMES.fa$communality))
parameters <- c("KMO","TLI","SRMR","RMSEA",paste("RMSEA ",(1-alpha)*100,"% CI (lower bound)",sep=""),paste("RMSEA ",(1-alpha)*100,"% CI (upper bound)",sep=""),HERMES.communalities,"Total communality")
values <- c(HERMES.KMO,HERMES.TLI,HERMES.SRMR,HERMES.RMSEA,HERMES.RMSEA.lower,HERMES.RMSEA.upper,as.numeric(HERMES.fa$communality),HERMES.total.communality)
HERMES.parameters <- data.frame(parameters=parameters,values=values)
write.table(HERMES.parameters,file="./Results/HERMES_parameters.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=FALSE)
HERMES.fa.scores.data.frame <- data.frame(Species=rownames(HERMES.file),ID=c(1:length(rownames(HERMES.file))),HERMES=HERMES.fa.scores)
write.table(HERMES.fa.scores.data.frame,file="./Results/HERMES_scores.txt",quote=FALSE,sep="\t",row.names=FALSE,col.names=c("Species","#","HERMES"))
