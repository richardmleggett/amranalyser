library(ggplot2)
library(reshape2)
library(plyr)
library(scales)
library(params)

#pdf("/Users/leggettr/Desktop/BAMBI/P8(200bp,80pc)_group(70pc).pdf", width=3.1, height=6.5);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/P8(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE, check.names=FALSE)
#highestv <- max(amrdata[][2:4], na.rm=TRUE)

#pdf("/Users/leggettr/Desktop/BAMBI/P8(100bp,70pc)_group(70pc).pdf", width=3.1, height=11.3);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/P8_(100bp70pc)_group(70pc)LH.txt", sep="\t",head=TRUE, check.names=FALSE)
#highestv <- max(amrdata[][2:4], na.rm=TRUE)

#pdf("/Users/leggettr/Desktop/BAMBI/MockWithSpike(200bp,80pc)_group(70pc).pdf", width=4, height=8);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/MockWithSpike(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE)

#pdf("/Users/leggettr/Desktop/BAMBI/BAMBI_EvenMock_nbc_11032019_subsample_lt3000_ss100000(200bp,80pc)_group(70pc).pdf", width=4.5, height=8);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/MockWithSpike(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE)

#pdf("/Users/leggettr/Desktop/BAMBI/BAMBI_EvenMock_nbc_11032019_subsample_lt3000_ss100000_t2(200bp,80pc)_group(70pc).pdf", width=4.5, height=8);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBIMockWithSpikeT2(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE)

#pdf("/Users/leggettr/Desktop/BAMBI/BAMBI_Mock50pc_nbc_11032019_subsample_lt3000_ss100000(200bp,80pc)_group(70pc).pdf", width=4.5, height=8);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/MockWithSpike50pc(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE)

#pdf("/Users/leggettr/Desktop/BAMBI/BAMBI_EvenMock_nbc_11032019_insilico_99999(200bp,80pc)_group(70pc).pdf", width=4.5, height=8);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/MockInsilico(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE)

#pdf("/Users/leggettr/Desktop/BAMBI/BAMBI_EvenMock_nbc_11032019_insilico_88888(200bp,80pc)_group(70pc).pdf", width=3.6, height=6.5);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/MockInsilico8(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE, check.names=FALSE)

pdf("/Users/leggettr/Desktop/BAMBI/P103M_Both(200bp,80pc)_group(70pc).pdf", width=4, height=9.6);
amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/P103M_Both_Spike(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE, check.names=FALSE)
#highestv <- max(amrdata[][2:5], na.rm=TRUE)

#pdf("/Users/leggettr/Desktop/BAMBI/P103M_40pc(200bp,80pc)_group(70pc).pdf", width=4, height=8);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/P103M_40pc(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE)
#highestv <- max(amrdata[][2:5], na.rm=TRUE)

#pdf("/Users/leggettr/Desktop/BAMBI/P103M_4pc(200bp,80pc)_group(70pc).pdf", width=4, height=8);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/P103M_4pc(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE)
#highestv <- max(amrdata[][2:5], na.rm=TRUE)

#pdf("/Users/leggettr/Desktop/BAMBI/P10_all_(200bp,80pc)_group(70pc).pdf", width=9, height=15);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/test_P10_all_(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE, check.names=FALSE)

# #e41a1c
#pdf("/Users/leggettr/Desktop/BAMBI/P10_group2_(200bp,80pc)_group(70pc).pdf", width=4, height=8.7);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/P10_group2_(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE, check.names=FALSE)

# #4daf4a
#pdf("/Users/leggettr/Desktop/BAMBI/P10_group3_(200bp,80pc)_group(70pc).pdf", width=3.6, height=6.6);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/P10_group3_(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE, check.names=FALSE)

# #377eb8
#pdf("/Users/leggettr/Desktop/BAMBI/P10_group1_(200bp,80pc)_group(70pc).pdf", width=4.5, height=2.3);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/P10_group1_(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE, check.names=FALSE)

#pdf("/Users/leggettr/Desktop/BAMBI/Flongle_(200bp,80pc)_group(70pc).pdf", width=5.1, height=6);
#amrdata <- read.delim("/Users/leggettr/Desktop/BAMBI/Flongle_(200bp,80pc)_group(70pc).txt", sep="\t",head=TRUE, check.names=FALSE)


amrdata.m <- melt(amrdata)
#amrdata.m <- ddply(amrdata.m, .(variable), transform, rescale=rescale(value))
p <- ggplot(amrdata.m, aes(variable, Gene)) +
     geom_tile(aes(fill = value), colour = "white") +
     #scale_fill_gradient(name = "Count", low = "#deebf7", high = "#3182bd", na.value = "#dddddd") +
     scale_fill_gradient(name = "Count", low = "#3182bd", high = "#3182bd", na.value = "#dddddd") +
     scale_x_discrete(position = "top") +
     scale_y_discrete(limits = rev(levels(amrdata.m$Gene))) +
     theme(legend.position="none", axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0))

#      scale_fill_gradient(name = "Count", low = "#3182bd", high = "#3182bd", na.value = "#dddddd") +
#     scale_fill_gradient(name = "Count", low = "#3182bd", high = "#3182bd", na.value = "white") +
#scale_fill_gradient(name = "Count", low = "#deebf7", high = "#3182bd", na.value = "white") +
  #theme(legend.position="right", axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle = 90, hjust = 0))
print(p)

garbage <- dev.off()

