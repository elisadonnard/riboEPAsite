library("RColorBrewer")
library("ggplot2")
library("scales")
ribosites = list.files(path = ".", pattern = "codon(E|P|A)site.*bed")
refsites = gsub("codonAsite","codonNsite", ribosites[grep("codonAsite",ribosites)])
for (i in 1:length(ribosites)) {
  sample=read.table(ribosites[i], stringsAsFactors = F)
  ref=read.table(refsites[i], stringsAsFactors = F)
  #
  yy=as.data.frame(table(sample$V8, sample$V9))
  xx=yy[which(yy$Freq!=0),]
  colnames(xx)=c("codon","aa","counts")
  #
  yyref=as.data.frame(table(ref$V8, ref$V9))
  xxref=yyref[which(yyref$Freq!=0),]
  colnames(xxref)=c("codon","aa","countsref")
  #
  merged=merge(xx, xxref, by=c(1,2))
  merged$freq=(merged$counts/sum(merged$counts))*100
  merged$freqref=(merged$countsref/sum(merged$countsref))*100
  merged$freqnorm=(merged$freq/merged$freqref)
  write.table(merged, file = gsub("bed$","tsv",ribosites[i]),sep="\t")
  #
  colourCount = length(unique(merged$aa))
  getPalette = colorRampPalette(brewer.pal(9, "Set1"))
  merged$codon <- factor(merged$codon, levels = merged$codon[order(merged$aa, merged$codon)])
  ggplot_alternative <- function() {ggplot(merged, aes(codon, freqnorm, fill=aa)) + 
      geom_bar(stat = "identity") + scale_fill_manual(values = getPalette(colourCount), name = "Amino acid") + 
      scale_y_continuous() + xlab("Codon") + ylab("normalized frequency") + 
      labs(title = paste0("Codon analysis riboprofiling ",ribosites[i])) + 
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size=7))
  }
  ggsave(filename = gsub("bed$","pdf",ribosites[i]), ggplot_alternative(), width=20, height=6)
}
