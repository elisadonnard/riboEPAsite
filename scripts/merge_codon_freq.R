library(ggplot2)
library(RColorBrewer)
library(reshape2)
cols=colorRampPalette(brewer.pal(7,"Accent"))(7)
sites = c("codonAsite","codonEsite","codonPsite")
for (i in sites) {
  ribosites = list.files(path = ".", pattern = "codon(E|P|A)site.*tsv")
  files=ribosites[grep(sites[i], ribosites)]
  name = gsub(sites[i],"",files)
  name=gsub("_offsetreads_","",name)
  name=gsub(".tsv","_",name)
  merged=read.table(files[1], header = T)
  colnames(merged)=gsub("^",name[1],colnames(merged))
  for (i in 2:length(files)) {
    fileb=read.table(files[i], header = T)
    colnames(fileb)=gsub("^",name[i],colnames(fileb))
    merged=merge(merged,fileb, by=c(1,2),all=T)
    merged[is.na(merged)]<-0
  }
  colnames(merged)=gsub(".*_codon","codon",colnames(merged))
  colnames(merged)=gsub(".*_aa","aa",colnames(merged))
  ### Select reference samples
  refsamples = scan("reference_samples", what = character())
  freqrefsamples = paste0(refsamples,"_freq")
  getcols = c("codon", "aa", refsamples)
  ### Add average overall freq for reference
  ref=merged[,getcols]
  ref$ref_mean=apply(ref[,3:length(ref)], 1, function(x) mean(x))
  ref$ref_sd=apply(ref[,3:length(ref)], 1, function(x) sd(x))
  ref=ref[,!(colnames(ref)%in%freqrefsamples)]
  merged=merge(ref, merged, by=c(1,2), all=T)
  ### Select reference samples
  refsamples = scan("reference_samples", what = character())
  freqnormrefsamples = paste0(refsamples,"_freqnorm")
  getcols = c("codon", "aa", refsamples)
  ### Add average 3downcodons relative freq for reference
  ref=merged[,getcols]
  ref$ref_3d_mean=apply(ref[,3:length(ref)], 1, function(x) mean(x))
  ref$ref_3d_sd=apply(ref[,3:length(ref)], 1, function(x) sd(x))
  ref=ref[,!(colnames(ref)%in%freqnormrefsamples)]
  merged=merge(ref, merged, by=c(1,2), all=T)
  ###
  melted=melt(merged, variable.name = "key", value.name = "values", 
              id.vars = c("codon","aa","ref_mean","ref_sd","ref_3d_mean","ref_3d_sd"))
  freqs=melted[grep("_freq$", melted$key),]
  freqs$relativeref=freqs$values/freqs$ref_mean
  freqs=freqs[!(freqs$key%in%freqrefsamples),]
  ggplot(freqs, aes(x=interaction(codon, aa), y=relativeref, colour=key)) + 
    geom_point(alpha=0.8) + theme_bw() + 
    theme(axis.text.x=element_text(angle = 270))
  ggsave(filename = paste0("freq_overall_",sites[i],".pdf"), height=6, width=12)
  ###
  freqnorms=melted[grep("_freqnorm$", melted$key),]
  freqnorms$relativeref_3d=freqnorms$values/freqnorms$ref_3d_mean
  freqnorms=freqnorms[!(freqnorms$key%in%freqnormrefsamples),]
  ggplot(freqnorms, aes(x=interaction(codon, aa), y=relativeref_3d, colour=key)) + 
    geom_point(alpha=0.8) + theme_bw() + 
    theme(axis.text.x=element_text(angle = 270))
  ggsave(paste0("freq_3downstream_",sites[i],".pdf"), height=6, width=12)
  ###
  write.table(merged, paste0("merged_",sites[i],".tsv"), sep="\t")
  write.table(freqs, paste0("freq_overall_",sites[i],".tsv"), sep="\t")
  write.table(freqnorms, paste0("freq_3downstream_",sites[i],".tsv"), sep="\t")
}