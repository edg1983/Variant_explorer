#Fix issue whene deNovo is in both affected
vars_df$inh_models[vars_df$inh_models == "deNovo,dominant"] <- "deNovo,deNovo"

#Get counts for recessive, dominant, denovo
inh_counts <- as.data.frame(vars_df %>% separate_rows(inh_models) %>% select(ID,inh_models) %>% group_by(ID) %>% count(inh_models) %>% filter(!is.na(inh_models)))
inh_counts <- spread(inh_counts,key="inh_models", value="n") %>% replace(is.na(.), 0)

#Compute comphet
mydf <- as.data.frame(vars_df %>% select(ID,Gene,inh_from) %>% filter(inh_from %in% c("m","d")) %>% group_by(Gene) %>% filter(sum(inh_from == "d") > 0 & sum(inh_from == "m") > 0))
comphet <- data.frame(d=numeric(),m=numeric(),Gene=character(), stringsAsFactors = F)

for (g in unique(mydf$Gene)) {
  comphet <- rbind(comphet, 
    mydf %>% filter(Gene == g) %>% split(.$inh_from) %>% map(.,1) %>% cross() %>% map_dfr(as_data_frame) %>% mutate(Gene=g) )
}
comphet <- as.data.frame(comphet)
comphet$ID <- paste0("comphet_",1:nrow(comphet))

affected <- c("GT_G179992K","GT_G179994V")
comphet$comphet <- 0

pb <- txtProgressBar(min = 0, max = nrow(comphet), style = 3)
for (n in 1:nrow(comphet)) {
  setTxtProgressBar(pb, n)
  comphet$comphet[n] <- sum(sum(vars_df[vars_df$ID %in% comphet[n,1:2],affected[1]])>=2, sum(vars_df[vars_df$ID %in% comphet[n,1:2],affected[2]])>=2)
}

inh_counts$comphet <- 0
inh_counts <- rbind(inh_counts, comphet %>% select(ID,comphet) %>% mutate(recessive=0,dominant=0,deNovo=0))
write.table(inh_counts, "~/shiny-server/Variant_explorer/example_data/V2.021Jea001.var2reg.segregation.tsv", sep="\t", row.names=F, quote=F)

colnames(comphet) <- c("var1","var2","Gene","ID","comphet")
write.table(comphet[,c(4,1,2,3)], "~/shiny-server/Variant_explorer/example_data/V2.021Jea001.var2reg.comphet.tsv", sep="\t", row.names=F, quote=F)

comphet_df <- as.data.frame(comphet %>% group_by(Gene) %>% mutate(Variants=paste(ID, collapse=","), Variants_n=n()) %>% select(Gene,Variants,Variants_n) %>% distinct() %>% mutate(Inh_model="comphet"))
