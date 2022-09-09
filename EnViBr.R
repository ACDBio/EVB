library(tidyverse)
library(dplyr)
library(ABAData)
library(ABAEnrichment)
library(ReactomePA)
library(enrichplot)
library(org.Hs.eg.db)
library(coldcuts)
library(fgsea)
library(msigdbr)
library(enrichR)
library(clusterProfiler)
#library(plyr)

data("dataset_adult")
#data("dataset_5_stages")
#data("dataset_dev_effect")
#setting working directory to package directory
# tryCatch({
#       setwd(getSrcDirectory()[1])
#     }, error = function(e) {
#       setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
#     })

brainseg_gz_filepath_half<-"./Data/ABA_Human_half/ABA_Human_half.nii.gz"
brainseg_ontology_filepath_half<-"./Data/ABA_Human_half/ABA_human_ontology.csv"
gtexpath<-"./Data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.tsv"
full_brainont_path<-"./Data/ABA_coldcuts_GTEx_ontology.tsv"
celldata_filepath<-'./Data/human_mca_smartseq.tsv'
gtex<-read_tsv(gtexpath)
n_enrichment_categories_toshow<-20
brainontology_data<-read_tsv(full_brainont_path)
load('./Data/atc_targets.rdata')
full_drugdb<-readRDS('./Data/DRUGBANK_DGIDB_CTD_interactions_long.rds')
cellinker<-read_tsv('./Data/human-sMOL_remapped.txt')
aca_ccuts_filepath<-'./Data/dataset_adult_incc.rds'
abaincc<-brainontology_data %>% 
  filter(coldcuts==TRUE) %>% 
  select(structure_id) %>% pull()
hugo<-read_tsv('./Data/HUGO_ensembl_extended.txt')
universe_d<-unique(hugo$`Approved symbol`)
#-----Functions-----
remap_genes<-function(genes){
  df<-tibble(genes=genes)
  df %>% 
    mutate(is_present=1) %>% 
    write_delim('genedf.csv', delim=',')
  system("python3 Remap_DF.py")
  df<-read_csv('mapped_df.csv')
  genes<-df$pipe_genesymbol
  system('rm mapped_df.csv')
  system('rm genedf.csv')
  genes
}


aba_enrich_genes<-function(genes){
  for_aba_enrich<-genes %>% 
    tibble() %>% 
    mutate(is_present=1)
  for_aba_enrich$gene<-for_aba_enrich$.
  for_aba_enrich<-for_aba_enrich %>% 
    dplyr::select(gene, is_present)
  for_aba_enrich<-as.data.frame(for_aba_enrich)
  enrichment_res<-ABAEnrichment::aba_enrich(genes=for_aba_enrich, dataset='adult')
  enrichment_res
}

enrich_pathways<-function(genes, showcats=20){
  genes_entrez<-bitr(genes, fromType="SYMBOL", 
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db) %>% 
    dplyr::select(ENTREZID) %>% 
    pull()
  react_erichments<-enrichPathway(gene=genes_entrez, pvalueCutoff = 1, qvalueCutoff = 1)
  kegg_enrichments<-enrichKEGG(genes_entrez, pvalueCutoff = 1, qvalueCutoff = 1)
  go_enrichments<-enrichGO(genes_entrez, org.Hs.eg.db, pvalueCutoff = 1, qvalueCutoff = 1)
  
  react_enrichments_plot<-dplyr::mutate(react_erichments, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory=showcats)
  
  kegg_enrichments_plot<-dplyr::mutate(kegg_enrichments, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory=showcats)
  
  go_enrichments_plot<-dplyr::mutate(go_enrichments, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory=showcats)
  
  results<-list()
  results$react<-react_erichments
  results$kegg<-kegg_enrichments
  results$go<-go_enrichments
  results$react_plot<-react_enrichments_plot
  results$kegg_plot<-kegg_enrichments_plot
  results$go_plot<-go_enrichments_plot
  results
}



list_to_t2g<-function(gene_set_list, universe=universe_d){
  gs_names<-c()
  gene_symbols<-c()
  for (gs in names(gene_set_list)){
    gene_symbols<-c(gene_symbols, gene_set_list[[gs]])
    gs_names<-c(gs_names, rep(gs, length(gene_set_list[[gs]])))
  }
  res=as.data.frame(tibble(gs_name=gs_names, gene_symbol=gene_symbols))
  #res$gs_name<-gsub('-','_', res$gs_name)
  #res$gs_name<-gsub(',','_', res2g$gs_name)
  
  add_cat<-tibble(gene_symbol=universe_d[!universe_d %in% res$gene_symbol]) %>%  #adding all genes from HUGO which didn't fall into the categories in Undefined category
    mutate(gs_name='Undefined') %>% 
    select(gs_name, gene_symbol)
  
  res<-rbind(res, add_cat)
}



enrich_genes_fgsea<-function(genes, expr_df){ #expr_df - a dataframe with first column 'hgnc_symbol' and other columns speifying expression levels across the analysed objects
  gene_sets<-c()
  gene_sets$gs<-genes
  enr_result<-'Undefined'
  started=0
  for (reg in colnames(expr_df)[c(2:length(colnames(expr_df)))]){
    
    print('Processing structure:')
    print(reg)
    expr_vec<-expr_df[[reg]]
    names(expr_vec)<-expr_df[['hgnc_symbol']]
    fgsea_res<-fgsea(pathways = gene_sets, 
                     stats = expr_vec,
                     scoreType='pos',
                     eps=0)
    
    fgsea_res<-fgsea_res %>% 
      mutate(target=reg) #here there was region=reg
    
    if (started==0){
      enr_result<-fgsea_res
    } else {
      enr_result<-rbind(enr_result, fgsea_res)
    }
    started=1
  }
  enr_result<-enr_result %>% 
    relocate(target) %>% 
    arrange(padj)
}

enrich_genes_geneset<-function(genes, expr_df, quantile){ #expr_df - a dataframe with first column 'hgnc_symbol' and other columns speifying expression levels across the analysed objects, quantile - a quantile to subset top expressed genes for gene set generation
  gene_sets<-c()
  gene_sets$gs<-genes
  enr_result<-'Undefined'
  
  gene_sets<-get_genesets_fromdf(expr_df, quantile, format2g='t2g')
  enrichment=enricher(gene = genes, TERM2GENE = gene_sets)
  results<-list()
  results$allres<-list(enrichment)
  results$df<-enrichment@result
  results$plot<-dplyr::mutate(enrichment, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory=showcats)
  results
}

get_genesets_fromdf<-function(df, quantile, format2g='t2g'){ #generate gene sets from a dataframe by getting a subset of genes with a specified expression quantile threshold
  if (format2g=='t2g'){
    cats=c()
    genes=c()
  } else {
    gene_sets<-list()
  }
  
  for (profile in colnames(df)[c(2:ncol(df))]){
    print('Processing:')
    print(profile)
    expr_df<-df[c('hgnc_symbol',profile)] %>% 
      drop_na()
    
    
    expr_df[[profile]]>quantile(expr_df[[profile]], quantile)
    sel_genes<-expr_df$hgnc_symbol[sel]
    if(format2g=='t2g'){
      genes<-c(genes, sel_genes)
      cats<-c(cats, rep(profile, length(sel_genes)))
    }else{
      gene_sets[[profile]]<-sel_genes
    }
  }
  
  if (format2g == 't2g'){
    gene_sets=as.data.frame(tibble(gs_name=cats, gene_symbol=genes) %>% 
                              drop_na())
  }
  gene_sets
}

enrich_msigdb<-function(genes, showcats=20){
  all_gene_sets = msigdbr(species = "Homo sapiens")
  msigdbr_t2g = all_gene_sets %>% dplyr::distinct(gs_name, gene_symbol) %>% as.data.frame()
  msigdbr_enrichment=enricher(gene = genes, TERM2GENE = msigdbr_t2g)
  print(msigdbr_enrichment)
  results<-list()
  results$allres<-list(msigdbr_enrichment)
  results$df<-msigdbr_enrichment@result
  results$plot<-dplyr::mutate(msigdbr_enrichment, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory=showcats)
  results
}

get_seg_default<-function(brainseg_gz_fp=brainseg_gz_filepath_half, brainseg_ontology_fp=brainseg_ontology_filepath_half){
  seg <- seg_draw(nifti_file = brainseg_gz_fp, 
                  ontology_file = brainseg_ontology_fp,
                  reference_space = "MNI152",
                  verbose = FALSE)
  seg
}


get_var_distribution_brainwise<-function(seg, brainstats, brainseg_ontology_fp=brainseg_ontology_filepath_half, brainseg_gz_fp=brainseg_gz_filepath_half, assay_name='default'){ #Brainstats - a df with variables of interest in rows and indices of brain regions corresponding to the IDs in segmentation ontology in columns. First row in brainstats represents variables.
  varids=brainstats[[1]]
  brainstats=as.data.frame(brainstats[c(2:ncol(brainstats))])
  rownames(brainstats)<-varids
  print(brainstats)
  ont<-read_csv(brainseg_ontology_fp)
  map<-ont$acronym
  names(map)<-ont$id
  map<-as.list(map)
  column_names=colnames(brainstats)
  
  coldata=as.data.frame(tibble(sample_id=column_names, structure_acronym=unlist(map[column_names])))
  rownames(coldata)<-coldata$sample_id
  print(unlist(map[column_names]))
  
  abaAssay <- new("segmentationAssay",
                  values = as.matrix(brainstats),
                  mapping = map,
                  sampledata = coldata)
  seg <- seg_assay_add(segmentation = seg, 
                       assay = abaAssay, 
                       name = assay_name)
  structures_chosen<-intersect(coldata$structure_acronym, seg@ontology$acronym)
  
  print(structures_chosen)
  seg <- seg_projection_add(name = assay_name, 
                            segmentation = seg, 
                            structures = structures_chosen, 
                            planes_chosen = "sagittal")
  seg
}




preprocess_enrichment_res_forcoldcuts<-function(enr_result, brainontology=brainontology_data, brainonsources=c('ABA_dataset_adult','GTEx')){ #Subsets enrichment results according to the specifiedfilters. Outputs a dataframe with columns according to coldcuts ids and rows - to metrics. Brainontology - the file with mappings of expression entity ids to coldcuts ids. It should have columns 'structure_id' (character), 'coldcuts_id' (character), 'Left_hemisphere' (bool),'Right_hemisphere' (Bool),'Undefined_hemisphere' (bool).
  enr_result<-enr_result %>% 
    dplyr::rename(structure_id=target)
  brainontology_dedup<-brainontology %>% 
    distinct(structure_id, .keep_all=TRUE)
  enr_result<-left_join(enr_result, brainontology_dedup)
  enr_result_left<-enr_result %>% 
    drop_na(coldcuts_id) %>% 
    filter(Left_hemisphere==TRUE | Undefined_hemisphere==TRUE)
  enr_result_right<-enr_result %>% 
    drop_na(coldcuts_id) %>% 
    filter(Right_hemisphere==TRUE | Undefined_hemisphere==TRUE)
  enr_result_all<-enr_result%>% 
    drop_na(coldcuts_id)
  print('Processing long to wide...')
  
  
  
  enr_result_left_minp<-long_to_wide_enrichment_df_forcoldcuts(df=enr_result_left, add_suffix = 'Left')
  enr_result_right_minp<-long_to_wide_enrichment_df_forcoldcuts(df=enr_result_right, add_suffix = 'Right')
  enr_result_all_minp<-long_to_wide_enrichment_df_forcoldcuts(df=enr_result_all, add_suffix = 'All')
  
  enr_result_left_med<-long_to_wide_enrichment_df_forcoldcuts(df=enr_result_left, add_suffix = 'Left', agg='median')
  enr_result_right_med<-long_to_wide_enrichment_df_forcoldcuts(df=enr_result_right, add_suffix = 'Right', agg='median')
  enr_result_all_med<-long_to_wide_enrichment_df_forcoldcuts(df=enr_result_all, add_suffix = 'All', agg='median')
  
  enr_result_left_mean<-long_to_wide_enrichment_df_forcoldcuts(df=enr_result_left, add_suffix = 'Left', agg='mean')
  enr_result_right_mean<-long_to_wide_enrichment_df_forcoldcuts(df=enr_result_right, add_suffix = 'Right', agg='mean')
  enr_result_all_mean<-long_to_wide_enrichment_df_forcoldcuts(df=enr_result_all, add_suffix = 'All', agg='mean')
  
  results<-c()
  results$left_hemi_minp<-enr_result_left_minp
  results$right_hemi_minp<-enr_result_right_minp
  results$all_minp<-enr_result_all_minp
  
  results$left_hemi_med<-enr_result_left_med
  results$right_hemi_med<-enr_result_right_med
  results$all_med<-enr_result_all_med
  
  results$left_hemi_mean<-enr_result_left_mean
  results$right_hemi_mean<-enr_result_right_mean
  results$all_mean<-enr_result_all_mean
  
  results$full_annotated_enrichment_data<-enr_result
  print('done.')
  results
}


long_to_wide_enrichment_df_forcoldcuts<-function(df, add_suffix='N', agg='minpadj'){
  if (agg=='minpadj'){
    df<-df %>%  #leaving only most enriched regions within each coldcuts id
      select(coldcuts_id, pval, padj, ES, NES) %>% 
      #mutate(neg_logpadj=-log10(padj)) %>% 
      group_by(coldcuts_id) %>% 
      filter(padj==min(padj)) %>% 
      group_by(coldcuts_id) %>% 
      summarise_all(.funs='max') %>% 
      distinct()
  } else if (agg=='median'){
    df<-df %>%  #leaving only most enriched regions within each coldcuts id
      select(coldcuts_id, pval, padj, ES, NES) %>% 
      #mutate(neg_logpadj=-log10(padj)) %>% 
      group_by(coldcuts_id) %>% 
      summarise_all(.funs='median')
  } else if (agg=='mean'){
    df<-df %>%  #leaving only most enriched regions within each coldcuts id
      select(coldcuts_id, pval, padj, ES, NES) %>% 
      #mutate(neg_logpadj=-log10(padj)) %>% 
      group_by(coldcuts_id) %>% 
      summarise_all(.funs='mean')
  }
  print(df)
  df<-df %>%
    pivot_longer(cols = -coldcuts_id, names_to = 'attribute') %>% 
    pivot_wider(names_from = coldcuts_id, values_from = value)
  #df$attribute<-paste(df$attribute, add_suffix, sep='_')
  df
}

preprocess_abaenrichment_output_forcoldcuts<-function(enr_result, brainontology=brainontology_data){
  
  brainontology_dedup<-brainontology %>% 
    distinct(structure_id, .keep_all=TRUE)
  
  enr_result<-left_join(enr_result, brainontology_dedup)
  enr_result$neg_logp_minfwer<-unlist(-log10(enr_result$min_FWER))
  enr_result$neg_logp_meanfwer<-unlist(-log10(enr_result$mean_FWER))
  
  print(enr_result)
  print(colnames(enr_result))
  enr_result_left<-enr_result %>% 
    drop_na(coldcuts_id) %>% 
    filter(Left_hemisphere==TRUE | Undefined_hemisphere==TRUE)
  enr_result_right<-enr_result %>% 
    drop_na(coldcuts_id) %>% 
    filter(Right_hemisphere==TRUE | Undefined_hemisphere==TRUE)
  enr_result_all<-enr_result%>% 
    drop_na(coldcuts_id)
  print('Processing long to wide...')  

  
  enr_result_left_minp<-long_to_wide_enrichment_df_forcoldcuts_abaenrichment(df=enr_result_left, add_suffix = 'minpadj', agg='minpadj')
  enr_result_right_minp<-long_to_wide_enrichment_df_forcoldcuts_abaenrichment(df=enr_result_right, add_suffix = 'minpadj', agg='minpadj')
  enr_result_all_minp<-long_to_wide_enrichment_df_forcoldcuts_abaenrichment(df=enr_result_all, add_suffix = 'minpadj', agg='minpadj')
  
  enr_result_left_med<-long_to_wide_enrichment_df_forcoldcuts_abaenrichment(df=enr_result_left, add_suffix = 'median', agg='median')
  enr_result_right_med<-long_to_wide_enrichment_df_forcoldcuts_abaenrichment(df=enr_result_right, add_suffix = 'median', agg='median')
  enr_result_all_med<-long_to_wide_enrichment_df_forcoldcuts_abaenrichment(df=enr_result_all, add_suffix = 'median', agg='median')
  
  enr_result_left_mean<-long_to_wide_enrichment_df_forcoldcuts_abaenrichment(df=enr_result_left, add_suffix = 'mean', agg='mean')
  enr_result_right_mean<-long_to_wide_enrichment_df_forcoldcuts_abaenrichment(df=enr_result_right, add_suffix = 'mean', agg='mean')
  enr_result_all_mean<-long_to_wide_enrichment_df_forcoldcuts_abaenrichment(df=enr_result_all, add_suffix = 'mean', agg='mean')
  
  results<-c()
  
  
  results$left_hemi<-plyr::rbind.fill(enr_result_left_minp, enr_result_left_med, enr_result_left_mean)
  results$right_hemi<-plyr::rbind.fill(enr_result_right_minp, enr_result_right_med, enr_result_right_mean)
  results$all<-plyr::rbind.fill(enr_result_all_minp, enr_result_all_med, enr_result_all_mean)
  
  results$full_annotated_enrichment_data<-enr_result
  print('done.')
  results  
}


long_to_wide_enrichment_df_forcoldcuts_abaenrichment<-function(df, add_suffix='N', agg='minpadj'){
  print(add_suffix)
  print(agg)
  if (agg=='minpadj'){
    df<-df %>%  #leaving only most enriched regions within each coldcuts id
      select(coldcuts_id, n_significant, mean_FWER, min_FWER, neg_logp_minfwer, neg_logp_meanfwer) %>% 
      group_by(coldcuts_id) %>% 
      filter(min_FWER==min(min_FWER)) %>% 
      group_by(coldcuts_id) %>% 
      summarise_all(.funs='max') %>% 
      distinct()
    print(df)  
    
  } else if (agg=='median'){
    df<-df %>%  #leaving only most enriched regions within each coldcuts id
      select(coldcuts_id, n_significant, mean_FWER, min_FWER, neg_logp_minfwer, neg_logp_meanfwer) %>% 
      group_by(coldcuts_id) %>% 
      summarise_all(.funs='median')
    print(df)  
  } else if (agg=='mean'){
    df<-df %>%  #leaving only most enriched regions within each coldcuts id
      select(coldcuts_id, n_significant, mean_FWER, min_FWER, neg_logp_minfwer, neg_logp_meanfwer) %>% 
      group_by(coldcuts_id) %>% 
      summarise_all(.funs='mean')
    print(df)  
  }
  
  
  df<-df %>%
    pivot_longer(cols = -coldcuts_id, names_to = 'attribute') %>% 
    pivot_wider(names_from = coldcuts_id, values_from = value)
  df$attribute<-paste(df$attribute, add_suffix, sep='_')
  df
}




subset_drugs_bygenes<-function(genes, full_drugdbd=full_drugdb){
  results<-c()
  results$df<-full_drugdbd %>% 
    dplyr::filter(pipe_genesymbol %in% genes)
  results$drugs<-unique(results$df$DRUG_Common_name)
  results
}



subset_metabolites_bygenes<-function(genes, showcats, cellinkerd=cellinker){
  results<-c()
  results$df<-cellinkerd %>% 
    filter(pipe_genesymbol %in% genes)
  results$small_mol_ligands<-unique(results$df$`ligand name`)
  t2g<-cellinkerd %>% 
    dplyr::select(`ligand name`,pipe_genesymbol) %>% 
    dplyr::rename(gs_name=`ligand name`, gene_symbol=pipe_genesymbol)
  
  tryCatch({
    results$enrichment_allres<-enricher(genes, TERM2GENE=t2g, pvalueCutoff = 1, minGSSize=1, qvalueCutoff=1)
    results$enrichment_df<-results$enrichment_allres@result
    results$plot<-dplyr::mutate(results$enrichment_allres, qscore = -log(p.adjust, base=10)) %>% 
      barplot(x="qscore", showCategory=showcats)
  }, error = function(e) {
    print('Lack of genetic overlap!')
  }
  )
  results
}



enrich_atc<-function(genes, showcats=20){
  atc_targets_t2g<-list_to_t2g(atc_targets)
  results<-c()
  
  tryCatch({
  res<-enricher(genes, TERM2GENE=atc_targets_t2g, pvalueCutoff = 1, minGSSize=1, qvalueCutoff=1)
  results$plot<-dplyr::mutate(res, qscore = -log(p.adjust, base=10)) %>% 
    barplot(x="qscore", showCategory=showcats)
  results$allres<-list(res)
  results$df<-res@result
  }, error = function(e) {
    print('Lack of genetic overlap!')
  }
  )
  results
}

get_expression_df<-function(sources=c('gtex','aba_ds_adult','aba_cells'), filter=FALSE,gtex_fp=gtexpath, celldata_fp=celldata_filepath, aba_ccts_rds_fp=aca_ccuts_filepath){
  resdf<-c()
  init=0
  for (source in sources){
    if (source=='gtex'){
      data<-read_tsv(gtex_fp)%>% 
        dplyr::select(-Name) %>% 
        dplyr::rename(hgnc_symbol=Description)
    } else if (source=='aba_ds_adult' && filter==FALSE){
      data<-dataset_adult %>% dplyr::select(hgnc_symbol, structure, signal) %>%
        group_by(hgnc_symbol, structure) %>%
        summarise_at('signal', .funs='mean') %>% 
        spread(key='structure', value='signal')
      
    } else if (source=='aba_cells'){
      data<-read_tsv(celldata_fp)
    } else if (source=='aba_ds_adult' && filter=='coldcuts'){
      data<-readRDS(aba_ccts_rds_fp)
      
    }
    
    if (init==0){
      resdf<-data
    } else {
      resdf<-left_join(resdf, data)
    }
   init=1
  }

  resdf<-resdf %>% 
    drop_na()
}
#----TESTING COMPONENTS----
# genes<-read_tsv('./Input_examples/gene_list.csv')$genes
# atc_targets_t2g<-list_to_t2g(atc_targets)
# hugo<-read_tsv('./Data/HUGO_ensembl_extended.txt')
# universe_d<-unique(hugo$`Approved symbol`)
# 
# atc_targets_t2g
# atc_targets_t2g$gs_name<-gsub('-','_', atc_targets_t2g$gs_name)
# atc_targets_t2g$gs_name<-gsub(',','_', atc_targets_t2g$gs_name)
# 
# 
# length(unique(atc_targets_t2g$gene_symbol))
# length(unique(universe_d))
# 
# 
# 
# atc_targets_t2g
# 
# add_cat<-tibble(gene_symbol=universe_d[!universe_d %in% atc_targets_t2g$gene_symbol]) %>% 
#   mutate(gs_name='Undefined') %>% 
#   select(gs_name, gene_symbol)
# 
# atc_targets_t2g<-rbind(atc_targets_t2g, add_cat)
# 
# res<-enricher(genes, TERM2GENE=atc_targets_t2g, pvalueCutoff = 10, minGSSize=0, qvalueCutoff=10, pAdjustMethod='none')
# res@result
# 
# 
# 
# 
# intersect(genes, unique(atc_targets_t2g$gene_symbol))
# res@geneSets
# 
# dplyr::mutate(res, qscore = -log(p.adjust, base=10)) %>% 
#   barplot(x="qscore", showCategory=20)

#----TESTING-END----
#----DEBUG----
# genes<-read_tsv('/home/biorp/Gitrepos/Psychiatry/ANHEDONIA/EnViBr/Input_examples/gene_list.csv')$genes
# expr_df<-get_expression_df(sources=c('gtex','aba_ds_adult'))
# enr_results<-enrich_genes_fgsea(genes, expr_df)
# print('Preprocessing enrichment results...')
# print(enr_results)
# enr_results<-preprocess_enrichment_res_forcoldcuts(enr_results, brainonsources=c('ABA_dataset_adult','GTEx'))
# view(enr_results)
# seg<-get_seg_default()
# print('Getting the assay for the left hemi...')
# enr_results$left_hemi
# 
# seg<-get_var_distribution_brainwise(seg=seg, brainstats=enr_results$left_hemi, assay_name='left')
# 
# 
# print('Getting the assay for the right hemi...')
# seg<-get_var_distribution_brainwise(seg=seg, brainstats=enr_results$right_hemi, assay_name='right')
# left_nes_plot<-seg_feature_plot(segmentation = seg,
#                                 assay = "left",
#                                 feature = "NES_Left",
#                                 smooth=FALSE)
# right_nes_plot<-seg_feature_plot(segmentation = seg,
#                                  assay = "right",
#                                  feature = "NES_Right",
#                                  smooth=FALSE)
#----DEBUG-end-----

#-----NEW FEATURE-----
# 
# genes<-read_tsv('/home/biorp/Gitrepos/Psychiatry/ANHEDONIA/EnViBr/Input_examples/gene_list.csv')$genes
# genes_vec<-rep(1, length(genes))
# names(genes_vec)<-genes
# abares<-aba_enrich(genes=genes_vec)
# # 
# abares<-abares$results
# abares$structure_id<-word(abares$structure_id, 2, sep=':')
# abares
# abares<-abares %>% 
#    filter(structure_id %in% abaincc)
# # 
# brainontology_dedup<-brainontology_data %>% 
#   distinct(structure_id, .keep_all=TRUE)
# enr_result<-left_join(abares, brainontology_dedup)
# enr_result$neg_logp_minfwer<-unlist(-log10(enr_result$min_FWER))
# enr_result$neg_logp_meanfwer<-unlist(-log10(enr_result$mean_FWER))
# # e
# # 
# # 
# # 
# view(abares)
# enr_results<-preprocess_abaenrichment_output_forcoldcuts(enr_result = abares)
# enr_results$full_annotated_enrichment_data
# cor.test(enr_results$full_annotated_enrichment_data$n_significant, enr_results$full_annotated_enrichment_data$neg_logp_minfwer)
# enr_results$right_hemi
# 
# 
# enr_results$full_annotated_enrichment_data
# seg_a<-get_seg_default()
# seg_a<-get_var_distribution_brainwise(seg=seg_a, brainstats=enr_results$left_hemi, assay_name='left')
# seg_a<-get_var_distribution_brainwise(seg=seg_a, brainstats=enr_results$right_hemi, assay_name='right')
# seg_a<-get_var_distribution_brainwise(seg=seg_a, brainstats=enr_results$all, assay_name='all')
# # 
# seg_feature_plot(segmentation = seg_a,
#                                   assay = "all",
#                                   feature = "n_significant_minpadj",
#                                   smooth=FALSE,
#                                   labelsize=6,
#                                   remove_axes=TRUE)
# 
# seg_feature_plot(segmentation = seg_a,
#                                    assay = "right",
#                                    feature = "neg_logp_minfwer_minpadj",
#                                    smooth=FALSE,
#                                    labelsize=6,
#                                    remove_axes=TRUE)
# seg_a@assays$left@values
# 
# 
# view(enr_results$full_annotated_enrichment_data)
# left_nsign_plot
# enr_results$left_hemi
#----NEW FEATURE-end----
# enr_results_a<-enr_results
# seg_a<-get_var_distribution_brainwise(seg=seg_a, brainstats=enr_results_a$left_hemi, assay_name='left')
# print('Getting the assay for the right hemi...')
# seg_a<-get_var_distribution_brainwise(seg=seg_a, brainstats=enr_results_a$right_hemi, assay_name='righ')
# print('Getting the assay for both hemi...')
# seg_a<-get_var_distribution_brainwise(seg=seg_a, brainstats=enr_results_a$all, assay_name='all')
# 
# nsign<-enr_results_a$full_annotated_enrichment_data %>% 
#   filter(coldcuts==TRUE) %>% 
#   select(n_significant) %>% 
#   pull()
# min_ns<-min(nsign)
# max_ns<-max(nsign)
# ns_span<-c(min_ns, max_ns)
# 
# 
# minfwer<-enr_results_a$full_annotated_enrichment_data %>% 
#   filter(coldcuts==TRUE) %>% 
#   select(neg_logp_minfwer) %>% 
#   pull()
# min_minfwer<-min(minfwer)
# max_minfwer<-max(minfwer)
# neg_logp_minfwer_span<-c(min_minfwer, max_minfwer)
# 
# 
# 
# meanfwer<-enr_results_a$full_annotated_enrichment_data %>% 
#   filter(coldcuts==TRUE) %>% 
#   select(neg_logp_meanfwer) %>% 
#   pull()
# min_meanfwer<-min(meanfwer)
# max_meanfwer<-max(meanfwer)
# neg_logp_meanfwer_span<-c(min_meanfwer, max_meanfwer)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# seg_feature_plot(segmentation = seg_a,
#                                   assay = "left",
#                                   feature = "n_significant_minpadj",
#                                   smooth=FALSE,
#                                   labelsize=6,
#                                   remove_axes=TRUE,
#                                   rng=ns_span)
# seg_feature_plot(segmentation = seg_a,
#                                    assay = "right",
#                                    feature = "n_significant_minpadj",
#                                    smooth=FALSE,
#                                    labelsize=6,
#                                    remove_axes=TRUE,
#                                    rng=ns_span)
# seg_feature_plot(segmentation = seg_a,
#                                  assay = "all",
#                                  feature = "n_significant_minpadj",
#                                  smooth=FALSE,
#                                  labelsize=6,
#                                  remove_axes=TRUE,
#                                  rng=ns_span)
# 
# 
# seg_feature_plot(segmentation = seg_a,
#                                     assay = "left",
#                                     feature = "neg_logp_minfwer_minpadj",
#                                     smooth=FALSE,
#                                     labelsize=6,
#                                     remove_axes=TRUE,
#                                     rng=neg_logp_minfwer_span)
# seg_feature_plot(segmentation = seg_a,
#                                      assay = "right",
#                                      feature = "neg_logp_minfwer_minpadj",
#                                      smooth=FALSE,
#                                      labelsize=6,
#                                      remove_axes=TRUE,
#                                      rng=neg_logp_minfwer_span)
# seg_feature_plot(segmentation = seg_a,
#                                    assay = "all",
#                                    feature = "neg_logp_minfwer_minpadj",
#                                    smooth=FALSE,
#                                    labelsize=6,
#                                    remove_axes=TRUE,
#                                    rng=neg_logp_minfwer_span)  
# 
# 
# seg_feature_plot(segmentation = seg_a,
#                                      assay = "left",
#                                      feature = "neg_logp_meanfwer_mean",
#                                      smooth=FALSE,
#                                      labelsize=6,
#                                      remove_axes=TRUE,
#                                      rng=neg_logp_minfwer_span)
# seg_feature_plot(segmentation = seg_a,
#                                       assay = "right",
#                                       feature = "neg_logp_meanfwer_mean",
#                                       smooth=FALSE,
#                                       labelsize=6,
#                                       remove_axes=TRUE,
#                                       rng=neg_logp_minfwer_span)
# seg_feature_plot(segmentation = seg_a,
#                                     assay = "all",
#                                     feature = "neg_logp_meanfwer_mean",
#                                     smooth=FALSE,
#                                     labelsize=6,
#                                     remove_axes=TRUE,
#                                     rng=neg_logp_minfwer_span)   
#Error in rownames: trying to get slot "values" from an object of a basic class ("NULL") with no slots