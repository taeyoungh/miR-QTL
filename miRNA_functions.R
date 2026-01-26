
# ======================================
# eQTL
# ======================================

plot.eqtl.raw <- function(plot.id, expr.bed, genotype.obj) {
  # generate a boxplot with linear regression line.
  # genotype.obj is an bngsnpr object.
  
  plot.y = subset(expr.bed, gene_id == plot.id$phenotype_id)[-(1:4)]
  plot.x = genotype.obj$genotypes[,which(genotype.obj$map$marker.ID==plot.id$variant_id)]
  
  # make sure that order should be same between x and y.
  idx = match(names(plot.y), genotype.obj$fam$sample.ID)
  plot.x = plot.x[idx]
  plot.cov = genotype.obj$fam[idx,]
  
  plot.df = data.frame(sample=plot.cov$sample.ID,
                       genotype=factor(plot.x, levels=c(0,1,2), labels = c(paste0(plot.id$allele2, plot.id$allele2), paste0(plot.id$allele2, plot.id$allele1), paste0(plot.id$allele1, plot.id$allele1))), 
                       expr=as.numeric(plot.y),
                       sex=as.factor(plot.cov$sex))
  
  p <- ggplot(plot.df, aes(x=genotype, y=expr))
  p <- p + geom_violin()
  p <- p + geom_jitter(aes(col=sex), width=0.1) + geom_smooth(aes(x=as.numeric(genotype), y=expr), method="lm")
  p <- p + ggtitle(paste(plot.id$variant_id, plot.id$phenotype_id, sep=" to "))
  print(p)
  
  return(plot.df)
}

prepare_eqtl_expr_miRNA <- function(rownames.sel, colnames.sel, expr.mat, allow.multi=F) {
  # require global variables: miRcount.id.expanded, miRcount.index
  
  miRcount = expr.mat
  miRcount.rowNames = rownames.sel
  miRcount.colNames = colnames.sel
  
  # clean name : to process a compound name such as "hsa-let-7a-5p/7c-5p". 
  # We need to split it to "hsa-let-7a-5p" and "hsa-let-7c-5p".
  # They can be mapped to different loci.
  
  miRcount.rowNames.expanded = miRcount.id.expanded[[miRcount.rowNames]]
  
  # find genomic loci
  miRcount.index.sel = subset(miRcount.index, name %in% miRcount.rowNames.expanded) %>% distinct(name.cleaned, .keep_all=T) # ignore redundancy due to SNP suffix
  
  #if (any(miRcount.index.sel$loci.num==0)) {cat(paste0(rownames.sel, " are not mapped to genome by name(s).")); return(NULL)}
  
  # allow multiple loci ?
  if (allow.multi) {
    miRcount.index.sel = subset(miRcount.index.sel, loci.num>0)
    if (nrow(miRcount.index.sel)==0) {cat(paste0(miRcount.rowNames, " are not mapped to genome by name(s).\n")); return(NULL)}
  } else {
    miRcount.index.sel = subset(miRcount.index.sel, loci.num==1)
    if (nrow(miRcount.index.sel)==0) {cat(paste0(miRcount.rowNames, " are not mapped or multi mapped to genome by name(s).\n")); return(NULL)}
  }
  
  # generate a bed format by considering multiple loci mapping.
  temp = strsplit(unlist(strsplit(miRcount.index.sel$loci.tss, split=",")), split=":")
  out = data.frame(chrom = gsub("chr", "", sapply(temp, "[[", 1)), 
                   start = as.numeric(sapply(temp, "[[", 2))-1,
                   end = as.numeric(sapply(temp, "[[", 2)),
                   gene_id = paste(unlist(strsplit(miRcount.index.sel$loci.id, split=",")), miRcount.index.sel$name.cleaned, sep="_"))
  
  # append with gene expression values 
  out = cbind(out, miRcount[rep(miRcount.rowNames, length(temp)), miRcount.colNames, drop=F], row.names = NULL)
  
  return(out)
}

run_exprPCA <- function(count.mat) { # Not USED ??
  cat("Assuming that input is count matrix\n")
  
  # Normalization and transformation
  cat("Performing normalization as implemented in DESeq2\n")
  deseq2.obj <- DESeq2::DESeqDataSetFromMatrix(countData = count.mat, colData = data.frame(Sample=colnames(count.mat)), design = ~1)
  deseq2.obj <- DESeq2::estimateSizeFactors(deseq2.obj)
  deseq2.obj <- DESeq2::varianceStabilizingTransformation(deseq2.obj)
  
  # PCA
  temp = SummarizedExperiment::assay(deseq2.obj)
  
  pca.obj <- prcomp(t(), center = T, scale = T)
  return(pca.obj)
}

runMASH <- function(r, input.strong, input.random) { # r is a condition vector of selection.
  # step 1: prepare data format
  cat("step 1: prepare data format\n")
  mash.input.strong = mash_set_data(sapply(r, function(i) { setNames(input.strong[[i]][["slope"]], input.strong[[i]]$phenotype_id)}), sapply(r, function(i) { setNames(input.strong[[i]][["slope_se"]], input.strong[[i]]$phenotype_id)}))
  
  mash.input.random = mash_set_data(sapply(r, function(i) { setNames(input.random[[i]][["slope"]], input.random[[i]]$phenotype_id)}), sapply(r, function(i) { setNames(input.random[[i]][["slope_se"]], input.random[[i]]$phenotype_id)}))
  
  # Accounting for correlations among measurements: In some settings measurements and tests in different conditions may be correlated with one another. For example, in eQTL applications this can occur due to sample overlap among the different conditions.
  v = estimate_null_correlation_simple(mash.input.random) 
  mash.input.random = mash_update_data(mash.input.random, V=v)
  mash.input.strong = mash_update_data(mash.input.strong, V=v)
  
  # step 2: Learn data-driven covariance matrices.
  cat("step 2: learn data-driven covariance matrices.\n")
  U.pca = cov_pca(mash.input.strong, 4) # 4: tissue number
  #U.ed = cov_ed(mash.input.strong, U.pca) # The ED step is helpful primarily because the previous step (PCA) will often estimate the covariances of the observed data (Bhat) whereas what is required is the covariances of the actual underlying effects (B), and this is what ED estimates. 
  
  #U.f = cov_flash(mash.input.strong)
  U.f = cov_flash(mash.input.strong, factors="nonneg", tag="non_neg")
  #U.ed = cov_ed(mash.input.strong, c(U.f, U.pca))
  
  # step 3: Fit the mash model to the random tests, to learn the mixture weights on all the different covariance matrices and scaling coefficients. This step must be performed using all the tests (or a large random subset), because this is where mash learns that many tests are null and corrects for it.
  # In general we recommend running mash with both data-driven and canonical covariances. 
  cat("step 3: fit the mash model to the random data.\n")
  
  U.c = cov_canonical(mash.input.random)
  #m = mash(mash.input.random, Ulist = c(U.ed, U.c), outputlevel = 1) # The outputlevel=1 option means that it will not compute posterior summaries for these tests (which saves time).
  
  # WARNING:: I skipped U.ed because of crash (memory shortage??)
  m = mash(mash.input.random, Ulist = c(U.f, U.c), outputlevel = 1) # The outputlevel=1 option means that it 
  
  # step 4: posterior
  cat("step 4: calculate posterior.\n")
  m2 = mash(mash.input.strong, g=get_fitted_g(m), fixg=TRUE)
  return(m2)
}

get_pairwise_sharing2 = function (m, factor = 0.5, lfsr_thresh = 0.05, FUN = identity, gene_type_sel, gene_annot) 
{ # to make a correlation plot for eQTL for selected gene_types
  print(all(sapply(m$result, function(x) identical(rownames(x), rownames(m$result[[1]])))))
  idx = rownames(m$result$PosteriorMean) %in% subset(temp1, gene_type==gene_type_sel)$gene_id
  m$result = lapply(m$result, function(mm) mm[idx, , drop=FALSE])
  lfsr = get_lfsr(m)
  R = ncol(lfsr)
  S = matrix(NA, nrow = R, ncol = R)
  for (i in 1:R) {
    for (j in i:R) {
      sig_i = get_significant_results(m, thresh = lfsr_thresh, 
                                      conditions = i)
      sig_j = get_significant_results(m, thresh = lfsr_thresh, 
                                      conditions = j)
      a = union(sig_i, sig_j)
      ratio = FUN(get_pm(m)[a, i])/FUN(get_pm(m)[a, j])
      S[i, j] = mean(ratio > factor & ratio < (1/factor))
    }
  }
  S[lower.tri(S, diag = FALSE)] = t(S)[lower.tri(S, diag = FALSE)]
  colnames(S) = row.names(S) = colnames(m$result$PosteriorMean)
  return(S)
}


# ======================================
# Host
# ======================================
library(Rfast) # recommended by coloc package

runColoc <- function(mir.host, baseDir, file.prefix, sample.size) {
  # baseDir and file.prefix for loading a tensorqtl parquet file and a genotype.

  # Read genotypes for variant position and ld calculation
  miRNA.nominal.genotype = bigsnpr::snp_attach(paste0(baseDir, file.prefix, "_miRNA_tsqtl-cis-nominal.rds")) 
  
  mir.host.split = split(mir.host, mir.host$chrom)
  temp.res = list()
  for (i in 1:length(mir.host.split)) { 
    i.res = mir.host.split[[i]]
    i.res$coloc.h4 = NA; i.res$coloc.n=NA;
    i.res$coloc.susie.miRNA = NA; i.res$coloc.susie.mRNA = NA; i.res$coloc.susie.h4 = NA;
    
    i.miRNA.nominal = arrow::read_parquet(paste0(baseDir, "/cis_nominal/", file.prefix, "_miRNA_tsqtl-cis-nominal.cis_qtl_pairs.", gsub("chr", "", i.res$chrom[1]),".parquet"))
    
    i.mRNA.nominal = arrow::read_parquet(paste0(baseDir, "/cis_nominal/", file.prefix, "_mRNA_tsqtl-cis-nominal.cis_qtl_pairs.", gsub("chr", "", i.res$chrom[1]),".parquet"))
    
    for (j in 1:nrow(i.res)) {
      print(paste0(file.prefix, "; ", i, "-th chromosome: ", j))
      coloc.miRNA = subset(i.miRNA.nominal, phenotype_id==i.res$phenotype_id[j])
      coloc.mRNA = subset(i.mRNA.nominal, phenotype_id==i.res$gene_id[j])
      
      coloc.input = left_join(coloc.miRNA, coloc.mRNA, by="variant_id", suffix = c("_miRNA", "_mRNA")) %>% filter(!is.na(slope_mRNA)) # select intersected variants only
      i.res$coloc.n[j] = nrow(coloc.input)
      
      # ld matrix
      coloc.ld = bigsnpr::snp_cor(miRNA.nominal.genotype$genotypes, ind.col= match(coloc.input$variant_id, miRNA.nominal.genotype$map$marker.ID))
      colnames(coloc.ld) = coloc.input$variant_id
      
      # get position
      coloc.input = coloc.input %>% left_join(miRNA.nominal.genotype$map[, c("marker.ID", "chromosome", "physical.pos", "allele1", "allele2")], by= c("variant_id" = "marker.ID")) %>% dplyr::rename(variant_chrom=chromosome, variant_pos=physical.pos)
      
      # prepare coloc inputs
      coloc.data1 = with(coloc.input, list(beta=slope_miRNA, varbeta=slope_se_miRNA^2, snp=variant_id, position=variant_pos, type="quant", sdY=1, N=sample.size, LD=as.matrix(coloc.ld)))
      
      coloc.data2 = with(coloc.input, list(beta=slope_mRNA, varbeta=slope_se_mRNA^2, snp=variant_id, position=variant_pos, type="quant", sdY=1, N=sample.size, LD=as.matrix(coloc.ld)))
      
      #coloc::plot_dataset(coloc.data1)
      #coloc::plot_dataset(coloc.data2)
      
      coloc.res = coloc::coloc.abf(coloc.data1, coloc.data2)
      i.res$coloc.h4[j] = coloc.res$summary["PP.H4.abf"] 
      
      # Susie
      susie = tryCatch({
        coloc.data1.susie = coloc::runsusie(coloc.data1)
        coloc.data2.susie = coloc::runsusie(coloc.data2)
      }, error = function(e) { print(e); return(-1) }) # Error in susie_suff_stat(XtX = XtX, Xty = Xty, n = n, yty = (n - 1) *  : The estimated prior variance is unreasonably large. This is usually caused by mismatch between the summary statistics and the LD matrix. Please check the input.
      # I tested: temp = susieR::susie_rss(bhat = coloc.input$slope_miRNA, shat = coloc.input$slope_se_miRNA, n = sample.size, R = as.matrix(coloc.ld), var_y = 1, L = 10, estimate_residual_variance = TRUE) # Estimating residual variance failed: the estimated value is negative
      # Can this hint of the function be an issue? "when covariates are included in the univariate regressions that produced the summary statistics, also consider removing these effects from X before computing R."...
      
      if (class(susie)=="susie") {
        coloc.susie.res = coloc::coloc.susie(coloc.data1.susie, coloc.data1.susie)
        if (!is.null(coloc.susie.res$summary)) {
          coloc.susie.res$summary = subset(coloc.susie.res$summary, PP.H4.abf>0.1)
          i.res$coloc.susie.miRNA[j] = paste(coloc.susie.res$summary$hit1, collapse=",")
          i.res$coloc.susie.mRNA[j] = paste(coloc.susie.res$summary$hit2, collapse=",")
          i.res$coloc.susie.h4[j] = paste(sprintf("%.6f", coloc.susie.res$summary$PP.H4.abf), collapse=",")
        } else { # SuSiE found no credible sets, meaning there is nothing to colocalize. This can happen if the signal in your data is too weak
          i.res$coloc.susie.miRNA[j] = "no_cs"
          i.res$coloc.susie.mRNA[j] = "no_cs"
          i.res$coloc.susie.h4[j] = "no_cs"
        }
      } 
    }
    temp.res[[i]] = i.res
  }
  out = do.call(rbind, temp.res)
  return(out)
}



# ======================================
# TWAS
# ======================================

prepare_fusion_genotype = function(cur_gene, full_bed) {
  # full_bed: plink bed
  # The global variable tsqtl.nominal is required
  
  cur.locus = subset(tsqtl.nominal, phenotype_id == cur_gene)
  
  # filter by 1000G SNP rsID
  cur.locus = subset(cur.locus, variant_id %in% cur.ldref)
  
  # extract genotypes with plink
  write.table(unique(cur.locus$variant_id), file=paste0(cur_gene, "_variantID.txt"), sep="\n", row.names=F, col.names=F, quote=F)
  
  cmd = paste0("plink --bfile ", full_bed, " --extract ", cur_gene, "_variantID.txt --make-bed --keep-allele-order --out ", cur_gene)
  cat(cmd)
  system(cmd)
}

prepare_fusion_expression = function(cur_gene) {
  # The global variable miRNA.bed is required
  
  cur.expr = subset(miRNA.bed, gene_id==cur_gene)
  
  temp = t(cur.expr[,-(1:4)]) %>% as.data.frame() %>% mutate(expr=.[,1]) %>% tibble::rownames_to_column("IID") %>% mutate(FID=IID) %>% dplyr::select(FID, IID, expr)
  
  write.table(temp, file=paste0(cur_gene, "_expr.txt"), sep="\t", row.names=F, col.names=T, quote=F)
}

# ======================================
# 8/15/2025
# ======================================

# Link to miRBase: check if miRBase.ids are found in miRge.id.
miRgeToMirBase <- function(miRge.id, miRBase.id) {
  
  target = ifelse(miRBase.id %in% miRge.id, miRBase.id, "Not")
  
  # did we find all targets? 
  remained = setdiff(miRge.id, target[target!="Not"])
  print(paste0(length(remained), " are not mapped as they are."))
  if (length(remained)==0) return(target)
  
  # if multiple miRNAs are merged with "/" (e.g. hsa-miR-129-1-3p/129-2-3p), split them.
  remained = sapply(remained, function(s) {gsub("/", paste0(",", substr(s,1,8)),s)}) # first, expand to hsa-miR-129-1-3p,hsa-miR-129-2-3p
  remained = unlist(strsplit(remained, split=",")) # then, split

  # update
  remained2 = c()
  for (i in 1:length(remained)) {
    idx = which(remained[i] == miRBase.id)
    if (length(idx)>0) {
      target[idx] = remained[i]  
    } else {
      remained2 = c(remained2, remained[i])
    }
  }
  print(paste0(length(remained2), " are not mapped after separating merged names"))
  if (length(remained2)==0) return(target)
  
  # assume that the arms are correct and check them.
  remained2 = sapply(strsplit(remained2, split="-"), function(s) paste(s[1:3], collapse="-"))
  
  # update
  remained3 = c()
  for (i in 1:length(remained2)) {
    idx = which(remained2[i] == miRBase.id)
    if (length(idx)>0) {
      target[idx] = remained2[i]  
    } else {
      remained3 = c(remained3, remained2[i])
    }
  }
  
  print(paste0(length(remained3), " are not mapped after removing arms: ", paste(remained3, collapse = ", ")))
  return(target)
}
