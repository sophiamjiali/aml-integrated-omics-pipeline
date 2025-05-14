# ------------------------------------------------------------------------------
# Script:       generate_id_mapping.R
# Purpose:      Generates mappings for methylation/expression data and cohorts
# Author:       Sophia Li
# Affiliation:  Sadreyev lab
# Date:         2025-05-12
#
# Version:      R version 4.4.2 (2024-10-31)
# ------------------------------------------------------------------------------

# ===| A. Generate Probe Annotations |==========================================

generate_probe_annotation <- function(methylation_cohorts, config) {

  assay <- config$methylation$assay
  
  # --- Fetch annotations for the correct assay ---
  if (assay == "Illumina EPIC") {
    return (generate_epic_annotation(methylation_cohorts, config))
    
  } else if (assay == "...") {
    # ... add further support here
    
  } else {
    stop ("The dataset references an unsupported assay: ", assay)
  }
}


# --- Sub-functions to generate annotations for the supported assays ---

generate_epic_annotation <- function(methylation_cohorts, config) {
  
  message("Building CpG probe annotations for the Illumina EPIC assay...")
  
  # --- Build probe annotations from the official EPIC annotation ---
  annotation <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)
  
  cpg_anno <- GRanges(
    seqnames = annotation$chr,
    ranges = IRanges(start = annotation$pos, width = 1),
    strand = annotation$strand,
    probeID = rownames(annotation)
  )

  # --- Standardize genomic features and seqinfo ---
  seqinfo(cpg_anno) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(cpg_anno)]
  genome(cpg_anno) <- "hg19"
  seqlevelsStyle(cpg_anno) <- "UCSC"
  
  # --- Filter annotations for probes present in the methylation data
  cpg_anno <- cpg_anno[cpg_anno$probeID %in% rownames(methylation_cohorts$grp1), ]
  
  # --- Filter annotation based on dataset configurations ---
  
  # A. Remove non-CpG sites based on configuration settings ---
  preproc <- config$methylation$preproc
  
  if (preproc$remove_non_cpg == TRUE) {
    message("Non-CpG sites successfully removed")
    cpg_anno <- cpg_anno[grepl("^cg", cpg_anno$probeID), ]
  }
  
  # B. Remove probes associated with sex chromosomes
  if (preproc$remove_sex_chr == TRUE) {
    message("Probes associated with sex chromosomes successfully removed")
    cpg_anno <- cpg_anno[!seqnames(cpg_anno) %in% c("chrX", "chrY")]
    cpg_anno <- dropSeqlevels(cpg_anno, setdiff(seqlevels(cpg_anno), 
                              seqlevelsInUse(cpg_anno)), pruning.mode = "coarse")
  }
  
  # C. Remove cross-reactive probes
  if (preproc$remove_cr_probes == TRUE) {
    message("Cross-reactive probes successfully removed")
    cross_reactive <- read.csv(file.path(
      paths$data_default, 
      "epic_cross_reactive_probes.csv")
    )
    cross_reactive_ids <- cross_reactive$PROBE
    cpg_anno <- cpg_anno[!cpg_anno$probeID %in% cross_reactive_ids]
  }
  
  # D. Remove SNP-overlapped probes
  if (preproc$remove_snp_probes == TRUE) {
    message("SNP-overlapped probes successfully removed")
    
    main_anno <- annotation[annotation$Name %in% cpg_anno$probeID, ]
    
    # remove probes with SNPs at CpG, SBE, or in probe body (standard QC)
    has_snp <- function(x) !is.na(x) & x != ""
    to_remove <- with(main_anno, has_snp(Probe_rs) | has_snp(CpG_rs) | has_snp(SBE_rs))
    snp_probes <- main_anno$Name[to_remove]
    
    cpg_anno <- cpg_anno[!cpg_anno$probeID %in% snp_probes]
  }
  
  # --- Save and return the processed CpG probe annotation ---
  processed_path <- file.path(
    paths$data_processed, 
    paste0(config$general$dataset_name, "_cpg_annotation.rds")
  )
  saveRDS(cpg_anno, processed_path)
  message("Wrote CpG probe annotation to: ", processed_path)
  
  return(cpg_anno)
}

# ===| B. Generate Promoter Annotations |=======================================

# --- Generate protein-coding gene promoter annotations ---
generate_promoters <- function(config) {
  
  genome <- config$general$genome_build
  
  # --- Fetch ensembl annotations for the defined genome build ---
  if (genome == "hg19") {
    
    message("Generating promoter annotations for genome build hg19...")
    
    # fetch ensembl protein-coding gene annotations for hg19
    ensdb <- EnsDb.Hsapiens.v75
    genes <- genes(ensdb, filter = ~ gene_biotype == "protein_coding")
    
    # define promoter regions as ±3kB TSS
    suppressWarnings({
      # suppress 'trim' warning, cannot be removed
      promoters <- promoters(genes, upstream = 3000, downstream = 3000)
    })
    
    # restore and standardize metadata
    seqinfo(promoters) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(promoters)]
    seqlevelsStyle(promoters) <- "UCSC"
    promoters <- keepStandardChromosomes(promoters, pruning.mode = "coarse")
    
    # remove sex chromosomes if defined in configurations
    if (config$methylation$preproc$remove_sex_chr == TRUE) {
      
      message("Promoters associated with sex chromosomes successfully removed")
      
      promoters <- promoters[!seqnames(promoters) %in% c("chrX", "chrY", "chrM")]
      promoters <- dropSeqlevels(
        promoters, 
        setdiff(seqlevels(promoters), seqlevelsInUse(promoters)), 
        pruning.mode = "coarse"
      )
      
      # restore seqinfo post-removal
      seqinfo(promoters) <- seqinfo(BSgenome.Hsapiens.UCSC.hg19)[seqlevels(promoters)]
    }
    
    message("Successfully generated promoter annotations")
    
  } else if (genome == "hg38") {
    
     # ... add further support here
    
  } else {
    stop ("The dataset references an unsupported genome build: ", genome)
  }
  
  # --- Save and return the processed CpG probe annotation ---
  processed_path <- file.path(
    paths$data_processed, 
    paste0(config$general$dataset_name, "_promoter_annotation.rds")
  )
  saveRDS(promoters, processed_path)
  message("Wrote promoter annotation to: ", processed_path)
  
  return(promoters)
}
# ===| C. Generate Background Genes |===========================================

# --- Generates background genes for pathway enrichment analysis ---
generate_background_genes <- function(cpg_anno, promoters, paths, config) {
  
  message("Generating background genes...")
  
  # --- Generate universe as all genes with atleast one CpG in its promoter ---
  overlaps <- findOverlaps(cpg_anno, promoters, ignore.strand = TRUE)
  background_genes <- unique(promoters$symbol[subjectHits(overlaps)])
  
  # --- Save processed results ---
  write_xlsx(as.data.frame(background_genes), file.path(paths$results_integration, paste0(
    config$general$dataset_name, "_background_genes.xlsx"
  )))
  
  message("Successfully generated background genes")
}