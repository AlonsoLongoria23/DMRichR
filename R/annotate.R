#' annotateRegions
#' @title Annotate DMRs and blocks
#' @description Annotate and tidy regions from \code{dmrseq::dmrseq()}
#' @param regions A \code{GRanges} object of DMRs, blocks, or background regions from \code{dmrseq::dmrseq()}
#' @param TxDb \code{TxDb} or \code{EnsDb} annotation package for genome of interest
#' @param annoDb Character specifying \code{OrgDb} annotation package for species of interest
#' @return A \code{tibble} of annotated regions
#' @rawNamespace import(ensembldb, except = c(select, filter))
#' @importFrom dplyr rename_with as_tibble case_when mutate select recode_factor distinct
#' @importFrom tidyselect any_of
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom ChIPseeker annotatePeak
#' @importFrom magrittr %>%
#' @importFrom glue glue glue_collapse
#' @importFrom GenomeInfoDb genome seqlevelsStyle
#' @importFrom plyranges join_overlap_left select mutate
#' @export annotateRegions
#' 
annotateRegions <- function(regions = sigRegions,
                            TxDb = TxDb,
                            annoDb = annoDb){
  
  genome <- TxDb %>%
    GenomeInfoDb::genome() %>%
    unique()
  
  print(glue::glue("Annotating {tidyRegions} regions from {genome} with gene symbols",
                   tidyRegions = length(regions)))
  
  if(is(TxDb, "EnsDb")){
    
    genome <- dplyr::case_when(GenomeInfoDb::genome(TxDb) == "GRCh38" ~ "hg38",
                               GenomeInfoDb::genome(TxDb) == "GRCm38" ~ "mm10",
                               GenomeInfoDb::genome(TxDb) == "Mmul_10" ~ "rheMac10",
                               GenomeInfoDb::genome(TxDb) == "Mmul_8.0.1" ~ "rheMac8",
                               GenomeInfoDb::genome(TxDb) == "Rnor_6.0" ~ "rn6",
                               GenomeInfoDb::genome(TxDb) == "GRCz11" ~ "danRer11",
                               GenomeInfoDb::genome(TxDb) == "GRCg6a" ~ "galGal6",
                               GenomeInfoDb::genome(TxDb) == "ARS-UCD1.2" ~ "bosTau9",
                               GenomeInfoDb::genome(TxDb) == "BDGP6.28" ~ "dm6",
                               GenomeInfoDb::genome(TxDb) == "Sscrofa11.1" ~ "susScr11",
                               GenomeInfoDb::genome(TxDb) == "CanFam3.1" ~ "canFam3") %>%
      unique()
                               }
  
  CpGs <- DMRichR::getCpGs(genome)
  
  regionsCpG <- regions %>% 
    plyranges::join_overlap_left(CpGs %>%
                                   dplyr::filter(type == "islands") %>% 
                                   dplyr::select(CpG.Island = type)) %>%
    unique() %>% 
    plyranges::join_overlap_left(CpGs %>%
                                   dplyr::filter(type == "shores") %>% 
                                   dplyr::select(CpG.Shore = type)) %>%
    unique() %>% 
    plyranges::join_overlap_left(CpGs %>%
                                   dplyr::filter(type == "shelves") %>% 
                                   dplyr::select(CpG.Shelf = type)) %>%
    unique() %>% 
    plyranges::join_overlap_left(CpGs %>%
                                   dplyr::filter(type == "inter") %>% 
                                   dplyr::select(Open.Sea = type)) %>%
    unique() %>% 
    dplyr::mutate(CpG.Island = dplyr::case_when(CpG.Island == "islands" ~ "Yes",
                                                    TRUE ~ "No"),
                      CpG.Shore = dplyr::case_when(CpG.Shore == "shores" ~ "Yes",
                                                   TRUE ~ "No"),
                      CpG.Shelf = dplyr::case_when(CpG.Shelf == "shelves" ~ "Yes",
                                                   TRUE ~ "No"),
                      Open.Sea = dplyr::case_when(Open.Sea == "inter" ~ "Yes",
                                                  TRUE ~ "No"))
  
  if(is(TxDb, "EnsDb")){
    GenomeInfoDb::seqlevelsStyle(regionsCpG) <- "Ensembl" # Work around for organism not supported
  }
  
  regionsCpG %>% 
    ChIPseeker::annotatePeak(TxDb = TxDb,
                             annoDb = annoDb,
                             overlap = "all",
                             verbose = FALSE) %>%
    dplyr::as_tibble() %>%
    dplyr::select(-tidyselect::any_of(c("strand",
                                        "index.start",
                                        "index.end",
                                        "index.width",
                                        "area",
                                        "geneId",
                                        "geneChr",
                                        "geneStart",
                                        "geneEnd",
                                        "geneLength",
                                        "geneStrand",
                                        "transcriptId",
                                        "transcriptBiotype",
                                        "ENTREZID"))) %>%
    dplyr::mutate(annotation = gsub(" \\(.*","", annotation)) %>%
    dplyr::rename_with(
      ~ dplyr::case_when(
        . == "seqnames" ~ "chr",
        . == "L" ~ "CpGs",
        . == "beta" ~ "betaCoefficient",
        . == "stat" ~ "statistic",
        . == "pval" ~ "p.value",
        . == "qval" ~ "q.value",
        . == "SYMBOL" ~ "geneSymbol",
        . == "GENENAME" ~ "gene",
        TRUE ~ .)) %>% 
    return()
}

#' DMReport
#' @title Create an html report of DMRs or blocks
#' @description Create an html report of significant regions from \code{dmrseq}
#' @param sigRegions \code{GRanges} object of significant regions (DMRs or blocks) from \code{dmrseq} that 
#' were annotated by \code{DMRichR::annotateRegions}
#' @param regions \code{GRanges} object of background regions from \code{dmrseq}
#' @param bs.filtered Filtered \code{bsseq} object from \code{processBismark()}
#' @param coverage Numeric of coverage samples were filtered for
#' @param name Character for html report name
#' @return Saves an html report of DMRs with genic annotations
#' @importFrom gt gt tab_header fmt_number fmt_scientific fmt_percent as_raw_html
#' @importFrom dplyr select mutate 
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importClassesFrom bsseq BSseq 
#' @export DMReport
#' 
DMReport <- function(sigRegions = sigRegions,
                     regions = regions,
                     bs.filtered = bs.filtered,
                     coverage = coverage,
                     name = "DMReport"){
  cat("\n","Preparing HTML report...")
  
  stopifnot(class(sigRegions) == c("tbl_df", "tbl", "data.frame"))
  
  required_columns <- c("chr",
                        "start",
                        "end",
                        "width",
                        "CpGs",
                        "betaCoefficient",
                        "statistic",
                        "p.value",
                        "q.value",
                        "difference",
                        "CpG.Island",
                        "CpG.Shore",
                        "CpG.Shelf",
                        "Open.Sea",
                        "annotation",
                        "distanceToTSS")
  optional_columns <- c("geneSymbol", "gene")

  missing_optional <- setdiff(optional_columns, colnames(sigRegions))
  if (length(missing_optional) > 0) {
    for (missing_col in missing_optional) {
      sigRegions[[missing_col]] <- rep(NA_character_, nrow(sigRegions))
    }
  }

  sigRegions %>%
    dplyr::select(tidyselect::all_of(required_columns),
                  tidyselect::any_of(optional_columns)) %>%
    dplyr::mutate(difference = difference/100) %>% 
    gt::gt() %>%
    gt::tab_header(
      title = name,
      subtitle = glue::glue("There are {tidySigRegions} regions \\
             ({tidyHyper}% hypermethylated, {tidyHypo}% hypomethylated) \\
             from {tidyRegions} background regions consisting of {tidyCpGs} CpGs \\
             assayed at {coverage}x coverage.
             On average, the DMRs are {avgLength} bp long and contain {avgCpGs} CpGs.", 
                            tidySigRegions = nrow(sigRegions),
                            tidyHyper = round(sum(sigRegions$statistic > 0) / nrow(sigRegions),
                                              digits = 2)*100,
                            tidyHypo = round(sum(sigRegions$statistic < 0) / nrow(sigRegions),
                                             digits = 2)*100,
                            tidyRegions = length(regions),
                            tidyCpGs = nrow(bs.filtered),
                            avgLength = mean(sigRegions$width) %>% round(),
                            avgCpGs = mean(sigRegions$CpGs) %>% round()
                            )) %>% 
    gt::fmt_number(
      columns = gt::vars("width", "CpGs"),
      decimals = 0
      ) %>% 
    gt::fmt_scientific(
      columns = gt::vars("p.value", "q.value"),
      decimals = 2
      ) %>%
    gt::fmt_percent(
      columns = gt::vars("difference"),
      drop_trailing_zeros = TRUE
      ) %>% 
    gt::as_raw_html(inline_css = FALSE) %>%
    write(glue::glue("{name}.html"))
  cat("Done", "\n")
}

#' getExons
#' @title Obtain exons for plotting
#' @description Obtain exon annotations from a \code{ensDb} object and format for \code{plotDMRs()} 
#' @param TxDb A \code{ensDb} object
#' @return A \code{GRanges} object of annotated exons for every gene with a symbol in the genome.
#' @rawNamespace import(ensembldb, except = c(select, filter))
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom BiocGenerics unlist
#' @importFrom plyranges mutate select filter
#' @importFrom GenomeInfoDb genome
#' @references Based on \code{annotatr::build_gene_annots()},
#'  see: \url{https://github.com/rcavalcante/annotatr/blob/master/R/build_annotations.R}
#' @export getExons
#' 
getExons <- function(TxDb = TxDb){
  
  stopifnot(is(TxDb, "EnsDb"))
  
  message('Building exons...')
  
  exons <- TxDb %>%
    ensembldb::cdsBy(by = "tx",
                     columns = c("tx_id", "gene_id", "symbol") #, # listColumns(TxDb)
                     # filter = GeneBiotypeFilter("protein_coding")
                     ) %>%
    BiocGenerics::unlist(use.names = FALSE) %>%
    dplyr::mutate(id = glue::glue("CDS:{seq_along(.)}"),
                      type = glue::glue("{unique(genome(TxDb))}_genes_cds")
                      ) %>%
    dplyr::select(id, tx_id, gene_id, symbol, type) # %>%
    # dplyr::filter(symbol != "")
  
  GenomeInfoDb::genome(exons) <- NA  
  ensembldb::seqlevelsStyle(exons) <- "UCSC"
  
  return(exons)
}

#' getCpGs
#' @title Obtain CpG island, CpG shore, CpG shelf, and open sea annotations
#' @description Obtain UCSC CpG islands and build CpG shore, CpG shelf, and open sea annotations.
#'  This function is based on \code{annotatr:::build_cpg_annots()}; however, 
#'  it obtains annotations for all genomes in the UCSC genome browser.
#' @param genome Character specifying the genome
#' @return A \code{GRanges} object of CpG island, CpG shore, CpG shelf, and open sea annotations
#' @importFrom GenomicRanges makeGRangesFromDataFrame mcols trim setdiff sort gaps
#' @importFrom GenomeInfoDb keepStandardChromosomes
#' @importFrom glue glue
#' @importFrom magrittr %>%
#' @importFrom plyranges stretch mutate select
#' @references Based on \code{annotatr:::build_cpg_annots()},
#'  see: \url{https://github.com/rcavalcante/annotatr/blob/master/R/build_annotations.R}
#' @export getCpGs
#' 

#' getGenes
#' @title Obtain exon annotations for plotting from a TxDb or EnsDb
#' @description Obtain exon/CDS annotations and format for plotDMRs()
#' @param TxDb A TxDb or EnsDb annotation object
#' @param annoDb Character string naming the OrgDb package, or an OrgDb object
#' @return A GRanges object of annotated exons for plotting
#' @rawNamespace import(ensembldb, except = c(select, filter))
#' @importFrom GenomicFeatures exonsBy
#' @importFrom AnnotationDbi select keytypes columns
#' @importFrom BiocGenerics unlist
#' @importFrom GenomeInfoDb genome
#' @importFrom glue glue
#' @importFrom dplyr mutate select
#' @export getGenes
getGenes <- function(TxDb = TxDb, annoDb = annoDb){

  stopifnot(is(TxDb, "TxDb") || is(TxDb, "EnsDb"))

  message("Building genes/exons for plotting...")

  if (is.character(annoDb)) {
    suppressPackageStartupMessages(require(annoDb, character.only = TRUE))
    annoDb_obj <- get(annoDb)
  } else {
    annoDb_obj <- annoDb
  }

  genome_name <- unique(GenomeInfoDb::genome(TxDb))
  genome_name <- genome_name[!is.na(genome_name)]
  if (length(genome_name) == 0) {
    genome_name <- "unknownGenome"
  } else {
    genome_name <- genome_name[1]
  }

  if (is(TxDb, "EnsDb")) {

    exons <- TxDb %>%
      ensembldb::cdsBy(
        by = "tx",
        columns = c("tx_id", "gene_id", "symbol")
      ) %>%
      BiocGenerics::unlist(use.names = FALSE) %>%
      dplyr::mutate(
        id = glue::glue("CDS:{seq_along(.)}"),
        type = glue::glue("{genome_name}_genes_cds")
      ) %>%
      dplyr::select(id, tx_id, gene_id, symbol, type)

    GenomeInfoDb::genome(exons) <- NA
    ensembldb::seqlevelsStyle(exons) <- "UCSC"

    return(exons)
  }

  exon_list <- GenomicFeatures::exonsBy(TxDb, by = "gene")
  exons <- BiocGenerics::unlist(exon_list, use.names = FALSE)

  gene_ids <- names(exon_list)
  exon_gene_id <- rep(gene_ids, lengths(exon_list))

  available_keytypes <- AnnotationDbi::keytypes(annoDb_obj)
  keytype_to_use <- if ("GID" %in% available_keytypes) "GID" else "ENTREZID"

  available_cols <- AnnotationDbi::columns(annoDb_obj)
  wanted_cols <- intersect(c("SYMBOL", "GENENAME"), available_cols)

  anno_tbl <- AnnotationDbi::select(
    annoDb_obj,
    keys = unique(exon_gene_id),
    keytype = keytype_to_use,
    columns = wanted_cols
  )

  anno_tbl <- anno_tbl[!duplicated(anno_tbl[[keytype_to_use]]), , drop = FALSE]
  rownames(anno_tbl) <- anno_tbl[[keytype_to_use]]

  symbol_vec <- if ("SYMBOL" %in% colnames(anno_tbl)) {
    anno_tbl[exon_gene_id, "SYMBOL", drop = TRUE]
  } else {
    rep(NA_character_, length(exon_gene_id))
  }

  tx_ids <- rep(NA_character_, length(exons))
  if ("tx_name" %in% names(S4Vectors::mcols(exons))) {
    tx_ids <- as.character(S4Vectors::mcols(exons)$tx_name)
  } else if ("tx_id" %in% names(S4Vectors::mcols(exons))) {
    tx_ids <- as.character(S4Vectors::mcols(exons)$tx_id)
  }

  S4Vectors::mcols(exons)$id <- glue::glue("exon:{seq_along(exons)}")
  S4Vectors::mcols(exons)$tx_id <- tx_ids
  S4Vectors::mcols(exons)$gene_id <- exon_gene_id
  S4Vectors::mcols(exons)$symbol <- symbol_vec
  S4Vectors::mcols(exons)$type <- glue::glue("{genome_name}_genes_exons")

  exons <- dplyr::select(exons, id, tx_id, gene_id, symbol, type)

  GenomeInfoDb::genome(exons) <- NA

  return(exons)
}

getCpGs <- function(genome = genome){

  if(genome == "Dpulex"){

    if(!file.exists("daphnia_pulex.chrom.sizes.txt")){
            download.file(url = "https://hgdownload.soe.ucsc.edu/hubs/GCF/021/134/715/GCF_021134715.1/GCF_021134715.1.chrom.sizes.txt",
                          destfile = "daphnia_pulex.chrom.sizes.txt")
    }
    
    chrom_info = read.csv(file = "daphnia_pulex.chrom.sizes.txt",
                          header = FALSE, sep = "\t",
                          col.names = c("chr","size"))
    chrom_info = rbind(chrom_info[order(chrom_info$chr[1:12]),],chrom_info[13,])
    
    message('Building CpG islands...')

    # CGI-Dpulex.txt must be in work directory
    islands = read.csv("CGI-Dpulex.txt",
                   header = TRUE,
                   sep = "\t",
                   colClasses = c(rep(NA,3), rep("NULL",5))) %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE, seqinfo = GenomeInfoDb::Seqinfo(seqnames = chrom_info$chr[1:12],
                                                                                                       seqlengths = chrom_info$size[1:12],
                                                                                                       isCircular = logical(12),
                                                                                                       genome = "Dpulex")) %>%
    # GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>% Dpulex not supported by this function
    dplyr::mutate(id = glue::glue("island:{seq_along(.)}"),
                      type = "islands")
  
  message('Building CpG shores...')
  
  shores <- islands %>% 
    plyranges::stretch(1000) %>% 
    GenomicRanges::trim() %>%
    GenomicRanges::setdiff(islands) %>%
    dplyr::mutate(id = glue::glue("shore:{seq_along(.)}"),
                      type = "shores")
  
  message('Building CpG shelves...')
  
  shelves <- shores %>% 
    plyranges::stretch(1000) %>% 
    GenomicRanges::trim() %>%
    GenomicRanges::setdiff(islands) %>%
    GenomicRanges::setdiff(shores) %>%
    dplyr::mutate(id = glue::glue("shelf:{seq_along(.)}"),
                      type = "shelves")
  
  message('Building inter-CpG-islands...')
  
  inter_cgi <- c(islands, shores, shelves) %>%
    GenomicRanges::sort() %>%
    GenomicRanges::gaps() %>%
    dplyr::mutate(id = glue::glue("inter:{seq_along(.)}"),
                      type = "inter")
  
  c(islands, shores, shelves, inter_cgi) %>%
    GenomicRanges::sort() %>%
    dplyr::mutate(tx_id = NA,
                      gene_id = NA,
                      symbol = NA) %>%
    dplyr::select(id, tx_id, gene_id, symbol, type) %>% 
    return()
    
} else if (genome %in% c("Tthymallus", "ThyArc1.0")) {
    
    # Check if required local files exist in the working directory
    if (genome == "ThyArc1.0") {
      cgi_file <- "CGI-Tarcticus.txt"
      size_file <- "Thymallus_arc_chr_sizes.txt"
        } else {
      cgi_file <- "CGI-Thymallus.txt"
      size_file <- "Thymallus_chr_sizes.txt"
        }

    # 2. Check if required local files exist
    if (!file.exists(size_file)) {
      stop(paste("Chromosome sizes file not found:", size_file))
    }
    if (!file.exists(cgi_file)) {
      stop(paste("CpG islands file not found:", cgi_file))
    }
    
# 3. Read chromosome info
    chrom_info <- read.delim(file = size_file, header = TRUE, sep = "\t")
    # Ensure consistent ordering
    chrom_info <- chrom_info[order(chrom_info$chr), ]

    message(glue::glue('Building CpG islands for {genome}...'))

    # 4. Read CpG islands from the local file
    islands <- read.delim(cgi_file, header = TRUE, sep = "\t") %>%
      GenomicRanges::makeGRangesFromDataFrame(
        keep.extra.columns = TRUE,
        seqinfo = GenomeInfoDb::Seqinfo(
          seqnames   = chrom_info$chr,
          seqlengths = chrom_info$size,
          isCircular = logical(nrow(chrom_info)),
          genome     = genome
        )
      ) %>%
      dplyr::mutate(id = glue::glue("island:{seq_along(.)}"),
                        type = "islands")

    message('Building CpG shores (4kb stretch)...')
    shores <- islands %>% 
      plyranges::stretch(4000) %>% 
      GenomicRanges::trim() %>%
      GenomicRanges::setdiff(islands) %>%
      dplyr::mutate(id = glue::glue("shore:{seq_along(.)}"),
                        type = "shores")

    message('Building CpG shelves (4kb stretch)...')
    shelves <- shores %>% 
      plyranges::stretch(4000) %>% 
      GenomicRanges::trim() %>%
      GenomicRanges::setdiff(islands) %>%
      GenomicRanges::setdiff(shores) %>%
      dplyr::mutate(id = glue::glue("shelf:{seq_along(.)}"),
                        type = "shelves")

    message('Building inter-CpG-islands...')
    inter_cgi <- c(islands, shores, shelves) %>%
      GenomicRanges::sort() %>%
      GenomicRanges::gaps() %>%
      dplyr::mutate(id = glue::glue("inter:{seq_along(.)}"),
                        type = "inter")

    # 5. Final selection and return
    c(islands, shores, shelves, inter_cgi) %>%
      GenomicRanges::sort() %>%
      dplyr::mutate(tx_id   = NA,
                        gene_id = NA,
                        symbol  = NA) %>%
      dplyr::select(id, tx_id, gene_id, symbol, type) %>% 
      return()

  } else {
    message('Building CpG islands...')
  
  islands <- readr::read_tsv(glue::glue("http://hgdownload.cse.ucsc.edu/goldenpath/{genome}/database/cpgIslandExt.txt.gz"),
                             col_names = c('chr','start','end'),
                             col_types = '-cii-------') %>%
    GenomicRanges::makeGRangesFromDataFrame(keep.extra.columns = TRUE) %>%
    GenomeInfoDb::keepStandardChromosomes(pruning.mode = "coarse") %>%
    dplyr::mutate(id = glue::glue("island:{seq_along(.)}"),
                      type = "islands")
  
  message('Building CpG shores...')
  
  shores <- islands %>% 
    plyranges::stretch(4000) %>% 
    GenomicRanges::trim() %>%
    GenomicRanges::setdiff(islands) %>%
    dplyr::mutate(id = glue::glue("shore:{seq_along(.)}"),
                      type = "shores")
  
  message('Building CpG shelves...')
  
  shelves <- shores %>% 
    plyranges::stretch(4000) %>% 
    GenomicRanges::trim() %>%
    GenomicRanges::setdiff(islands) %>%
    GenomicRanges::setdiff(shores) %>%
    dplyr::mutate(id = glue::glue("shelf:{seq_along(.)}"),
                      type = "shelves")
  
  message('Building inter-CpG-islands...')
  
  inter_cgi <- c(islands, shores, shelves) %>%
    GenomicRanges::sort() %>%
    GenomicRanges::gaps() %>%
    dplyr::mutate(id = glue::glue("inter:{seq_along(.)}"),
                      type = "inter")
  
  c(islands, shores, shelves, inter_cgi) %>%
    GenomicRanges::sort() %>%
    dplyr::mutate(tx_id = NA,
                      gene_id = NA,
                      symbol = NA) %>%
    dplyr::select(id, tx_id, gene_id, symbol, type) %>% 
    return()
  }
  
}
