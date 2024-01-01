#' annotationDatabases
#' @title Load annotation databases
#' @description Assigns Bioconductor annotation databases (BSgenome, TxDb, org.db).
#' @param genome Character string of genome symbol (i.e. "hg38").
#' @param EnsDb A logical indicating whether Ensembl annotations should be used. 
#' @return BSgenome, TxDb, org.db for genome of interest are loaded and
#'  assigned to the global environment.
#' @import BiocManager
#' @importFrom dplyr case_when
#' @importFrom glue glue
#' @import AnnotationHub
#' @importFrom utils installed.packages
#' @rawNamespace import(ensembldb, except = c(select, filter))
#' @export annotationDatabases
#' 
annotationDatabases <- function(genome = genome,
                                EnsDb = FALSE){
  packages <- dplyr::case_when(genome == "hg38" ~ c("BSgenome.Hsapiens.UCSC.hg38.masked",
                                                    "TxDb.Hsapiens.UCSC.hg38.knownGene",
                                                    "org.Hs.eg.db"),
                               genome == "hg19" ~ c("BSgenome.Hsapiens.UCSC.hg19.masked",
                                                    "TxDb.Hsapiens.UCSC.hg19.knownGene",
                                                    "org.Hs.eg.db"),
                               genome == "mm10" ~ c("BSgenome.Mmusculus.UCSC.mm10.masked",
                                                    "TxDb.Mmusculus.UCSC.mm10.knownGene",
                                                    "org.Mm.eg.db"),
                               genome == "mm9" ~ c("BSgenome.Mmusculus.UCSC.mm9.masked",
                                                   "TxDb.Mmusculus.UCSC.mm9.knownGene",
                                                   "org.Mm.eg.db"),
                               genome == "rheMac10" ~ c("BSgenome.Mmulatta.UCSC.rheMac10",
                                                        "TxDb.Mmulatta.UCSC.rheMac10.refGene",
                                                        "org.Mmu.eg.db"),
                               genome == "rheMac8" ~ c("BSgenome.Mmulatta.UCSC.rheMac8",
                                                       "TxDb.Mmulatta.UCSC.rheMac8.refGene",
                                                       "org.Mmu.eg.db"),
                               genome == "rn6" ~ c("BSgenome.Rnorvegicus.UCSC.rn6",
                                                   "TxDb.Rnorvegicus.UCSC.rn6.refGene",
                                                   "org.Rn.eg.db"),
                               genome == "danRer11" ~ c("BSgenome.Drerio.UCSC.danRer11",
                                                        "TxDb.Drerio.UCSC.danRer11.refGene",
                                                        "org.Dr.eg.db"),
                               genome == "galGal6" ~ c("BSgenome.Ggallus.UCSC.galGal6",
                                                       "TxDb.Ggallus.UCSC.galGal6.refGene",
                                                       "org.Gg.eg.db"),
                               genome == "bosTau9" ~ c("BSgenome.Btaurus.UCSC.bosTau9",
                                                       "TxDb.Btaurus.UCSC.bosTau9.refGene",
                                                       "org.Bt.eg.db"),
                               genome == "panTro6" ~ c("BSgenome.Ptroglodytes.UCSC.panTro6",
                                                       "TxDb.Ptroglodytes.UCSC.panTro6.refGene",
                                                       "org.Pt.eg.db"),
                               genome == "dm6" ~ c("BSgenome.Dmelanogaster.UCSC.dm6",
                                                   "TxDb.Dmelanogaster.UCSC.dm6.ensGene",
                                                   "org.Dm.eg.db"),
                               genome == "susScr11" ~ c("BSgenome.Sscrofa.UCSC.susScr11",
                                                        "TxDb.Sscrofa.UCSC.susScr11.refGene",
                                                        "org.Ss.eg.db"),
                               genome == "canFam3" ~ c("BSgenome.Cfamiliaris.UCSC.canFam3.masked",
                                                       "TxDb.Cfamiliaris.UCSC.canFam3.refGene",
                                                       "org.Cf.eg.db"),
                               # TAIR10 is an "annotation release" based on TAIR9.
                               genome == "TAIR10" ~ c("BSgenome.Athaliana.TAIR.TAIR9",
                                                      "TxDb.Athaliana.BioMart.plantsmart28",
                                                      "org.At.tair.db"),
                               genome == "TAIR9" ~ c("BSgenome.Athaliana.TAIR.TAIR9",
                                                     "TxDb.Athaliana.BioMart.plantsmart28",
                                                     "org.At.tair.db"),
                               genome == "Dpulex" ~ c("BSgenome.Dpulex.NCBI.ASM2113471v1",
                                                      "TxDb.Dpulex.NCBI.ASM2113471v1.knownGene",
                                                      "org.Dpulex.eg.db")
  )
  
  new.packages <- packages[!(packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)>0){
    glue::glue("Installing {new.packages}")
    
    if(genome == "Dpulex"){
      options(timeout = 200)
        if("BSgenome.Dpulex.NCBI.ASM2113471v1" %in% new.packages){
          # BSgenome.Dpulex.NCBI.ASM2113471v1-seed must be in working directory
          # Download fna.gz sequence and export as 2bit file
          if(!file.exists("daphnia_pulex.fasta.gz")){
            download.file(url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/134/715/GCF_021134715.1_ASM2113471v1/GCF_021134715.1_ASM2113471v1_genomic.fna.gz", destfile = "daphnia_pulex.fasta.gz")
          }
          pulex_fasta = Biostrings::readDNAStringSet("daphnia_pulex.fasta.gz")
          pulex_2bit = file.path(getwd(), "daphnia_pulex.2bit")
          rtracklayer::export.2bit(pulex_fasta, pulex_2bit)
          
          BSgenome::forgeBSgenomeDataPkg("BSgenome.Dpulex.NCBI.ASM2113471v1-seed", verbose = TRUE)
          system('R CMD build BSgenome.Dpulex.NCBI.ASM2113471v1/')
          system('R CMD check BSgenome.Dpulex.NCBI.ASM2113471v1_1.0.0.tar.gz')
          system('R CMD INSTALL BSgenome.Dpulex.NCBI.ASM2113471v1_1.0.0.tar.gz')
        } 
        
        if("TxDb.Dpulex.NCBI.ASM2113471v1.knownGene" %in% new.packages) {
          
          # Gene Transfer Format (GTF) file
          if(!file.exists("daphnia_pulex.gtf.gz")){
            download.file(url = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/134/715/GCF_021134715.1_ASM2113471v1/GCF_021134715.1_ASM2113471v1_genomic.gtf.gz",
                          destfile = "daphnia_pulex.gtf.gz")
          }
          # Chromosome data
          if(!file.exists("daphnia_pulex.chrom.sizes.txt")){
            download.file(url = "https://hgdownload.soe.ucsc.edu/hubs/GCF/021/134/715/GCF_021134715.1/GCF_021134715.1.chrom.sizes.txt",
                          destfile = "daphnia_pulex.chrom.sizes.txt")
          }
          
          chrom_info = read.csv(file = "daphnia_pulex.chrom.sizes.txt",
                                header = FALSE, sep = "\t",
                                col.names = c("chr","size"))
          chrom_info = rbind(chrom_info[order(chrom_info$chr[1:12]),],chrom_info[13,])
          
          seqinfo_Dpulex = GenomeInfoDb::Seqinfo(seqnames = chrom_info$chr,
                                                 seqlengths = chrom_info$size,
                                                 isCircular = logical(13),
                                                 genome = "Dpulex")
          
          # Build metadata dataframe
          name = c("Resource URL", "Type of Gene ID", "exon_nrow", "cds_nrow")
          value = c("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/021/134/715/GCF_021134715.1_ASM2113471v1/", "Entrez Gene ID", "159649", "113453")
          
          pulex_metadata = data.frame(name, value)
          
          TxDb_pulex = GenomicFeatures::makeTxDbFromGFF(file = "daphnia_pulex.gtf.gz",
                                                        dataSource = "NCBI",
                                                        organism = "Daphnia pulex",
                                                        taxonomyId = 6669,
                                                        chrominfo = seqinfo_Dpulex,
                                                        metadata = pulex_metadata)
          GenomicFeatures::makeTxDbPackage(TxDb_pulex,
                                           version = "1.0",
                                           maintainer = "Wassim Salam <wassimsalam49@gmail.com>",
                                           author = "Wassim Salam",
                                           destDir = ".",
                                           pkgname = "TxDb.Dpulex.NCBI.ASM2113471v1.knownGene")
          install.packages("TxDb.Dpulex.NCBI.ASM2113471v1.knownGene/", repos = NULL, type = "source", quiet = TRUE)
          
        }
        
        if("org.Dpulex.eg.db" %in% new.packages){
          org.Dp.eg.db = AnnotationHub::AnnotationHub()[["AH115573"]]
          GIDkeys = keys(org.Dp.eg.db,"GID")
          DpSym = AnnotationDbi::select(org.Dp.eg.db,
                                        keys = GIDkeys,
                                        columns = c("GID", "ENTREZID", "SYMBOL", "GENENAME"))
          
          DpChr = na.omit(AnnotationDbi::select(org.Dp.eg.db,
                                                keys = GIDkeys,
                                                columns = c("GID", "CHR")))
          
          DpGO = na.omit(AnnotationDbi::select(org.Dp.eg.db,
                                               keys = GIDkeys,
                                               columns = c("GID", "GO", "EVIDENCE")))
          
          AnnotationForge::makeOrgPackage(gene_info = DpSym, chromosome = DpChr, go = DpGO,
                                          version = "1.0",
                                          maintainer = "Wassim Salam <wassimsalam49@gmail.com>",
                                          author = "Wassim Salam <wassimsalam49@gmail.com>",
                                          outputDir = ".",
                                          tax_id = "6669",
                                          genus = "Daphnia",
                                          species = "pulex",
                                          goTable = "go")
          
          install.packages("org.Dpulex.eg.db/", repos = NULL, type = "source", quiet = TRUE)
        } 
        
      devtools::install_github("wassimsalam01/annotatr", force = TRUE)
      devtools::install_github("wassimsalam01/dmrseq", force = TRUE)
      devtools::install_github("wassimsalam01/ChIPseeker", force = TRUE)
    }else{
      suppressMessages(BiocManager::install(new.packages, update = FALSE, ask = FALSE, quiet = TRUE))
    }
    cat("Done", "\n")
  }
  print(glue::glue("Loading {packages}"))
  stopifnot(suppressMessages(sapply(packages, require, character.only = TRUE)))
  
  if(genome == "hg38"){
    assign("goi", BSgenome.Hsapiens.UCSC.hg38.masked, envir = parent.frame())
    assign("TxDb", TxDb.Hsapiens.UCSC.hg38.knownGene, envir = parent.frame())
    assign("annoDb", "org.Hs.eg.db", envir = parent.frame())
  }else if(genome == "hg19"){
    assign("goi", BSgenome.Hsapiens.UCSC.hg19.masked, envir = parent.frame()) 
    assign("TxDb", TxDb.Hsapiens.UCSC.hg19.knownGene, envir = parent.frame())
    assign("annoDb", "org.Hs.eg.db", envir = parent.frame())
  }else if(genome == "mm10"){
    assign("goi", BSgenome.Mmusculus.UCSC.mm10.masked, envir = parent.frame())
    assign("TxDb", TxDb.Mmusculus.UCSC.mm10.knownGene, envir = parent.frame())
    assign("annoDb", "org.Mm.eg.db", envir = parent.frame())
  }else if(genome == "mm9"){
    assign("goi", BSgenome.Mmusculus.UCSC.mm9.masked, envir = parent.frame())
    assign("TxDb", TxDb.Mmusculus.UCSC.mm9.knownGene, envir = parent.frame())
    assign("annoDb", "org.Mm.eg.db", envir = parent.frame())
  }else if(genome == "rheMac10"){
    assign("goi", BSgenome.Mmulatta.UCSC.rheMac10, envir = parent.frame())
    assign("TxDb", TxDb.Mmulatta.UCSC.rheMac10.refGene, envir = parent.frame())
    assign("annoDb", "org.Mmu.eg.db", envir = parent.frame())
  }else if(genome == "rheMac8"){
    assign("goi", BSgenome.Mmulatta.UCSC.rheMac8, envir = parent.frame())
    assign("TxDb", TxDb.Mmulatta.UCSC.rheMac8.refGene, envir = parent.frame())
    assign("annoDb", "org.Mmu.eg.db", envir = parent.frame())
  }else if(genome == "rn6"){
    assign("goi", BSgenome.Rnorvegicus.UCSC.rn6, envir = parent.frame())
    assign("TxDb", TxDb.Rnorvegicus.UCSC.rn6.refGene, envir = parent.frame())
    assign("annoDb", "org.Rn.eg.db", envir = parent.frame())
  }else if(genome == "danRer11"){
    assign("goi", BSgenome.Drerio.UCSC.danRer11, envir = parent.frame())
    assign("TxDb", TxDb.Drerio.UCSC.danRer11.refGene, envir = parent.frame())
    assign("annoDb", "org.Dr.eg.db", envir = parent.frame())
  }else if(genome == "galGal6"){
    assign("goi", BSgenome.Ggallus.UCSC.galGal6, envir = parent.frame())
    assign("TxDb", TxDb.Ggallus.UCSC.galGal6.refGene, envir = parent.frame())
    assign("annoDb", "org.Gg.eg.db", envir = parent.frame())
  }else if(genome == "bosTau9"){
    assign("goi", BSgenome.Btaurus.UCSC.bosTau9, envir = parent.frame())
    assign("TxDb", TxDb.Btaurus.UCSC.bosTau9.refGene, envir = parent.frame())
    assign("annoDb", "org.Bt.eg.db", envir = parent.frame())
  }else if(genome == "panTro6"){
    assign("goi", BSgenome.Ptroglodytes.UCSC.panTro6, envir = parent.frame())
    assign("TxDb", TxDb.Ptroglodytes.UCSC.panTro6.refGene, envir = parent.frame())
    assign("annoDb", "org.Pt.eg.db", envir = parent.frame())
  }else if(genome == "dm6"){
    assign("goi", BSgenome.Dmelanogaster.UCSC.dm6, envir = parent.frame())
    assign("TxDb", TxDb.Dmelanogaster.UCSC.dm6.ensGene, envir = parent.frame())
    assign("annoDb", "org.Dm.eg.db", envir = parent.frame())
  }else if(genome == "susScr11"){
    assign("goi", BSgenome.Sscrofa.UCSC.susScr11, envir = parent.frame())
    assign("TxDb", TxDb.Sscrofa.UCSC.susScr11.refGene, envir = parent.frame())
    assign("annoDb", "org.Ss.eg.db", envir = parent.frame())
  }else if(genome == "canFam3"){
    assign("goi", BSgenome.Cfamiliaris.UCSC.canFam3.masked, envir = parent.frame())
    assign("TxDb", TxDb.Cfamiliaris.UCSC.canFam3.refGene, envir = parent.frame())
    assign("annoDb", "org.Cf.eg.db", envir = parent.frame())
  }else if(genome %in% c("TAIR9", "TAIR10")){
    assign("goi", BSgenome.Athaliana.TAIR.TAIR9, envir = parent.frame())
    assign("TxDb", TxDb.Athaliana.BioMart.plantsmart28, envir = parent.frame())
    assign("annoDb", "org.At.tair.db", envir = parent.frame())
  }else if(genome == "Dpulex"){
    assign("goi", BSgenome.Dpulex.NCBI.ASM2113471v1, envir = parent.frame())
    assign("TxDb", TxDb.Dpulex.NCBI.ASM2113471v1.knownGene, envir = parent.frame())
    assign("annoDb", "org.Dpulex.eg.db", envir = parent.frame()) 
  }else{
    stop(glue::glue("{genome} is not supported, please choose either hg38, hg19, mm10, mm9, \\
    rheMac10, rheMac8, rn6, danRer11, galGal6, bosTau9, panTro6, dm6, susScr11, canFam3, TAIR10, \\
    TAIR9 or Dpulex [Case Sensitive]"))
  }
  
  if(EnsDb == TRUE){
    if(genome %in% c("hg38", "mm10", "rheMac10", "rheMac8", "rn6", "danRer11", "galGal6",
                     "bosTau9", "dm6", "susScr11", "canFam3")){
      
      print(glue::glue("EnsemblDb annotations for {genome} will be used."))
      # AnnotationHub::display(AnnotationHub::query(AnnotationHub::AnnotationHub(), c("EnsDb", "Pan troglodytes")))
      
      ahCode <- dplyr::case_when(genome == "hg38" ~ "AH83216",
                                 genome == "mm10" ~ "AH83247",
                                 genome == "rheMac10" ~ "AH83244",
                                 genome == "rheMac8" ~ "AH73903",
                                 genome == "rn6" ~ "AH83318",
                                 genome == "danRer11" ~ "AH83189",
                                 genome == "galGal6" ~ "AH83209",
                                 genome == "bosTau9" ~ "AH83145",
                                 genome == "dm6" ~ "AH83185",
                                 genome == "susScr11" ~ "AH83340",
                                 genome == "canFam3" ~ "AH78741")
      
      print(glue::glue("Your AnnotationHub code is {ahCode}."))
      assign("TxDb", AnnotationHub::AnnotationHub()[[ahCode]], envir = parent.frame())
      
    }else if(genome %in% c("hg19", "mm9", "panTro6", "TRAIR9", "TAIR10")){
      
      print(glue::glue("EnsemblDB annotations for {genome} are not available, \\
                   the default Bioconductor TxDb will be used."))
    }
  }
  
}
