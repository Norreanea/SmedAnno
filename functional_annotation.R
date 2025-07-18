#!/usr/bin/env Rscript

# functional_annotation.R
# Script to perform functional annotation and identify fragmented/chimeric genes

# ---------------------------
# Function to install and load Libraries
# ---------------------------
# install_and_load <- function(packages) {
#   for (pkg in packages) {
#     if (!requireNamespace(pkg, quietly = TRUE)) {
#       if (pkg %in% rownames(installed.packages())) {
#         library(pkg, character.only = TRUE)
#       } else {
#         if (pkg %in% c("rtracklayer","GenomicFeatures", "GenomicRanges", "UniProt.ws","reactome.db")) {
#           if (!requireNamespace("BiocManager", quietly = TRUE)) {
#             install.packages("BiocManager", repos = "http://cran.us.r-project.org")
#           }
#           BiocManager::install(pkg, ask = FALSE)
#         } else {
#           install.packages(pkg, repos = "http://cran.us.r-project.org")
#         }
#         library(pkg, character.only = TRUE)
#       }
#     } else {
#       library(pkg, character.only = TRUE)
#     }
#   }
# }
# 
# ---------------------------
# Define required libraries
# ---------------------------
required_packages <- c(
  "rtracklayer",
  "dplyr",
  "data.table",
  "stringr",
  "readr",
  #"rhmmer",
  "tidyr",
  "GenomicFeatures",
  "UniProt.ws",
  "GenomicRanges",
  "reactome.db",
  "httr",
  "jsonlite",
  #"GeneAnswers",
  "purrr"

)
# ---------------------------
# Special handling for rhmmer installation --> time consuming
# ---------------------------
# # Define a personal library path inside the user's writable home directory
# personal_lib_path <- file.path(Sys.getenv("HOME"), "R_libs")
# dir.create(personal_lib_path, showWarnings = FALSE, recursive = TRUE)
# # Ensure the personal library is in the R search path
# .libPaths(c(personal_lib_path, .libPaths()))
# 
# # Check if rhmmer is installed. If not, install it from GitHub.
# if (!requireNamespace("remotes", quietly = TRUE)) {
#   install.packages(
#     "remotes",
#     repos = "https://cloud.r-project.org",
#     lib   = personal_lib_path,
#     dependencies = TRUE
#   )
# }
# library(remotes)                         
# 
# if (!requireNamespace("rhmmer", quietly = TRUE)) {
#   install_github(
#     "arendsee/rhmmer",
#     lib          = personal_lib_path,
#     dependencies = TRUE,
#     upgrade      = "never"
#   )
# }
#library(rhmmer)
# 
# ---------------------------
# Code bellow was directly taken from https://github.com/arendsee/rhmmer/blob/master/R/parse.R
# ---------------------------
read_domtblout <- function(file){
  .parse_hmmer_output(file, 'domtblout')
}


.parse_hmmer_output <- function(file, type){
  
  col_types <- if(type == 'tblout'){
    readr::cols(
      domain_name         = readr::col_character(),
      domain_accession    = readr::col_character(),
      query_name          = readr::col_character(),
      query_accession     = readr::col_character(),
      sequence_evalue     = readr::col_double(),
      sequence_score      = readr::col_double(),
      sequence_bias       = readr::col_double(),
      best_domain_evalue  = readr::col_double(),
      best_domain_score   = readr::col_double(),
      best_domain_bis     = readr::col_double(),
      domain_number_exp   = readr::col_double(),
      domain_number_reg   = readr::col_integer(),
      domain_number_clu   = readr::col_integer(),
      domain_number_ov    = readr::col_integer(),
      domain_number_env   = readr::col_integer(),
      domain_number_dom   = readr::col_integer(),
      domain_number_rep   = readr::col_integer(),
      domain_number_inc   = readr::col_character()
    )
  } else if(type == 'domtblout'){
    readr::cols(
      domain_name         = readr::col_character(),
      domain_accession    = readr::col_character(),
      domain_len          = readr::col_integer(),
      query_name          = readr::col_character(),
      query_accession     = readr::col_character(),
      qlen                = readr::col_integer(),
      sequence_evalue     = readr::col_double(),
      sequence_score      = readr::col_double(),
      sequence_bias       = readr::col_double(),
      domain_N            = readr::col_integer(),
      domain_of           = readr::col_integer(),
      domain_cevalue      = readr::col_double(),
      domain_ievalue      = readr::col_double(),
      domain_score        = readr::col_double(),
      domain_bias         = readr::col_double(),
      hmm_from            = readr::col_integer(),
      hmm_to              = readr::col_integer(),
      ali_from            = readr::col_integer(),
      ali_to              = readr::col_integer(),
      env_from            = readr::col_integer(),
      env_to              = readr::col_integer(),
      acc                 = readr::col_double()
    )
  }
  
  N <- length(col_types$cols)
  
  # the line delimiter should always be just "\n", even on Windows
  lines <- readr::read_lines(file, lazy=FALSE, progress=FALSE)
  
  table <- sub(
    pattern = sprintf("(%s).*", paste0(rep('\\S+', N), collapse=" +")),
    replacement = '\\1',
    x=lines,
    perl = TRUE
  ) %>%
    gsub(pattern="  *", replacement="\t") %>%
    paste0(collapse="\n") %>%
    readr::read_tsv(
      col_names=names(col_types$cols),
      comment='#',
      na='-',
      col_types = col_types,
      lazy=FALSE,
      progress=FALSE
    )
  
  if(type == 'domtblout'){
    table$description <- lines[!grepl("^#", lines, perl=TRUE)] %>%
      sub(
        pattern = sprintf("%s *(.*)", paste0(rep('\\S+', N), collapse=" +")),
        replacement = '\\1',
        perl = TRUE
      )
  }
  
  table
}
# # Install and load libraries
#install_and_load("rhmmer")

# Load all required libraries (from Docker), stopping if one is missing
for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    stop(paste("Package", pkg, "is not installed. Please add it to environment.yml and rebuild the Docker image."))
  }
}
# ---------------------------
# Parse command-line arguments
# ---------------------------

args <- commandArgs(trailingOnly = TRUE)

# Expected arguments:
# --annotation_dir PATH
# --functional_dir PATH
# --genome_ref PATH
# --output_dir PATH

# Initialize variables
annotation_dir <- ""
functional_dir <- ""
genome_ref <- ""
output_dir <- ""
transdecoder_gff_nucl <- ""  
transdecoder_gff_mito <- "" 

# Thresholds with default values
max_distance <- 1000
large_intron_threshold <- 100000
blast_evalue_threshold <- 1e-5
blast_identity_threshold <- 25
cevalue_threshold <- 1e-5
bit_score_threshold <- 10
coverage_threshold <- 50
interpro_evalue_threshold <- 1e-5 

# --- Parse arguments loop ---
if (length(args) > 0) {
  for (i in seq(1, length(args), by = 2)) {
    arg <- args[i]
    value <- args[i + 1]
    if (arg == "--annotation_dir") {
      annotation_dir <- value
    } else if (arg == "--functional_dir") {
      functional_dir <- value
    } else if (arg == "--genome_ref") {
      genome_ref <- value
    } else if (arg == "--output_dir") {
      output_dir <- value
    } else if (arg == "--maxDistance") {
      max_distance <- as.numeric(value)
    } else if (arg == "--largeIntronThreshold") {
      large_intron_threshold <- as.numeric(value)
    } else if (arg == "--blastEvalue") {
      blast_evalue_threshold <- as.numeric(value)
    } else if (arg == "--blastIdentity") {
      blast_identity_threshold <- as.numeric(value)
    } else if (arg == "--pfamCEvalue") {
      cevalue_threshold <- as.numeric(value)
    } else if (arg == "--pfamBitScore") {
      bit_score_threshold <- as.numeric(value)
    } else if (arg == "--pfamCoverage") {
      coverage_threshold <- as.numeric(value)
    } else if (arg == "--interproEvalue") {
      interpro_evalue_threshold <- as.numeric(value)
    } else if (arg == "--transdecoder_gff_nucl") {
      transdecoder_gff_nucl <- value
    } else if (arg == "--transdecoder_gff_mito") {
      transdecoder_gff_mito <- value
    }
  }
}
# 
# # ---------------------------
# # Validate arguments
# # ---------------------------
# # Make genome_ref optional by not requiring it
# if (annotation_dir == "" || functional_dir == "" || output_dir == "") {
#   stop("Error: Arguments --annotation_dir, --functional_dir, and --output_dir must be provided.")
# }
# 

# ---------------------------
# Import GTF annotations
# ---------------------------

genomic_gtf_path<- file.path(annotation_dir, "sanitized.final.gtf")

if (!file.exists(genomic_gtf_path)) {
  stop(paste("Error: GTF file not found at", genomic_gtf_path))
}
# Import StringTie annotations
genomic_gr  <- rtracklayer::import.gff(genomic_gtf_path)

# ---------------------------
# Map TransDecoder ORFs to genomic coordinates
# ---------------------------
map_transdecoder_to_genome <- function(td_gff_path, txdb, genomic_gr) {
  message(paste("Mapping TransDecoder file:", td_gff_path))
  td_gr <- rtracklayer::import.gff(td_gff_path)
  # Filter for mRNA, CDS, and UTR features
  features_to_map <- td_gr[td_gr$type %in% c("mRNA", "CDS", "five_prime_UTR", "three_prime_UTR")]
  if (length(features_to_map) == 0) {
    return(GRanges()) # Return empty object if no features found
  }
  # The seqnames in TransDecoder's GFF are the transcript IDs
  names(features_to_map) <- seqnames(features_to_map)
  # Get exon structures from the main GTF, named by transcript ID
  exons_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
  # Keep only the transcripts that have ORF predictions to map
  exons_for_mapping <- exons_by_tx[names(exons_by_tx) %in% names(features_to_map)]
  # Map the transcript-relative coordinates to genomic coordinates
  mapped_features <- mapFromTranscripts(features_to_map, exons_for_mapping)
  # Create a lookup table from the original genomic GTF (transcript -> gene, oId)
  tx_to_gene_map <- as.data.frame(mcols(genomic_gr[genomic_gr$type == 'transcript',])) %>%
    dplyr::select(transcript_id, gene_id, oId) %>%
    dplyr::distinct()
  # Get the transcript_id for each mapped feature
  # The name of the original feature is the transcript_id
  original_feature_names <- names(features_to_map[mapped_features$xHits])
  mapped_features$transcript_id <- original_feature_names
  # Preserve the original TransDecoder info
  mapped_features$ID <- features_to_map[mapped_features$xHits]$ID
  mapped_features$Parent <- features_to_map[mapped_features$xHits]$Parent
  mapped_features$type <- features_to_map[mapped_features$xHits]$type
  # Convert to data.frame to merge in the gene_id and oId from our lookup table
  names(mapped_features)<- NULL
  mapped_features_df <- as.data.frame(mapped_features)
  mapped_features_df <- dplyr::left_join(mapped_features_df, tx_to_gene_map, by = "transcript_id")
  # Convert back to GRanges, which preserves all our new columns
  final_mapped_gr <- makeGRangesFromDataFrame(mapped_features_df, keep.extra.columns = TRUE)
  # Add source and clean up
  mcols(final_mapped_gr)$source <- "TransDecoder"
  # The original feature type is preserved from the TransDecoder GFF
  final_mapped_gr <- final_mapped_gr[
    with(final_mapped_gr, order(seqnames, start)),
  ]
  return(final_mapped_gr)
}

# Create a TxDb object from our genomic GTF for mapping
txdb <- makeTxDbFromGRanges(genomic_gr)
all_mapped_cds <- GRangesList()

# Map nuclear and mitochondrial ORFs if they exist
if (transdecoder_gff_nucl != "" && file.exists(transdecoder_gff_nucl)) {
  all_mapped_cds$nuclear <- map_transdecoder_to_genome(transdecoder_gff_nucl, txdb, genomic_gr)
}
if (transdecoder_gff_mito != "" && file.exists(transdecoder_gff_mito)) {
  all_mapped_cds$mitochondrial <- map_transdecoder_to_genome(transdecoder_gff_mito, txdb, genomic_gr)
}

# Unlist into a single GRanges object
final_mapped_cds <- unlist(all_mapped_cds)

# Combine the original genomic features with the newly mapped CDS features
# Remove any old/incorrect CDS features from the original GTF
stringtie_anno <- genomic_gr[genomic_gr$type != "CDS"]
# Append correctly mapped CDS features
stringtie_anno <- c(stringtie_anno, final_mapped_cds)
stringtie_anno <- stringtie_anno[
  with(stringtie_anno, order(seqnames, start)),
]
names(stringtie_anno)<- NULL
message("TransDecoder ORF mapping complete.")
# ---------------------------
# Conditionally import reference annotations
# ---------------------------
# Initialize ref_anno as NULL
ref_anno <- NULL

# Flag to indicate if reference annotation is available
ref_available <- FALSE

# Check if genome_ref argument is provided
if (genome_ref != "") {
  # genome_ref is expected to be the full path to the GTF file
  ref_anno_path <- genome_ref
  
  # Validate that the provided genome_ref is a file and has a .gtf extension
  if (file.exists(ref_anno_path) && grepl("\\.gtf$", ref_anno_path, ignore.case = TRUE)) {
    ref_available <- TRUE
    # Import reference annotations
    ref_anno <- rtracklayer::import.gff(ref_anno_path)
    message(paste("Reference annotation loaded from:", ref_anno_path))
  } else {
    warning(paste("Warning: Reference GTF file not found or invalid at", ref_anno_path, "- Skipping reference-based steps."))
  }
} else {
  warning("Warning: No --genome_ref provided - Skipping reference-based steps.")
}

# ---------------------------
# Import functional annotations
# ---------------------------
orfs_path <- file.path(functional_dir, "transcripts_with_orfs.txt")
blastp_path <- file.path(functional_dir, "blastp_results.out")
blastx_path <- file.path(functional_dir, "blastx_results.out")
pfam_path <- file.path(functional_dir, "pfam_results.out")
interproscan_path <- file.path(functional_dir, "interpro.tsv")
#pfam2go_path <- file.path(functional_dir, "pfam2go.txt")  # Ensure this file exists

if (!file.exists(orfs_path) & !file.exists(blastp_path) & !file.exists(blastx_path) & !file.exists(pfam_path) & !file.exists(interproscan_path)) {
  print("Functional annotation files are missing. Skipping integration.")
}

any_functional_data_exists <- FALSE

# Define empty data frame structures
empty_orfs <- data.frame(V1=character(), V2=character(), V3=character(), V4=character(), stringsAsFactors=FALSE)
empty_blast <- data.frame(transcript_id=character(), hit_id=character(), percent_identity=numeric(), 
                          alignment_length=integer(), mismatches=integer(), gap_opens=integer(), 
                          query_start=integer(), query_end=integer(), subject_start=integer(), 
                          subject_end=integer(), evalue=numeric(), bit_score=numeric(), stringsAsFactors=FALSE)
empty_pfam <- data.frame(domain_name=character(), domain_accession=character(), query_name=character(), 
                         query_accession=character(), sequence_evalue=numeric(), sequence_score=numeric(), 
                         sequence_bias=numeric(), domain_cevalue=numeric(), domain_ievalue=numeric(), 
                         domain_score=numeric(), domain_bias=numeric(), hmm_from=integer(), hmm_to=integer(), 
                         ali_from=integer(), ali_to=integer(), env_from=integer(), env_to=integer(), 
                         acc=numeric(), description=character(), domain_len=integer(), stringsAsFactors=FALSE)
empty_interpro <- data.frame(transcript_id=character(), MD5=character(), length=integer(), analysis=character(),
                             signature_accession=character(), description=character(), query_start=integer(), 
                             query_end=integer(), evalue=numeric(), status=character(), date=character(),
                             interpro_accession=character(), interpro_description=character(), GO=character(),
                             pathways=character(), stringsAsFactors=FALSE)

# Load functional annotations
# Load ORFs if file exists
orfs <- if (file.exists(orfs_path)) {
  message("Loading ORFs data.")
  any_functional_data_exists <<- TRUE
  read.table(orfs_path, stringsAsFactors = FALSE)
} else {
  empty_orfs
}

# Load BLASTp if file exists
uniprot_blastp <- if (file.exists(blastp_path)) {
  message("Loading BLASTp data.")
  any_functional_data_exists <<- TRUE
  read.table(blastp_path, stringsAsFactors = FALSE)
} else {
  empty_blast
}

# Load BLASTx if file exists
uniprot_blastx <- if (file.exists(blastx_path)) {
  message("Loading BLASTx data.")
  any_functional_data_exists <<- TRUE
  read.table(blastx_path, stringsAsFactors = FALSE)
} else {
  empty_blast
}

# Load Pfam if file exists
pfam <- if (file.exists(pfam_path)) {
  message("Loading Pfam data.")
  any_functional_data_exists <<- TRUE
  read_domtblout(pfam_path)
} else {
  empty_pfam
}

# Load InterProScan if file exists
interproscan <- if (file.exists(interproscan_path)) {
  message("Loading InterProScan data.")
  any_functional_data_exists <<- TRUE
  read.table(interproscan_path, stringsAsFactors = FALSE, sep = "\t", header = FALSE, quote = "")
} else {
  empty_interpro
}
#pfam2go <- read.table(pfam2go_path, sep = ";", stringsAsFactors = FALSE)


################################################################################
#---------------Step 10: Integrate functional annotation into gtf---------------
################################################################################

#------------------------FILTER DATA--------------------------------------------
# Prepare transcripts
transcripts <- as.data.frame(stringtie_anno)
# transcripts <- transcripts[,c("seqnames","start","end","width","strand","source","type","score","phase",
#                               "cmp_ref","gene_id","gene_name","ref_gene_id","transcript_id","exon_number","cmp_ref_gene" )]

gene_trans <- na.omit(unique(transcripts[,c("gene_id","transcript_id")]))

if (any_functional_data_exists) {
  
  message("Step 10: Integrating functional annotations...")
  
  #------------------------FILTER DATA--------------------------------------------
  # Process ORFs data
  if (nrow(orfs) > 0) {
    orfs <- orfs %>%
      tidyr::separate(V4, into = c("transcript_id", "orf_info"), sep = ":") %>%
      dplyr::select(transcript_id, orf_info)
  }
  
  # Process BLASTp results
  if (nrow(uniprot_blastp) > 0) {
    colnames(uniprot_blastp) <- c("transcript_id", "hit_id", "percent_identity", "alignment_length",
                                  "mismatches", "gap_opens", "query_start", "query_end",
                                  "subject_start", "subject_end", "evalue", "bit_score")
    uniprot_blastp <- uniprot_blastp %>%
      mutate(transcript_id = stringr::str_replace(transcript_id, "^((?:[^.]+\\.){2}[^.]+).*", "\\1"))
  }
  
  # Process BLASTx results
  if (nrow(uniprot_blastx) > 0) {
    colnames(uniprot_blastx) <- c("transcript_id", "hit_id", "percent_identity", "alignment_length",
                                  "mismatches", "gap_opens", "query_start", "query_end",
                                  "subject_start", "subject_end", "evalue", "bit_score")
  }
  
  # Process Pfam results
  pfam_df <- as.data.frame(pfam)
  if (nrow(pfam_df) > 0) {
    pfam_df <- pfam_df %>%
      mutate(transcript_id = stringr::str_replace(query_name, "^((?:[^.]+\\.){2}[^.]+).*", "\\1"))
  }
  
  # Process Interproscan results
  if (nrow(interproscan) > 0) {
    colnames(interproscan) <- c("transcript_id", "MD5", "length", "analysis",
                                "signature_accession", "description", "query_start", "query_end",
                                "evalue", "status", "date", "interpro_accession","interpro_description","GO","pathways")
    interproscan <- interproscan %>%
      mutate(transcript_id = stringr::str_replace(transcript_id, "^((?:[^.]+\\.){2}[^.]+).*", "\\1"))
  }
  
  #Filter Pfam data
  pfam_filtered <- data.frame()
  if (nrow(pfam_df) > 0) {
    pfam_df <- pfam_df %>%
      dplyr::mutate(coverage = (ali_to - ali_from + 1) / domain_len * 100)
    #cevalue_threshold <- 1e-5
    #bit_score_threshold <- 10
    #coverage_threshold <- 50
    
    pfam_filtered <- pfam_df %>%
      dplyr::filter(domain_cevalue < cevalue_threshold, domain_score > bit_score_threshold, coverage > coverage_threshold)
  }
  
  # Filter BLASTp hits
  signif_blastp <- data.frame()
  if (nrow(uniprot_blastp) > 0) {
    #blast_evalue_threshold <- 1e-5
    #blast_identity_threshold <- 25
    
    signif_blastp <- uniprot_blastp %>%
      dplyr::filter(evalue < blast_evalue_threshold & percent_identity > blast_identity_threshold) %>%
      mutate(
        UniProt_accession = stringr::str_extract(hit_id, "(?<=\\|)[A-Z0-9]+(?=\\|)"),
        UniProt_entry_name = stringr::str_extract(hit_id, "(?<=\\|)[A-Z0-9_]+$")
      )
  }
  
  # Filter BLASTx hits
  signif_blastx <- data.frame()
  if (nrow(uniprot_blastx) > 0) {
    #blast_evalue_threshold <- 1e-5
    #blast_identity_threshold <- 25
    
    signif_blastx <- uniprot_blastx %>%
      dplyr::filter(evalue < blast_evalue_threshold & percent_identity > blast_identity_threshold) %>%
      mutate(
        UniProt_accession = stringr::str_extract(hit_id, "(?<=\\|)[A-Z0-9]+(?=\\|)"),
        UniProt_entry_name = stringr::str_extract(hit_id, "(?<=\\|)[A-Z0-9_]+$")
      )
  }
  
  # Filter Interproscan hits
  interproscan_filtered <- data.frame()
  if (nrow(interproscan) > 0) {
    interproscan_filtered <- interproscan %>%
      dplyr::filter(interpro_accession != "-" & !is.na(interpro_accession)) %>%
      dplyr::filter(evalue != "-" & !is.na(evalue) & as.numeric(evalue) < interpro_evalue_threshold)
    # Parse the pathway column 
    # This separates each Database:ID pair into its own row
    pathway_map <- interproscan_filtered %>%
      dplyr::select(transcript_id, pathways) %>%
      tidyr::separate_rows(pathways, sep = "\\|") %>%
      dplyr::filter(pathways != "" & !is.na(pathways)) %>%
      # Separate "Database" from "ID"
      tidyr::separate(pathways, into = c("db", "id"), sep = ":", extra = "merge")
    
    print("Parsed Pathway IDs:")
    print(pathway_map)
    
    # Query Reactome for pathway names
    # Filter for only the Reactome entries
    reactome_ids <- pathway_map %>%
      dplyr::filter(db == "Reactome") %>%
      dplyr::pull(id) %>%
      unique()
    
    if(length(reactome_ids) > 0) {
      reactome_names <- AnnotationDbi::mapIds(
        reactome.db,
        keys = reactome_ids,
        column = "PATHNAME",
        keytype = "PATHID"
      )
      
      # Create a mapping table for the results
      reactome_results <- data.frame(
        id = paste0("R-HSA-", names(reactome_names)),
        pathway_name = reactome_names,
        stringsAsFactors = FALSE
      )
      
      print("Fetched Reactome Names:")
      print(reactome_results)
    }
    
    reactome_results$id <- row.names(reactome_results)
    reactome_results <- merge(pathway_map, reactome_results, by = "id", all.x = TRUE)
    reactome_results <- unique(reactome_results)
    reactome_results <- reactome_results %>%
      dplyr::mutate(
        # Create a new, clean column without the organism prefix
        pathway_name_clean = stringr::str_remove(pathway_name, "^.*: ")
      )
    reactome_results <- reactome_results[,-4]
    reactome_results <- unique(reactome_results)
    
    
    reactome_results_combined <- reactome_results %>%
      dplyr::filter(!is.na(pathway_name_clean) & pathway_name_clean != "") %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::summarize(all_pathways = paste(unique(pathway_name_clean), collapse = ", "), .groups = 'drop') %>% 
      dplyr::left_join(gene_trans, by = "transcript_id")
  }
  
  # Get best hits
  best_blastp <- if (nrow(signif_blastp) > 0) {
    signif_blastp %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(desc(percent_identity), evalue) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(transcript_id, best_SwissProt_blastp_hit = UniProt_entry_name)
  } else { data.frame(transcript_id=character(), best_SwissProt_blastp_hit=character()) }
  
  best_blastx <- if (nrow(signif_blastx) > 0) {
    signif_blastx %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(desc(percent_identity), evalue) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(transcript_id, best_SwissProt_blastx_hit = UniProt_entry_name)
  } else { data.frame(transcript_id=character(), best_SwissProt_blastx_hit=character()) }
  
  best_pfam <- if (nrow(pfam_filtered) > 0) {
    pfam_filtered %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(desc(coverage), domain_cevalue, desc(domain_score)) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(transcript_id, best_Pfam_hit = description)
  } else { data.frame(transcript_id=character(), best_Pfam_hit=character()) }
  
  best_interpro <- if (nrow(interproscan_filtered) > 0) {
    interproscan_filtered %>%
      dplyr::group_by(transcript_id) %>%
      dplyr::arrange(evalue) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup() %>%
      dplyr::select(transcript_id, best_interpro_hit = interpro_description)
  } else { data.frame(transcript_id=character(), best_interpro_hit=character()) }
  
  # Annotate transcripts with ORF information
  transcripts$has_ORF <- transcripts$transcript_id %in% orfs$transcript_id
  transcripts$has_ORF <- ifelse(transcripts$type=="gene", NA, transcripts$has_ORF)
  
  # Merge best hits with transcripts
  transcripts <- transcripts %>%
    dplyr::left_join(best_blastp, by = "transcript_id") %>%
    dplyr::left_join(best_blastx, by = "transcript_id") %>%
    dplyr::left_join(best_pfam, by = "transcript_id")  %>% 
    dplyr::left_join(best_interpro, by = "transcript_id")
  
  # Combine ALL hits into one dataframe
  pfam_hits <- if (nrow(pfam_filtered) > 0) pfam_filtered %>% dplyr::select(transcript_id, hit = description) else data.frame(transcript_id=character(), hit=character())
  blastx_hits <- if (nrow(signif_blastx) > 0) signif_blastx %>% dplyr::select(transcript_id, hit = UniProt_entry_name) else data.frame(transcript_id=character(), hit=character())
  blastp_hits <- if (nrow(signif_blastp) > 0) signif_blastp %>% dplyr::select(transcript_id, hit = UniProt_entry_name) else data.frame(transcript_id=character(), hit=character())
  interpro_hits <- if (nrow(interproscan_filtered) > 0) interproscan_filtered %>% dplyr::select(transcript_id, hit = interpro_description) else data.frame(transcript_id=character(), hit=character())
  interpro_GOs <- if (nrow(interproscan_filtered) > 0) interproscan_filtered %>% dplyr::select(transcript_id, GO = GO) else data.frame(transcript_id=character(), GO=character()) 
  interpro_GOs <- interpro_GOs[interpro_GOs$GO != "-",] %>% dplyr::group_by(transcript_id) %>%
    dplyr::summarize(all_GOs = paste(unique(GO), collapse = ", "), .groups = 'drop') %>% dplyr::left_join(gene_trans, by = "transcript_id")
  #interpro_pathways <- if (nrow(interproscan_filtered) > 0) interproscan_filtered %>% dplyr::select(transcript_id, hit = pathways) else data.frame(transcript_id=character(), hit=character())
  
  # MODIFIED: bind_rows now includes interpro_hits
  all_hits_combined <- dplyr::bind_rows(pfam_hits, blastx_hits, blastp_hits, interpro_hits) %>%
    dplyr::filter(!is.na(hit) & hit != "") %>%
    dplyr::group_by(transcript_id) %>%
    dplyr::summarize(all_hits = paste(unique(hit), collapse = ", "), .groups = 'drop') %>%
    dplyr::left_join(interpro_GOs[,-3], by = "transcript_id") %>%
    dplyr::left_join(reactome_results_combined[,-3], by = "transcript_id") 
    #dplyr::summarize(all_pathways = paste(unique(all_hits), collapse = ", "), .groups = 'drop')
  
  
  # Merge the all_hits_combined with transcripts
  transcripts_with_hits <- transcripts %>%
    dplyr::left_join(all_hits_combined, by = "transcript_id")
    
  interpro_GOs_genes <- interpro_GOs %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarize(all_GOs = paste(unique(all_GOs), collapse = ", "), .groups = 'drop')
  
  reactome_results_combined_genes <- reactome_results_combined %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarize(all_pathways = paste(unique(all_pathways), collapse = ", "), .groups = 'drop')
  
  # Process hits at the gene level
  #gene_trans <- na.omit(unique(transcripts[,c("gene_id","transcript_id")]))
  all_hits_combined_genes <- dplyr::bind_rows(pfam_hits, blastx_hits, blastp_hits, interpro_hits) %>%
    dplyr::left_join(gene_trans, by = "transcript_id") %>%
    dplyr::filter(!is.na(gene_id)) %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarize(all_hits = paste(unique(hit), collapse = ", "), .groups = 'drop') %>%
    dplyr::left_join(interpro_GOs_genes, by = "gene_id") %>%
    dplyr::left_join(reactome_results_combined_genes, by = "gene_id") 
  
  # Filter out transcripts with no functional information
  functional_transcripts <- transcripts_with_hits
  #   filter(type == "gene" | (type != "gene" & (has_ORF == TRUE | !is.na(best_SwissProt_blastp_hit) | !is.na(best_SwissProt_blastx_hit) | !is.na(best_Pfam_hit))))
  
  genes_tb <- as.data.frame(table(functional_transcripts$gene_id))
  functional_transcripts <- functional_transcripts[functional_transcripts$gene_id %in% genes_tb$Var1[genes_tb$Freq > 1],]
  
  # Aggregate functional annotations at the gene level
  functional_genes_unique <- functional_transcripts %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarize(
      best_SwissProt_blastp_hit = paste(unique(na.omit(best_SwissProt_blastp_hit)), collapse = ", "),
      best_SwissProt_blastx_hit = paste(unique(na.omit(best_SwissProt_blastx_hit)), collapse = ", "),
      best_Pfam_hit = paste(unique(na.omit(best_Pfam_hit)), collapse = ", "),
      best_InterPro_hit = paste(unique(na.omit(best_interpro_hit)), collapse = ", "),
      has_ORF = any(has_ORF, na.rm = TRUE),
      all_GOs = paste(unique(na.omit(all_GOs)), collapse = ", "),
      all_pathways = paste(unique(na.omit(all_pathways)), collapse = ", "),
      .groups = 'drop'
    )
  
  # Update gene-level annotations
  functional_genes_annotated <- functional_transcripts %>%
    dplyr::filter(type == "gene") %>%
    dplyr::select(-any_of(c("best_SwissProt_blastp_hit", "best_SwissProt_blastx_hit", "best_Pfam_hit", "best_InterPro_hit","all_hits","has_ORF","all_GOs","all_pathways"))) %>%
    dplyr::left_join(functional_genes_unique, by = "gene_id")
  
  functional_non_genes_subset <- functional_transcripts %>%
    dplyr::filter(type != "gene")
  
  functional_transcripts_updated <<- dplyr::bind_rows(functional_genes_annotated, functional_non_genes_subset) %>%
    dplyr::arrange(match(gene_id, functional_transcripts$gene_id))
  
  # Merge gene-level all_hits into the final data frame
  functional_transcripts_updated <<- functional_transcripts_updated %>%
    dplyr::left_join(all_hits_combined_genes, by = "gene_id", suffix = c("", "_gene")) %>%
    dplyr::mutate(all_hits = if_else(type == "gene", all_hits_gene, all_hits)) %>%
    dplyr::mutate(all_GOs = if_else(type == "gene", all_GOs_gene, all_GOs)) %>%
    dplyr::mutate(all_pathways = if_else(type == "gene", all_pathways_gene, all_pathways)) %>%
    dplyr::select(-c(all_hits_gene, all_GOs_gene, all_pathways_gene))
  
} else {
  message("Step 10: No functional annotation files found - Skipping functional annotation integration.")
  # If skipping, ensure the main data frame has the necessary columns for downstream steps, filled with NA
  functional_transcripts_updated <- transcripts
  functional_transcripts_updated$best_SwissProt_blastp_hit <- NA_character_
  functional_transcripts_updated$best_SwissProt_blastx_hit <- NA_character_
  functional_transcripts_updated$best_Pfam_hit <- NA_character_
  functional_transcripts_updated$best_InterPro_hit <- NA_character_
  functional_transcripts_updated$all_hits <- NA_character_
  functional_transcripts_updated$has_ORF <- NA
  functional_transcripts_updated$all_GOs <- NA_character_
  functional_transcripts_updated$all_pathways <- NA_character_
}


################################################################################
#---------------Step 11: Find overlapping genes and transcripts-----------------
################################################################################
# Convert functional_transcripts to GRanges
message("Step 11: Finding overlapping genes and transcripts...")
functional_transcripts_gr <- GenomicRanges::makeGRangesFromDataFrame(
  functional_transcripts_updated,
  keep.extra.columns = TRUE,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)


# Subset to genes only
functional_genes_gr <- functional_transcripts_gr[functional_transcripts_gr$type == "gene"]
# Assign gene_id as names of the GRanges object
names(functional_genes_gr) <- functional_genes_gr$gene_id

# Convert GRanges to data.table
functional_genes_dt <- as.data.table(as.data.frame(functional_genes_gr))
# Select necessary columns
functional_genes_dt <- functional_genes_dt[, .(seqnames, start, end, gene_id)]
# Set keys for fast overlap queries
setkey(functional_genes_dt, seqnames, start, end)
# Perform overlaps using data.table's foverlaps
overlaps_dt <- foverlaps(functional_genes_dt, functional_genes_dt, nomatch = NULL)
# Remove self-overlaps if any
overlaps_dt <- overlaps_dt[gene_id != i.gene_id]
overlaps_dt <- overlaps_dt[, c("gene_id", "i.gene_id")]
colnames(overlaps_dt) <- c("gene_id", "overlapped_predicted_gene")
overlaps_dt_uniq <- overlaps_dt %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarize(
    overlapped_predicted_gene = {
      hits <- unique(na.omit(overlapped_predicted_gene))
      if(length(hits) > 0) paste(hits, collapse = "; ") else NA_character_
    },
    .groups = 'drop' # To ungroup after summarization
  )

# Merge overlaps with functional_transcripts_updated
functional_transcripts_updated <- functional_transcripts_updated %>%
  dplyr::left_join(
    overlaps_dt_uniq,
    by = "gene_id"
  )

# Ensure that every transcript within a gene has at least one overlapping exonic region with any other transcript in the same gene
check_all_transcripts_overlap <- function(gene_exons_gr) {
  # Get unique transcript IDs within the gene
  transcripts <- unique(mcols(gene_exons_gr)$transcript_id)
  # If only one transcript, no overlaps to check
  if (length(transcripts) <= 1) {
    return(NA)
  }
  # Create a list of GRanges for each transcript
  transcript_list <- split(gene_exons_gr, mcols(gene_exons_gr)$transcript_id)
  # Convert to GRangesList
  transcript_grl <- GRangesList(transcript_list)
  # Find overlaps within the GRangesList
  overlaps <- rtracklayer::findOverlaps(transcript_grl, transcript_grl)
  # Remove self-overlaps
  overlaps <- overlaps[queryHits(overlaps) != subjectHits(overlaps)]
  # Create a list where each transcript maps to the transcripts it overlaps with
  overlap_map <- split(subjectHits(overlaps), queryHits(overlaps))
  # Check for each transcript if it has at least one overlap
  transcripts_with_overlap <- names(overlap_map)
  # All transcripts that have at least one overlap
  transcripts_flagged <- transcripts_with_overlap
  # Compare with all transcripts to see if any lack overlaps
  all_transcripts_flagged <- transcripts_flagged %in% transcripts
  # If all transcripts have at least one overlap, return TRUE
  if (length(transcripts_flagged) == length(transcripts)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# Convert to GRanges 
if (!inherits(functional_transcripts_updated, "GRanges")) {
  functional_transcripts_updated <- GenomicRanges::makeGRangesFromDataFrame(
    functional_transcripts_updated,
    keep.extra.columns = TRUE,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand"
  )
}

# Subset to exon features if available
exons_gr <- functional_transcripts_updated[functional_transcripts_updated$type == "exon"]
# Subset to gene features
genes_gr <- functional_transcripts_updated[functional_transcripts_updated$type == "gene"]

# Extract gene IDs
gene_ids <- unique(mcols(genes_gr)$gene_id)

# Initialize a dataframe to store overlap information
transcript_overlap_df <- data.frame(
  gene_id = character(),
  has_all_transcripts_overlapping = logical(),
  stringsAsFactors = FALSE
)

# Iterate over each gene to check for overlapping transcripts
for (gene in gene_ids) {
  # Subset exons for the current gene
  gene_exons <- exons_gr[mcols(exons_gr)$gene_id == gene]
  # Check for overlaps
  overlaps_flag <- check_all_transcripts_overlap(gene_exons)
  # Append to the dataframe
  transcript_overlap_df <- rbind(
    transcript_overlap_df,
    data.frame(
      gene_id = gene,
      has_all_transcripts_overlapping  = overlaps_flag,
      stringsAsFactors = FALSE
    )
  )
}

if (ref_available) {  # Only perform if reference annotation is available
  #------Detect genes that share genomic regions, indicating potential overlaps---
  # Interpretation of columns added by gffcompare: ref_gene_id - unique identifier of the matched reference gene, predicted gene has been mapped to an exact reference gene ID
  # cmp_ref_gene - gene name or ID in the reference annotation that the predicted gene overlaps with or aligns to
  
  
  ref_anno_gr <- GenomicRanges::makeGRangesFromDataFrame(
    ref_anno,
    keep.extra.columns = TRUE,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand"
  )
  
  # Subset to genes only
  ref_anno_gr <- ref_anno_gr[ref_anno_gr$type == "gene"]
  ref_anno_dt <- as.data.table(as.data.frame(ref_anno_gr))
  # Select necessary columns
  ref_anno_dt <- ref_anno_dt[, .(seqnames, start, end, gene_id)]
  # Set keys for fast overlap queries
  setkey(ref_anno_dt, seqnames, start, end)
  # Perform overlaps using data.table's foverlaps
  overlaps_with_ref_dt <- foverlaps(functional_genes_dt, ref_anno_dt, nomatch = NULL)
  # Remove self-overlaps if any
  overlaps_with_ref_dt <- overlaps_with_ref_dt[gene_id != i.gene_id]
  overlaps_with_ref_dt <- overlaps_with_ref_dt[, c("gene_id", "i.gene_id")]
  colnames(overlaps_with_ref_dt) <- c("overlapped_ref_gene", "gene_id")
  
  overlaps_with_ref_dt_uniq <- overlaps_with_ref_dt %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarize(
      overlapped_ref_gene = {
        hits <- unique(na.omit(overlapped_ref_gene))
        if(length(hits) > 0) paste(hits, collapse = "; ") else NA_character_
      },
      .groups = 'drop' # To ungroup after summarization
    )
  
  # Merge overlaps with functional_transcripts_updated
  functional_transcripts_updated <- functional_transcripts_updated %>%
    dplyr::left_join(
      overlaps_with_ref_dt_uniq,
      by = "gene_id"
    )
  
} else {
  # If reference annotation is not available, partially skip Step 11
  message("Step 11: Reference annotation not available - Skipping overlap detection with reference annotation.")
}
# Merge overlap information back to functional_transcripts_updated
functional_transcripts_updated <- as.data.frame(functional_transcripts_updated)
functional_transcripts_updated <- functional_transcripts_updated %>%
  dplyr::left_join(transcript_overlap_df, by = "gene_id")
################################################################################
#-----------Step 12: Find reversed duplicates-----------------------------------
################################################################################
if (ref_available) {  # Only perform if reference annotation is available
  message("Step 12: Finding reversed duplicates...")
  ref_anno_dt <- as.data.table(as.data.frame(ref_anno_gr))
  ref_anno_dt <- ref_anno_dt[, .(gene_id, strand_ref = as.character(strand))]
  head(ref_anno_dt)
  
  functional_transcripts_updated <- functional_transcripts_updated %>%
    dplyr::left_join(ref_anno_dt, by = c("overlapped_ref_gene" = "gene_id"))
  
  functional_transcripts_updated <- functional_transcripts_updated %>%
    dplyr::group_by(overlapped_ref_gene) %>%
    dplyr::mutate(
      # Flag if the predicted gene's strand matches the reference gene's strand
      same_strand = (strand == strand_ref),
      # Determine if any predicted genes in the group are on the same strand
      has_same_strand = any(same_strand, na.rm = TRUE)
    ) %>%
    # Retain only predicted genes on the same strand if they exist; otherwise, keep those on the opposite strand
    dplyr::filter((has_same_strand & same_strand) | (!has_same_strand)) %>%
    dplyr::ungroup()
  
} else {
  # If reference annotation is not available, skip Step 12
  message("Step 12: Reference annotation not available - Skipping reversed duplicates analysis.")
}
################################################################################
#-----------Step 13.1: Find fragmented and chimeric genes (with reference)--------
################################################################################
if (ref_available) {  # Only perform if reference annotation is available
  # One to many: one reference gene matches to many predictions
  message("Step 13.1: Finding fragmented and chimeric genes (with reference)...")
  fragmented_genes <- functional_transcripts_updated[functional_transcripts_updated$type == "gene", ] %>%
    dplyr::group_by(overlapped_ref_gene) %>%
    dplyr::filter(n() > 1) %>%
    dplyr::distinct(overlapped_ref_gene) %>%
    dplyr::pull(overlapped_ref_gene)
  
  functional_transcripts_updated$is_potentially_fragmented_reference <- functional_transcripts_updated$overlapped_ref_gene %in% fragmented_genes
  
  # Many to one: many reference genes match to one predicted gene
  chimeric_genes <- functional_transcripts_updated[functional_transcripts_updated$type == "gene", ] %>%
    dplyr::group_by(gene_id) %>%
    dplyr::summarize(n_ref_genes = n_distinct(overlapped_ref_gene), .groups = 'drop') %>%
    dplyr::filter(n_ref_genes > 1) %>%
    dplyr::pull(gene_id)
  
  functional_transcripts_updated$is_potentially_chimeric_reference <- functional_transcripts_updated$overlapped_ref_gene %in% chimeric_genes
  
} else {
  # If reference annotation is not available, skip Step 13.1
  message("Step 13.1: Reference annotation not available - Skipping fragmented and chimeric genes analysis with reference.")
}
################################################################################
#-----------Step 13.2: Find fragmented and chimeric genes (de novo)-------------
################################################################################
message("Step 13.2: Finding fragmented and chimeric genes (de novo)...")
# ----------------------------------------
# Identify fragmented genes de novo
# ----------------------------------------
# Subset to gene entries
#head(functional_transcripts_updated)
genes_df <- functional_transcripts_updated %>%
  dplyr::filter(type == "gene")

# Convert to GRanges
genes_gr <- GenomicRanges::makeGRangesFromDataFrame(
  genes_df,
  keep.extra.columns = TRUE,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)

# Sort the GRanges object
genes_gr_sorted <- sort(genes_gr)
# Define the maximum distance to consider genes as fragmented (e.g., 1000 bp)
#max_distance <- 1000
# Initialize cluster IDs
cluster_id <- 1
clusters <- rep(NA, length(genes_gr_sorted))
# Assign clusters based on proximity and strand
prev_gene <- NULL
for (i in seq_along(genes_gr_sorted)) {
  current_gene <- genes_gr_sorted[i]
  
  if (is.null(prev_gene)) {
    clusters[i] <- cluster_id
    prev_gene <- current_gene
  } else {
    # Check if current gene is on the same strand and within max_distance
    if (as.character(strand(current_gene)) == as.character(strand(prev_gene)) &&
        start(current_gene) - end(prev_gene) <= max_distance) {
      clusters[i] <- cluster_id
    } else {
      cluster_id <- cluster_id + 1
      clusters[i] <- cluster_id
    }
    prev_gene <- current_gene
  }
}

# Add cluster information to genes_gr_sorted
genes_gr_sorted$fragment_cluster <- clusters
# Convert back to data frame for further processing
genes_df_sorted <- as.data.frame(genes_gr_sorted)
# Identify clusters with more than one gene
fragmented_clusters <- genes_df_sorted %>%
  dplyr::group_by(fragment_cluster) %>%
  dplyr::filter(n() > 1) %>%
  dplyr::pull(fragment_cluster) %>%
  unique()

functional_annotations_clusters <- genes_df_sorted %>%
  dplyr::filter(fragment_cluster %in% fragmented_clusters) %>%
  dplyr::group_by(fragment_cluster, gene_id) %>%
  # Reframe to aggregate functional_entries
  dplyr::reframe(
    functional_entries = unique(
      na.omit(
        unlist(str_split(all_hits, ","), use.names = FALSE)
      )
    ),
    .groups = 'drop'
  ) %>%
  dplyr::group_by(fragment_cluster) %>%
  dplyr::reframe(
    functional_entries = list(functional_entries),
    .groups = 'drop'
  )
# Function to compute common functional entries within a cluster
compute_common_entries <- function(entries_list) {
  # Remove empty strings and NA
  entries_list_clean <- lapply(entries_list, function(x) x[x != "" & !is.na(x)])
  # Unlist the data to get a character vector
  functional_entries <- unlist(entries_list_clean)
  # Trim leading and trailing whitespace
  functional_entries <- stringr::str_trim(functional_entries)
  # Remove empty strings
  functional_entries <- functional_entries[functional_entries != ""]
  # Identify duplicates
  duplicate_flags <- duplicated(functional_entries)
  # Extract duplicated entries
  common <- functional_entries[duplicate_flags]
  # # Check if any gene has no functional entries
  # if(any(sapply(entries_list_clean, length) == 0)) {
  #   return(character(0))  # No common entries if any gene lacks annotations
  # }
  # # Compute intersection of functional entries across all genes in the cluster
  # common <- Reduce(intersect, entries_list_clean)
  # Remove any empty strings that might have slipped through
  common <- common[common != "" & !is.na(common)]
  return(common)
}

functional_annotations_clusters <- functional_annotations_clusters %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    common_entries = list(compute_common_entries(functional_entries))
  ) %>%
  dplyr::ungroup()

clusters_to_flag <- functional_annotations_clusters %>%
  dplyr::rowwise() %>%
  dplyr::filter(length(common_entries) > 0 ) %>%
  dplyr::pull(fragment_cluster)

# Flag genes in 'genes_df_sorted' accordingly
genes_df_sorted <- genes_df_sorted %>%
  dplyr::mutate(
    is_potentially_fragmented_de_novo = ifelse(fragment_cluster %in% clusters_to_flag, TRUE, FALSE)
  )

# Update the main dataframe with fragmentation information
functional_transcripts_updated <- functional_transcripts_updated %>%
  dplyr::left_join(
    genes_df_sorted %>% 
      dplyr::select(gene_id, is_potentially_fragmented_de_novo),
    by = "gene_id"
  )
# ----------------------------------------
# Identify chimeric genes based on structural anomalies
# ----------------------------------------

intron_info <- functional_transcripts_updated %>%
  dplyr::filter(type == "intron") %>%
  dplyr::arrange(gene_id, transcript_id) %>%
  dplyr::group_by(gene_id, transcript_id) %>%
  dplyr::ungroup()

# Define a threshold for unusually large introns (e.g., >100,000 bp)
#large_intron_threshold <- 100000

# Identify genes with large introns
genes_with_large_introns <- intron_info %>%
  dplyr::filter(!is.na(width) & width > large_intron_threshold) %>%
  dplyr::pull(gene_id) %>%
  unique()

# Flag genes with large introns as potentially chimeric
functional_transcripts_updated <- functional_transcripts_updated %>%
  dplyr::mutate(
    is_potentially_chimeric_structure_de_novo = gene_id %in% genes_with_large_introns
  )

# ----------------------------------------
# Identify chimeric genes based on functional annotations
# ----------------------------------------
# Identify if any transcript within a gene has functional annotations.
# Count the number of transcripts (transcripts_nr) per gene and identify functional entries that appear exactly transcripts_nr times, indicating they are present in all transcripts.
# Only flag a gene as chimeric if there are functional annotations present but no common entries across all transcripts.
# Prevent genes with overlapping transcripts (i.e., transcripts that share exons) from being flagged as chimeric
# Ensure that genes with all transcripts lacking annotations are not flagged as chimeric.

# Examine only genes with > 1 transcript

# Function to compute common entries across all transcripts
compute_common_entries_all <- function(entries_list) {
  # Clean each transcript's entries
  entries_list_clean <- lapply(entries_list, function(x) {
    # Remove NA and empty strings
    x <- x[!is.na(x) & x != ""]
    # Trim whitespace
    x <- stringr::tr_trim(x)
    # Return unique entries within the transcript
    unique(x)
  })
  # Check if there are any transcripts left after cleaning
  if(length(entries_list_clean) == 0) {
    return(character(0))  # No entries to compare
  }
  # Check if any transcripts have functional entries
  # any_non_empty <- any(sapply(entries_list_clean, length) > 0)
  # # If no transcripts have functional entries, return NULL or a special value
  # if(!any_non_empty){
  #   return(NULL)  # Indicates no functional annotations present
  # }
  # Count the number of transcripts
  transcripts_nr <- length(entries_list_clean)
  # Create a frequency table of all entries
  freq_table <- table(unlist(entries_list_clean))
  # Identify entries that appear exactly 'transcripts_nr' times
  common <- names(freq_table)[freq_table == transcripts_nr]
  #common <- Reduce(intersect, entries_list_clean)
  
  return(common)
}

# functional_annotations_transcripts <- functional_transcripts_updated %>%
#   dplyr::filter(type == "transcript") %>%
#   dplyr::group_by(gene_id) %>%
#   dplyr::filter(n() > 1) %>%
#   # Keep genes with >1 transcript
#   dplyr::ungroup()   %>%
#   dplyr::group_by(transcript_id, gene_id) %>%
#   dplyr::reframe(
#     functional_entries = pmap(
#       list(all_hits),
#       ~ {
#         combined_entries <- c(..1, ..2, ..3)
#         split_entries <- combined_entries %>%
#           str_split(",") %>%
#           unlist() %>%
#           str_trim()
#         clean_entries <- split_entries[!is.na(split_entries) & split_entries != ""]
#         unique(clean_entries)
#       }
#     )
#   )
# functional_transcripts_updated
# list(unlist(functional_transcripts_updated$all_hits))
# tidyr::separate_rows(functional_transcripts_updated$all_hits, sep = "\\|")
# unlist(str_split(functional_transcripts_updated$all_hits, ","), use.names = FALSE)


#unique(functional_transcripts_updated$gene_id)
# Filter transcripts with functional annotations
#head(functional_transcripts_updated)
functional_annotations_transcripts <- functional_transcripts_updated %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::filter(n() > 1) %>%
  # Keep genes with >1 transcript
  dplyr::ungroup()   %>%
  dplyr::group_by( gene_id) %>%
  dplyr::summarize(
    functional_entries = list(unlist(all_hits)),
    #functional_entries = unlist(str_split(functional_transcripts_updated$all_hits, ","), use.names = FALSE),
    transcripts_nr = n(),
    .groups = 'drop'
  )  %>%
  dplyr::mutate(
    common_entries = purrr::map(functional_entries, compute_common_entries_all),
    # Determine if any transcript has functional annotations
    has_any_annotations = purrr::map_lgl(functional_entries, ~ any(sapply(.x, function(fe) {
      fe_clean <- fe[!is.na(fe) & fe != ""]
      length(fe_clean) > 0
    }))),
    # Conditionally flag as chimeric based on common functional entries
    is_potentially_chimeric_function_de_novo = if_else(
      has_any_annotations,
      # No common entries across all transcripts
      (purrr::map_lgl(common_entries, ~ is.null(.x) || length(.x) == 0)),
      FALSE
    )
  ) %>%
  # Merge overlap information
  dplyr::left_join(transcript_overlap_df, by = "gene_id") %>%
  # Update the flag to exclude genes with overlapping transcripts
  dplyr::mutate(
    is_potentially_chimeric_function_de_novo = if_else(
      is_potentially_chimeric_function_de_novo & has_all_transcripts_overlapping == FALSE,
      TRUE,
      FALSE
    )
  )


#functional_annotations_transcripts[1,2][[1]]
# Update the main dataframe with fragmentation information
functional_transcripts_updated <- functional_transcripts_updated %>%
  dplyr::left_join(
    functional_annotations_transcripts %>% 
      dplyr::select(gene_id, is_potentially_chimeric_function_de_novo),
    by = "gene_id"
  )

# ---------------------------
# Final GTF Formatting and Export
# ---------------------------
# In functional_annotation.R

# ---------------------------
# Final GTF Formatting and Export
# ---------------------------

# Convert to data frame for final manipulation with dplyr
final_df <- as.data.frame(functional_transcripts_updated)

# --- Clean up feature types and generate informative IDs ---
final_df_cleaned <- final_df %>%
  # 1. Remove redundant 'transcript' rows with generic IDs
  dplyr::filter(!(type == "transcript" & grepl("^nbis-transcript-", ID))) %>%
  # 2. Standardize all transcript-level features ('RNA', 'mRNA') to 'transcript'
  dplyr::mutate(type = if_else(type %in% c("RNA", "mRNA"), "transcript", type)) %>%
  # 3. Generate clean, hierarchical IDs
  dplyr::group_by(transcript_id) %>%
  dplyr::mutate(
    ID = case_when(
      type == "gene"       ~ gene_id,
      type == "transcript" ~ transcript_id,
      type == "exon"       ~ paste0(transcript_id, ".exon.", exon_number),
      # Keep original IDs for CDS and UTR features, which are already hierarchical
      TRUE                 ~ ID 
    )
  ) %>%
  dplyr::ungroup()

# --- Remove functional annotations from all sub-features ---
annotation_cols <- c("best_SwissProt_blastp_hit", "best_SwissProt_blastx_hit", 
                     "best_Pfam_hit", "best_InterPro_hit", "all_hits", "has_ORF", 
                     "all_GOs", "all_pathways")
# Ensure columns exist before trying to modify them
cols_to_clean <- intersect(annotation_cols, names(final_df_cleaned))

sub_feature_types <- c("exon", "intron", "CDS", "five_prime_UTR", "three_prime_UTR")

final_df_cleaned <- final_df_cleaned %>%
  dplyr::mutate(across(all_of(cols_to_clean), ~if_else(type %in% sub_feature_types, NA, .)))


# --- Convert back to GRanges and Export ---
final_granges_for_export <- GenomicRanges::makeGRangesFromDataFrame(
  final_df_cleaned,
  keep.extra.columns = TRUE,
  seqnames.field = "seqnames",
  start.field = "start",
  end.field = "end",
  strand.field = "strand"
)

final_gtf_path <- file.path(functional_dir, "final_annotation.gtf")
rtracklayer::export(final_granges_for_export, final_gtf_path, format = "gtf")

cat("Functional annotation completed successfully. Final GTF saved at:", final_gtf_path, "\n")

# ---------------------------
# Print Summary and Recommendations
# ---------------------------
print_summary_and_advice <- function(final_df, ref_is_available) {
  
  # Ensure flag columns exist, if not, create them as FALSE
  flag_cols <- c("overlapped_predicted_gene", "is_potentially_fragmented_reference", 
                 "is_potentially_chimeric_reference", "is_potentially_fragmented_de_novo", 
                 "is_potentially_chimeric_structure_de_novo", "is_potentially_chimeric_function_de_novo")
  for(col in flag_cols){
    if(!col %in% names(final_df)){
      final_df[[col]] <- NA
    }
  }
  
  # Filter for gene-level entries only for counting
  gene_df <- final_df %>%
    dplyr::filter(type == "gene")
  
  total_genes <- dplyr::n_distinct(gene_df$gene_id)
  total_transcripts <- final_df %>%
    dplyr::filter(type == "transcript") %>%
    dplyr::summarise(n = dplyr::n_distinct(transcript_id)) %>%
    dplyr::pull(n)
  
  overlapping_genes <- dplyr::n_distinct(gene_df$gene_id[!is.na(gene_df$overlapped_predicted_gene)])
  fragmented_de_novo <- dplyr::n_distinct(gene_df$gene_id[which(gene_df$is_potentially_fragmented_de_novo == TRUE)])
  chimeric_structural <- dplyr::n_distinct(gene_df$gene_id[which(gene_df$is_potentially_chimeric_structure_de_novo == TRUE)])
  chimeric_functional <- dplyr::n_distinct(gene_df$gene_id[which(gene_df$is_potentially_chimeric_function_de_novo == TRUE)])
  
  # Reference-based stats
  fragmented_ref <- if(ref_is_available) dplyr::n_distinct(gene_df$gene_id[which(gene_df$is_potentially_fragmented_reference == TRUE)]) else "N/A"
  chimeric_ref <- if(ref_is_available) dplyr::n_distinct(gene_df$gene_id[which(gene_df$is_potentially_chimeric_reference == TRUE)]) else "N/A"
  reversed_duplicates_info <- if(ref_is_available) "Resolved where possible" else "N/A"
  
  cat("\n\n--- SmedAnno final annotation summary ---\n\n")
  cat(sprintf("Total genes annotated: %d\n", total_genes))
  cat(sprintf("Total transcripts annotated: %d\n\n", total_transcripts))
  
  cat("--- Quality control flags ---\n")
  cat(sprintf("Potentially overlapping genes: %s\n", overlapping_genes))
  cat(sprintf("Potentially fragmented genes (de novo): %s\n", fragmented_de_novo))
  cat(sprintf("Potentially chimeric genes (structural anomaly): %s\n", chimeric_structural))
  cat(sprintf("Potentially chimeric genes (functional conflict): %s\n\n", chimeric_functional))
  
  cat("--- Reference-based QC flags ---\n")
  cat(sprintf("Potentially fragmented genes (vs. reference): %s\n", fragmented_ref))
  cat(sprintf("Potentially chimeric genes (vs. reference): %s\n", chimeric_ref))
  cat(sprintf("Reversed duplicates: %s\n\n", reversed_duplicates_info))
  
  cat("--- Recommendations and next steps ---\n\n")
  cat("The automated annotation is complete, but manual inspection of flagged genes is critical for a high-quality final gene set.\n\n")
  
  cat("1. For POTENTIALLY FRAGMENTED and OVERLAPPING genes:\n")
  cat("   These genes are candidates for merging into a single, contiguous gene model. Analysis should include:\n")
  cat("   - Visual Inspection: Load  final GTF and the alignment BAM files into a genome browser (e.g., IGV).\n     Look for continuous read coverage spanning the gap between the flagged genes. This is strong evidence for merging.\n")
  cat("   - Promoter/Terminator Analysis: As noted in the manuscript, extract the genomic sequence between the flagged genes.\n     Use tools like Promoter 2.0 (for promoters) and ARNold (for terminators) to check for regulatory elements.\n     The ABSENCE of these elements in the intervening sequence supports a merge.\n\n")
  
  cat("2. For POTENTIALLY CHIMERIC genes:\n")
  cat("   These genes may be incorrect fusions of two or more distinct genes. Your analysis should include:\n")
  cat("   - Structural anomalies (Large introns): For genes flagged with large introns, inspect them in a genome browser.\n     A very large intron with poor or non-existent read coverage often indicates an incorrect fusion. Check if other valid genes exist within that 'intron'.\n")
  cat("   - Functional conflicts: This is a strong indicator of a chimeric gene. The gene model likely fuses exons from two different functional genes.\n     Identify the conflicting annotations (e.g., from the 'all_hits' column) and use a genome browser to locate the likely fusion point and manually split the gene model.\n\n")
  
  cat("3. For REFERENCE-BASED Flags (if applicable):\n")
  cat("   - 'Fragmented (vs. reference)': This occurs when multiple predicted genes map to a single reference gene. This is a high-priority case for merging. Use the reference gene as a guide for the correct merged structure.\n")
  cat("   - 'Chimeric (vs. reference)': This occurs when a single predicted gene spans multiple distinct reference genes. This is strong evidence of a fusion event in your assembly and the gene should be split.\n")
  cat("   - 'Reversed duplicates': The pipeline attempts to resolve these by prioritizing predictions on the correct strand. Any loci where only an opposite-strand gene remains should be carefully reviewed for potential genuine antisense transcripts or annotation errors.\n\n")
  
  cat("FINAL ADVICE: After performing manual reviews and edits, you should finalize gene models. This ensures the most accurate representation of gene structures, lengths, and functional annotations for downstream analyses.\n\n")
}

print_summary_and_advice(as.data.frame(functional_transcripts_updated), ref_available)
