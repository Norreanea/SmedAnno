#!/usr/bin/env Rscript

# functional_annotation.R
# Script to perform functional annotation and identify fragmented/chimeric genes

# ---------------------------
# Function to install and load Libraries
# ---------------------------
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      if (pkg %in% rownames(installed.packages())) {
        library(pkg, character.only = TRUE)
      } else {
        if (pkg %in% c("GenomicFeatures", "UniProt.ws")) {
          if (!requireNamespace("BiocManager", quietly = TRUE)) {
            install.packages("BiocManager", repos = "http://cran.us.r-project.org")
          }
          BiocManager::install(pkg, ask = FALSE)
        } else {
          install.packages(pkg, repos = "http://cran.us.r-project.org")
        }
        library(pkg, character.only = TRUE)
      }
    } else {
      library(pkg, character.only = TRUE)
    }
  }
}

# ---------------------------
# Define required libraries
# ---------------------------
required_packages <- c(
  "rtracklayer",
  "dplyr",
  "data.table",
  "stringr",
  "rhmmer",
  "tidyr",
  "GenomicFeatures",
  "UniProt.ws",
  "GenomicRanges"
)

# Install and load libraries
install_and_load(required_packages)

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

# Parse arguments
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
  }
}

# ---------------------------
# Validate arguments
# ---------------------------
# Make genome_ref optional by not requiring it
if (annotation_dir == "" || functional_dir == "" || output_dir == "") {
  stop("Error: Arguments --annotation_dir, --functional_dir, and --output_dir must be provided.")
}


# ---------------------------
# Import GTF annotations
# ---------------------------
stringtie_anno_path <- file.path(annotation_dir, "corrected_with_introns.gtf")

if (!file.exists(stringtie_anno_path)) {
  stop(paste("Error: GTF file not found at", stringtie_anno_path))
}
# Import StringTie annotations
stringtie_anno <- import.gff(stringtie_anno_path)
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
    ref_anno <- import.gff(ref_anno_path)
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
pfam2go_path <- file.path(functional_dir, "pfam2go.txt")  # Ensure this file exists

if (!file.exists(orfs_path) || !file.exists(blastp_path) || !file.exists(blastx_path) || !file.exists(pfam_path)) {
  stop("Error: One or more functional annotation files are missing.")
}

# Load functional annotations
orfs <- read.table(orfs_path, stringsAsFactors = FALSE)
uniprot_blastp <- read.table(blastp_path, stringsAsFactors = FALSE)
uniprot_blastx <- read.table(blastx_path, stringsAsFactors = FALSE)
pfam <- rhmmer::read_domtblout(pfam_path)
pfam2go <- read.table(pfam2go_path, sep = ";", stringsAsFactors = FALSE)


################################################################################
#---------------Step 10: Integrate functional annotation into gtf---------------
################################################################################

#------------------------FILTER DATA--------------------------------------------
# Prepare transcripts
transcripts <- as.data.frame(stringtie_anno)
transcripts <- transcripts[,c("seqnames","start","end","width","strand","source","type","score","phase",
                              "cmp_ref","gene_id","gene_name","ref_gene_id","transcript_id","exon_number","cmp_ref_gene" )]
# Process ORFs data
orfs <- orfs %>%
  tidyr::separate(V4, into = c("transcript_id", "orf_info"), sep = ":") %>%
  dplyr::select(transcript_id, orf_info)

# Process BLASTp results
colnames(uniprot_blastp) <- c("transcript_id", "hit_id", "percent_identity", "alignment_length",
                              "mismatches", "gap_opens", "query_start", "query_end",
                              "subject_start", "subject_end", "evalue", "bit_score")

# Remove protein ID from "transcript_id"
uniprot_blastp <- uniprot_blastp %>%
  mutate(transcript_id = stringr::str_replace(transcript_id, "^((?:[^.]+\\.){2}[^.]+).*", "\\1"))

# Process BLASTx results
colnames(uniprot_blastx) <- c("transcript_id", "hit_id", "percent_identity", "alignment_length",
                              "mismatches", "gap_opens", "query_start", "query_end",
                              "subject_start", "subject_end", "evalue", "bit_score")

pfam_df <- as.data.frame(pfam)
pfam_df <- pfam_df %>%
  mutate(transcript_id = stringr::str_replace(query_name, "^((?:[^.]+\\.){2}[^.]+).*", "\\1"))

#Filter Pfam data
# Calculate alignment coverage (%)
pfam_df <- pfam_df %>%
  mutate(
    coverage = (ali_to - ali_from + 1) / domain_len * 100
  )
# Define thresholds
cevalue_threshold <- 1e-5    # Conditional E-value < 1e-5
#bit_score_threshold <- 20    # Bit score > 20: real threshold
bit_score_threshold <- 10 #just for testing purpose
coverage_threshold <- 50     # Coverage > 50%

# Apply filtering
pfam_filtered <- pfam_df %>%
  filter(
    domain_cevalue < cevalue_threshold,
    domain_score > bit_score_threshold,
    coverage > coverage_threshold
  )
# Select relevant columns
pfam_info <- pfam_filtered %>%
  dplyr::select(
    transcript_id,
    domain_accession,
    description
  )
#Add into gtf annotation: has_ORF (TRUE/FALSE), best_SwissProt_blastp_hit, best_SwissProt_blastx_hit, PFAM
# Annotate transcripts with ORF information
transcripts$has_ORF <- transcripts$transcript_id %in% orfs$transcript_id
transcripts$has_ORF <- ifelse(transcripts$type=="gene",NA,transcripts$has_ORF)
# Filter BLAST hits based on thresholds
# Define significance thresholds
#blast_evalue_threshold <- 1e-5 #correct threshold
blast_evalue_threshold <- 1 #just for testing purpose
blast_identity_threshold <- 25 #the percentage of amino acids in the alignment that are identical between the query (our transcript) and the subject (the database protein).
#uniprot_blastp[uniprot_blastp$percent_identity>49,]
#head(uniprot_blastp)
signif_blastp <- uniprot_blastp %>%
  dplyr::filter(evalue < blast_evalue_threshold & percent_identity > blast_identity_threshold) #%>%
# dplyr::group_by(transcript_id) %>%
# dplyr::summarize(num_blastp_hits = dplyr::n())

# Filter BLASTx hits based on thresholds
signif_blastx <- uniprot_blastx %>%
  dplyr::filter(evalue < blast_evalue_threshold & percent_identity > blast_identity_threshold) #%>%
# dplyr::group_by(transcript_id) %>%
# dplyr::summarize(num_blastx_hits = dplyr::n())

# Select relevant columns
signif_blastp_info <- signif_blastp %>%
  dplyr::select(
    transcript_id,
    hit_id
  )
# Select relevant columns
signif_blastx_info <- signif_blastx %>%
  dplyr::select(
    transcript_id,
    hit_id
  )

colnames(pfam_info)=c("transcript_id","Pfam_domain_accession", "Pfam_description" )
# Parse hit_id to extract Primary Accession and UniProt ID
signif_blastp_info <- signif_blastp_info %>%
  mutate(
    # Extract Primary Accession 
    UniProt_accession = str_extract(hit_id, "(?<=\\|)[A-Z0-9]+(?=\\|)"),
    # Extract UniProt ID 
    UniProt_entry_name = str_extract(hit_id, "(?<=\\|)[A-Z0-9_]+$")
  )
signif_blastx_info <- signif_blastx_info %>%
  mutate(
    UniProt_accession = str_extract(hit_id, "(?<=\\|)[A-Z0-9]+(?=\\|)"),
    UniProt_entry_name = str_extract(hit_id, "(?<=\\|)[A-Z0-9_]+$")
  )
signif_blastp <- signif_blastp %>%
  mutate(
    UniProt_accession = str_extract(hit_id, "(?<=\\|)[A-Z0-9]+(?=\\|)"),
    UniProt_entry_name = str_extract(hit_id, "(?<=\\|)[A-Z0-9_]+$")
  )
signif_blastx <- signif_blastx %>%
  mutate(
    UniProt_accession = str_extract(hit_id, "(?<=\\|)[A-Z0-9]+(?=\\|)"),
    UniProt_entry_name = str_extract(hit_id, "(?<=\\|)[A-Z0-9_]+$")
  )

best_blastp <- signif_blastp %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::arrange(desc(percent_identity), evalue) %>%
  slice(1) %>%  # Select the top hit based on the sorting
  dplyr::ungroup() %>%
  dplyr::select(transcript_id, best_SwissProt_blastp_hit = UniProt_entry_name)

best_blastx <- signif_blastx %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::arrange(desc(percent_identity), evalue) %>%
  slice(1) %>%  
  dplyr::ungroup() %>%
  dplyr::select(transcript_id, best_SwissProt_blastx_hit = UniProt_entry_name)

best_pfam <- pfam_filtered %>%
  dplyr::group_by(transcript_id) %>%
  dplyr::arrange(desc(coverage), domain_cevalue, desc(domain_score)) %>%
  dplyr::slice(1) %>%  
  dplyr::ungroup() %>%
  dplyr::select(transcript_id, best_Pfam_hit = description)
# Merge with transcripts
transcripts <- transcripts %>%
  left_join(best_blastp, by = "transcript_id") %>%
  left_join(best_blastx, by = "transcript_id") %>%
  left_join(best_pfam, by = "transcript_id")


# Extract and rename hits from pfam_filtered
pfam_hits <- pfam_filtered %>%
  dplyr::select(transcript_id, description) %>%
  dplyr::rename(hit = description)
# Extract and rename hits from signif_blastx
blastx_hits <- signif_blastx %>%
  dplyr::select(transcript_id, UniProt_entry_name) %>%
  dplyr::rename(hit = UniProt_entry_name)
# Extract and rename hits from signif_blastp
blastp_hits <- signif_blastp %>%
  dplyr::select(transcript_id, UniProt_entry_name) %>%
  dplyr::rename(hit = UniProt_entry_name)

# Combine all hits into one dataframe
all_hits_combined <- bind_rows(
  pfam_hits,
  blastx_hits,
  blastp_hits
) %>%
  # Remove any NA hits
  dplyr::filter(!is.na(hit)) %>%
  # Group by transcript_id
  dplyr::group_by(transcript_id) %>%
  # Concatenate unique hits separated by commas
  dplyr::summarize(all_hits = paste(unique(hit), collapse = ", "), .groups = 'drop')

# Merge the all_hits_combined with transcripts
transcripts_with_hits <- transcripts %>%
  left_join(all_hits_combined, by = "transcript_id")
# Add gene_id into all_hits_combined
gene_trans <- na.omit(unique(transcripts[,c("gene_id","transcript_id")]))
all_hits_combined_genes <- all_hits_combined %>%
  left_join(gene_trans, by = "transcript_id")
all_hits_combined_genes <- all_hits_combined_genes[,-1]
all_hits_combined_genes <- bind_rows(all_hits_combined_genes) %>%
  # Remove any NA hits
  dplyr::filter(!is.na(all_hits)) %>%
  # Group by transcript_id
  dplyr::group_by(gene_id) %>%
  # Concatenate unique hits separated by commas
  dplyr::summarize(all_hits = paste(unique(all_hits), collapse = ", "), .groups = 'drop')

# Remove redundant MSTRG (optional) - remove transcripts with no functional information

functional_transcripts <- transcripts_with_hits[transcripts_with_hits$type=="gene" | 
                                                  (transcripts_with_hits$type!="gene" & (transcripts_with_hits$has_ORF==TRUE | 
                                                                                           !(is.na(transcripts_with_hits$best_SwissProt_blastp_hit)) | 
                                                                                           !(is.na(transcripts_with_hits$best_SwissProt_blastx_hit)) |
                                                                                           !(is.na(transcripts_with_hits$best_Pfam_hit)))), ]
#Identify genes without any transcripts (after filtering in previous step)
genes_tb <- as.data.frame(table(functional_transcripts$gene_id))
functional_transcripts <- functional_transcripts[functional_transcripts$gene_id %in% genes_tb$Var1[genes_tb$Freq>1],]


#Add functional description for genes
functional_genes <- functional_transcripts[functional_transcripts$type=="transcript",c("gene_id","transcript_id","best_SwissProt_blastp_hit","best_SwissProt_blastx_hit","best_Pfam_hit")]
functional_genes <- functional_genes[,-2]
# Aggregate functional annotations at the gene level
functional_genes_unique <- functional_genes %>%
  group_by(gene_id) %>%
  summarize(
    best_SwissProt_blastp_hit = {
      hits <- unique(na.omit(best_SwissProt_blastp_hit))
      if(length(hits) > 0) paste(hits, collapse = ", ") else NA_character_
    },
    best_SwissProt_blastx_hit = {
      hits <- unique(na.omit(best_SwissProt_blastx_hit))
      if(length(hits) > 0) paste(hits, collapse = ", ") else NA_character_
    },
    best_Pfam_hit = {
      hits <- unique(na.omit(best_Pfam_hit))
      if(length(hits) > 0) paste(hits, collapse = ", ") else NA_character_
    },
    .groups = 'drop' # To ungroup after summarization
  )

# View the aggregated data
print(functional_genes_unique)

# Separate the gene subset 
functional_genes_subset <- functional_transcripts %>%
  filter(type == "gene") 

n <- ncol(functional_genes_subset)
#remove "best_SwissProt_blastp_hit","best_SwissProt_blastx_hit","best_Pfam_hit", "all_hits" columns --> we will add it later
functional_genes_subset <- functional_genes_subset[,-c(n,n-1,n-2,n-3)]
# Perform a left_join to add functional annotations to the gene subset
functional_genes_annotated <- functional_genes_subset %>%
  left_join(functional_genes_unique, by = "gene_id") 

# Separate the non-gene subset
functional_non_genes_subset <- functional_transcripts %>%
  filter(type != "gene")


# Recombine the annotated gene subset with the non-gene subset
functional_transcripts_updated <- bind_rows(functional_genes_annotated, functional_non_genes_subset)

#functional_transcripts_updated[functional_transcripts_updated$gene_id=="MSTRG.127",]
# Arrange the data frame to preserve original order 
functional_transcripts_updated <- functional_transcripts_updated %>%
  arrange(match(gene_id, functional_transcripts$gene_id))

# Merge the all_hits_combined with genes
# Merge gene-level all_hits into transcripts_with_hits
functional_transcripts_updated <- functional_transcripts_updated %>%
  # Perform a left join on gene_id
  dplyr::left_join(all_hits_combined_genes, by = "gene_id", suffix = c("", "_gene")) %>%
  
  # Conditionally replace all_hits where type == "gene"
  dplyr::mutate(
    all_hits = if_else(
      type == "gene",
      all_hits_gene,  # Use gene-level all_hits
      all_hits        # Retain existing transcript-level all_hits
    )
  ) %>%
  
  # Remove the auxiliary all_hits_gene column
  dplyr::select(-all_hits_gene)


################################################################################
#---------------Step 11: Find overlapping genes and transcripts-----------------
################################################################################
# Convert functional_transcripts to GRanges
functional_transcripts_gr <- makeGRangesFromDataFrame(
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
  group_by(gene_id) %>%
  summarize(
    overlapped_predicted_gene = {
      hits <- unique(na.omit(overlapped_predicted_gene))
      if(length(hits) > 0) paste(hits, collapse = "; ") else NA_character_
    },
    .groups = 'drop' # To ungroup after summarization
  )

# Merge overlaps with functional_transcripts_updated
functional_transcripts_updated <- functional_transcripts_updated %>%
  left_join(
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
  overlaps <- findOverlaps(transcript_grl, transcript_grl)
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
  functional_transcripts_updated <- makeGRangesFromDataFrame(
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
  
  ref_anno_gr <- makeGRangesFromDataFrame(
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
    group_by(gene_id) %>%
    summarize(
      overlapped_ref_gene = {
        hits <- unique(na.omit(overlapped_ref_gene))
        if(length(hits) > 0) paste(hits, collapse = "; ") else NA_character_
      },
      .groups = 'drop' # To ungroup after summarization
    )
  
  # Merge overlaps with functional_transcripts_updated
  functional_transcripts_updated <- functional_transcripts_updated %>%
    left_join(
      overlaps_with_ref_dt_uniq,
      by = "gene_id"
    )
  
} else {
  # If reference annotation is not available, partially skip Step 11
  message("Step 11: Reference annotation not available - Skipping overlapping genes.")
}
# Merge overlap information back to functional_transcripts_updated
functional_transcripts_updated <- as.data.frame(functional_transcripts_updated)
functional_transcripts_updated <- functional_transcripts_updated %>%
  left_join(transcript_overlap_df, by = "gene_id")
################################################################################
#-----------Step 12: Find reversed duplicates-----------------------------------
################################################################################
if (ref_available) {  # Only perform if reference annotation is available
  ref_anno_dt <- as.data.table(as.data.frame(ref_anno_gr))
  ref_anno_dt <- ref_anno_dt[, .(gene_id, strand_ref = as.character(strand))]
  head(ref_anno_dt)
  
  functional_transcripts_updated <- functional_transcripts_updated %>%
    left_join(ref_anno_dt, by = c("overlapped_ref_gene" = "gene_id"))
  
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
#-----------Step 13: Find fragmented and chimeric genes (with reference)--------
################################################################################
if (ref_available) {  # Only perform if reference annotation is available
  # One to many: one reference gene matches to many predictions
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
#-----------Step 13: Find fragmented and chimeric genes (no reference)----------
################################################################################
# ----------------------------------------
# Identify fragmented genes de novo
# ----------------------------------------
# Subset to gene entries
#head(functional_transcripts_updated)
genes_df <- functional_transcripts_updated %>%
  dplyr::filter(type == "gene")

# Convert to GRanges
genes_gr <- makeGRangesFromDataFrame(
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
max_distance <- 1000
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
  functional_entries <- str_trim(functional_entries)
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
  mutate(
    is_potentially_fragmented_de_novo = ifelse(fragment_cluster %in% clusters_to_flag, TRUE, FALSE)
  )

# Update the main dataframe with fragmentation information
functional_transcripts_updated <- functional_transcripts_updated %>%
  left_join(
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
large_intron_threshold <- 100000

# Identify genes with large introns
genes_with_large_introns <- intron_info %>%
  dplyr::filter(!is.na(width) & width > large_intron_threshold) %>%
  pull(gene_id) %>%
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
    x <- str_trim(x)
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

functional_annotations_transcripts <- functional_transcripts_updated %>%
  dplyr::filter(type == "transcript") %>%
  dplyr::group_by(gene_id) %>%
  dplyr::filter(n() > 1) %>%
  # Keep genes with >1 transcript
  dplyr::ungroup()   %>%
  dplyr::group_by( gene_id) %>%
  dplyr::summarize(
    functional_entries = list(unlist(all_hits)),
    transcripts_nr = n(),
    .groups = 'drop'
  )  %>%
  dplyr::mutate(
    common_entries = map(functional_entries, compute_common_entries_all),
    # Determine if any transcript has functional annotations
    has_any_annotations = map_lgl(functional_entries, ~ any(sapply(.x, function(fe) {
      fe_clean <- fe[!is.na(fe) & fe != ""]
      length(fe_clean) > 0
    }))),
    # Conditionally flag as chimeric based on common functional entries
    is_potentially_chimeric_function_de_novo = if_else(
      has_any_annotations,
      # No common entries across all transcripts
      (map_lgl(common_entries, ~ is.null(.x) || length(.x) == 0)),
      FALSE
    )
  ) %>%
  # Merge overlap information
  left_join(transcript_overlap_df, by = "gene_id") %>%
  # Update the flag to exclude genes with overlapping transcripts
  mutate(
    is_potentially_chimeric_function_de_novo = if_else(
      is_potentially_chimeric_function_de_novo & has_all_transcripts_overlapping == FALSE,
      TRUE,
      FALSE
    )
  )


#functional_annotations_transcripts[1,2][[1]]
# Update the main dataframe with fragmentation information
functional_transcripts_updated <- functional_transcripts_updated %>%
  left_join(
    functional_annotations_transcripts %>% 
      dplyr::select(gene_id, is_potentially_chimeric_function_de_novo),
    by = "gene_id"
  )


# Continue with your existing R code, ensuring all file paths are constructed using `functional_dir` and `annotation_dir`.

# ---------------------------
# Save the Updated GTF
# ---------------------------
final_gtf_path <- file.path(functional_dir, "final_annotation.gtf")
# Convert to GRanges for export 
if (!inherits(functional_transcripts_updated, "GRanges")) {
  functional_transcripts_updated <- makeGRangesFromDataFrame(
    functional_transcripts_updated,
    keep.extra.columns = TRUE,
    seqnames.field = "seqnames",
    start.field = "start",
    end.field = "end",
    strand.field = "strand"
  )
}
export(functional_transcripts_updated, final_gtf_path, format = "gtf")
#export(functional_transcripts_updated, file = final_gtf_path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

cat("Functional annotation completed successfully. Final GTF saved at:", final_gtf_path, "\n")
