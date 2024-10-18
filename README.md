# Genome-annotation-pipeline

Make the script executable:
```bash
chmod +x rna_seq_pipeline.sh
```

Required inputs:
- **rRNA Reference FASTA:** Path to rrna.fa.
- **Genome reference FASTA:** Path to genome.fa.
- **Genome rnnotation GTF:** Path to reference.gtf.
- **BLAST databases:** Ensure BLAST NR and SwissProt databases are downloaded and formatted.
- **PFAM database:** Download and prepare the PFAM database for pfam_scan.pl.
- **rRNA STAR Index:** Ensure you have a STAR index for rRNA references.
- **Genome STAR Index:** Ensure you have a STAR index for the genome.

To Run All Steps:
```bash
./rna_seq_pipeline.sh --genomeDir /path/to/genomeDir \
                      --rrnaRef /path/to/rrna.fa \
                      --genomeRef /path/to/genome.fa \
                      --genomeGTF /path/to/reference.gtf \
                      --blastDB_NR /path/to/blastdb_nr \
                      --blastDB_SwissProt /path/to/blastdb_swissprot \
                      --pfamDB /path/to/pfam_database \
                      --dataDir /path/to/data \
                      --outputDir /path/to/output \
                      --threads 8 \
                      --all
```

To Run Specific Steps (e.g., Steps 2, 4, and 6):
```bash
./rna_seq_pipeline.sh --genomeDir /path/to/genomeDir \
                      --rrnaRef /path/to/rrna.fa \
                      --genomeRef /path/to/genome.fa \
                      --genomeGTF /path/to/reference.gtf \
                      --blastDB_NR /path/to/blastdb_nr \
                      --blastDB_SwissProt /path/to/blastdb_swissprot \
                      --pfamDB /path/to/pfam_database \
                      --dataDir /path/to/data \
                      --outputDir /path/to/output \
                      --threads 8 \
                      --steps 2,4,6
```

To Display Help:
```bash
./rna_seq_pipeline.sh --help
```
