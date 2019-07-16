# USEARCH 11 

## List of all commands

    - _E.g._ `s`

#### Sequence alignment commands

 - See also: [Output formats](https://www.drive5.com/usearch/manual/output_files.html),  [Accept hit options](https://www.drive5.com/usearch/manual/accept_options.html), [Alignment parameters](https://www.drive5.com/usearch/manual/aln_params.html)

 - **allpairs_global**  Align all pairs in FASTx file using global alignment
    -  _E.g._ `usearch -allpairs_global pairs.fasta -alnout results.aln`
  - **allpairs_local**  Align all pairs in FASTx file using local alignment
    - _E.g._ `usearch -allpairs_local seqs.fasta -alnout results.aln`
  - **pairs_global**  Align pairs of sequences in FASTx file using global alignment
    - _E.g._ `usearch -pairs_global pairs.fasta -alnout results.aln`
  - **pairs_local**  Align pairs of sequences in FASTx file using local alignment
    - _E.g._ `usearch -pairs_local pairs.fasta -alnout results.aln`
 

#### Chimera detection and filtering
- **annot**  Annotate OTU sequences as known (mock or large ref. db.), chimeric etc.
    - _E.g._ `usearch -threads 8 -annot mock_reads.fq -knowndb mock_ref.fa -db silva.udb -tabbedout annot.txt -fastqout annot.fq`
- **uchime2_ref**  Chimera search using UCHIME2 algorithm
  - Better using _unoise3_ or _cluster\_otus_
  - _E.g._ `usearch -uchime_ref reads.fasta -db 16s_ref.udb -uchimeout out.txt -strand plus -mode sensitive`
  - [Read more](https://www.drive5.com/usearch/manual/cmd_uchime2_ref.html)
- **uchime3_denovo**  Chimera search using UCHIME3 de-novo algorithm
  - Better using _unoise3_ or _cluster\_otus_
  - [Read more](https://www.drive5.com/usearch/manual/cmd_uchime3_denovo.html)
- **unoise3**  Denoise amplicon reads
  - _E.g._ `usearch -unoise3 uniques.fa -zotus zotus.fa -tabbedout unoise3.txt`
  - [Read more](https://www.drive5.com/usearch/manual/cmd_unoise3.html)

#### Sequence, tree and graph-based clustering
- **closed_ref**  Make OTU table using closed-reference clustering
   - Not recommended
- **cluster_aggd**  Cluster distance matrix using agglomerative clustering
    - _E.g._ `usearch -cluster_aggd mx.txt -treeout clusters.tree -clusterout clusters.txt   -id 0.80 -linkage min`
    - See _calc\_distmx_
- **cluster_edges**  Find connected components of graph (single-linkage clustering)
- **cluster_fast**  Cluster sequences using UCLUST
- **cluster_otus**  Cluster sequences using UPARSE
- **cluster_smallmem**  Cluster sequencees using UCLUST
- **cluster_tree**  Construct clusters from tree using distance cutoff

#### Distance matrices
- **calc_distmx**  Calculate sparse distance matrix
- **calc_lcr_probs**  Calculate Lowest Common Rank probabilities from dist. matrix with taxonomy
- **distmx_split_identity**  Split distance matrix into test/training pair for CVI
- **tree2distmx**  Calculate distance matrix implied by tree

#### Commands for diversity analysis
- **alpha_di**v  Calculate alpha diversity metric(s) from OTU table
- **alpha_div_rare**  Calculate alpha diversity metric(s) from OTU table with rarefaction
- **alpha_div_sig**  Statistical significance of alpha diversity correlation with metadata
- **beta_div**  Calculate beta diversity metric(s) from OTU table

#### Commands for reads in FASTQ format
- **fastq_chars**  Report frequencies of Q score ASCII characters in FASTQ file
- **fastq_eestats**  Report quality/e.e. per position for reads in FASTQ file
- **fastq_eestats2**  Report number of reads retained at difference length and e.e. cutoffs
- **fastq_filter**  Filter reads in FASTQ file by e.e. and other criteria
- **fastq_join**  Concatenate forward (R1) and reverse (R2) paired reads
- **fastq_mergepairs**  Assemble (merge) paired reads
- **fastq_sra_splitpairs**  Recover paired reads from SRA interleaved or concatenated format

#### Commands for sequences in FASTx format (FASTA and FASTQ)
- **allpairs_global**  Align all pairs in FASTx file using global alignment
- **allpairs_local**  Align all pairs in FASTx file using local alignment
- **fasta_explode**  De-unique FASTA file with size=nnn annotations
- **fasta_stripgaps**  Remove gap symbols from FASTA file
- **fastx2qiime**  Convert sample labels from usearch to QIIME format
- **fastx_demux**  Assign reads to samples (demultiplex)
- **fastx_findorfs**  Identify ORFs in nucleotide sequences
- **fastx_get_sample_names**  Extract sample names from FASTx file
- **fastx_getlabels**  Extract sequence labels from FASTx file
- **fastx_getseq**  Extract sequence(s) matching label
- **fastx_getseqs**  Extract sequence(s) matching labels
- **fastx_getsubseq**  Extract subsequence given label, start, stop
- **fastx_info**  Report summary information about a FASTx file
- **fastx_learn**  Estimate error rates from amplicon reads
- **fastx_mask**  Mask low-complexity sequence
- **fastx_relabel**  Re-label sequences in FASTx file with prefix plus sequential number
- **fastx_revcomp**  Reverse-complement nucleotide sequence
- **fastx_split**  Divide sequences in FASTx file into given number of files
- **fastx_strip_annots**  Remove usearch-style annotations (name=xxx) from FASTx file
- **fastx_subsample**  Extract random sub-sample from FASTx file
- **fastx_syncpairs**  Sort forward and reverse reads into the same order
- **fastx_trim_primer**  Remove primer-binding sequence from FASTx file
- **fastx_truncate**  Truncate sequences in FASTx file
- **fastx_uniques**  Identify unique sequences in FASTx file (dereplicate)
- **fastx_uniques_persample**  Identify unique sequences per sample in FASTx file
- **filter_lowc**  Filter low-complexity sequences from FASTx file
- **filter_phix**  Remove PhiX spike sequences from FASTx file
- **sortbylength**  Sort sequences in FASTx file by decreasing length
- **sortbysize**  Sort sequences in FASTx file by decreasing size=nnn

#### Machine learning and finding informative OTUs
 - **forest_classify**  Classify data using random forest
 - **forest_train**  Train random forest
 - **otutab_core**  Identify core microbiome in OTU table
 - **otutab_forest_classify**  Classify samples using random forest
 - **otutab_forest_train**  Train random forest on OTU table
 - **otutab_select**  Identify OTUs which are informative (predictive of metadata)

#### Miscellaneous commands
 - **search_16s**  Identify 16S sequences in chromosomes or contigs
   - `usearch -makeudb_usearch gg97.fa -wordlength 13 -output gg97.udb`
   - `usearch -udb2bitvec gg97.udb -output gg97.bitvec`
   - ``
 - **udb2bitvec**  Create database for search_16s command

#### Commands for OTU analysis and denoising
 - **alpha_div**  Calculate alpha diversity metric(s) from OTU table
 - **alpha_div_rare**  Calculate alpha diversity metric(s) from OTU table with rarefaction
 - **alpha_div_sig**  Statistical significance of alpha diversity correlation with metadata
 - **annot**  Annotate OTU sequences as known (mock or large ref. db.), chimeric etc.
 - **beta_div**  Calculate beta diversity metric(s) from OTU table
 - **closed_ref**  Make OTU table using closed-reference clustering
 - **cluster_aggd**  Cluster distance matrix using agglomerative clustering
 - **cluster_otus**  Cluster sequences using UPARSE
 - **fastx_learn**  Estimate error rates from amplicon reads
 - **filter_lowc**  Filter low-complexity sequences from FASTx file
 - **filter_phix**  Remove PhiX spike sequences from FASTx file
 - **nbc_tax**  Predict taxonomy using RDP Naive Bayesian Classifier algorithm
 - **otutab**  Generate OTU table
 - **otutab2biom**  Convert OTU table from tabbed to biom (json) format
 - **otutab_binary**  Convert OTU table with counts to presence(1)/absence(0)
 - **otutab_core**  Identify core microbiome in OTU table
 - **otutab_counts2freqs**  Convert counts to frequencies in OTU table
 - **otutab_forest_classify**  Classify samples using random forest
 - **otutab_forest_train**  Train random forest on OTU table
 - **otutab_group**  Sum subsets of samples in OTU table
 - **otutab_merge**  Merge two or more OTU tables
 - **otutab_octave**  Generate octave plot visualizing OTU abundance distribution
 - **otutab_otu_subset**  Extract subset of OTUs from OTU table
 - **otutab_otus**  Extract OTU names from OTU table
 - **otutab_rare**  Rarefy OTU table so that samples have same number of reads
 - **otutab_sample_subset**  Extract subset of samples from OTU table
 - **otutab_samples**  Extract sample names from OTU table
 - **otutab_select**  Identify OTUs which are informative (predictive of metadata)
 - **otutab_sortotus**  Sort OTU table in order of decreasing OTU size
 - **otutab_stats**  Report summary information about OTU table
 - **otutab_trim**  Trim OTU table to remove small counts, OTU and/or samples
 - **otutab_xtalk**  Estimate and filter cross-talk in OTU table
 - **qiimemap2otutab**  Convert QIIME map file to OTU table
 - **search_oligodb**  Search for matches to short nucleotide sequences, e.g. primers
 - **search_pcr**  In-silico PCR, search for matches to pairs of primers in database
 - **search_pcr2**  In-silico PCR, search for matches to primer pair
 - **search_phix**  Search for matches to PhiX sequence
 - **sinaps**  Predict traits
 - **sintax**  Predict taxonomy using SINTAX algorithm
 - **sintax_summary**  Generate summary report from sintax output
 - **tabbed2otutab**  Convert read mapping file (read+OTU) to OTU table
 - **uchime2_ref**  Chimera search using UCHIME2 algorithm
 - **uchime3_denovo**  Chimera search using UCHIME3 de-novo algorithm
 - **unbias**  Correct abundance bias in OTU table
 - **unoise3**  Denoise amplicon reads
 - **uparse_ref**  Classify sequences derived from mock community sample

#### OTU table commands
 - **alpha_div**  Calculate alpha diversity metric(s) from OTU table
 - **alpha_div_rare**  Calculate alpha diversity metric(s) from OTU table with rarefaction
 - **alpha_div_sig**  Statistical significance of alpha diversity correlation with metadata
 - **beta_div**  Calculate beta diversity metric(s) from OTU table
 - **otutab**  Generate OTU table
 - **otutab2biom**  Convert OTU table from tabbed to biom (json) format
 - **otutab_binary**  Convert OTU table with counts to presence(1)/absence(0)
 - **otutab_core**  Identify core microbiome in OTU table
 - **otutab_counts2freqs**  Convert counts to frequencies in OTU table
 - **otutab_forest_classify**  Classify samples using random forest
 - **otutab_forest_train**  Train random forest on OTU table
 - **otutab_group**  Sum subsets of samples in OTU table
 - **otutab_merge**  Merge two or more OTU tables
 - **otutab_octave**  Generate octave plot visualizing OTU abundance distribution
 - **otutab_otu_subset**  Extract subset of OTUs from OTU table
 - **otutab_otus**  Extract OTU names from OTU table
 - **otutab_rare**  Rarefy OTU table so that samples have same number of reads
 - **otutab_sample_subset**  Extract subset of samples from OTU table
 - **otutab_samples**  Extract sample names from OTU table
 - **otutab_select**  Identify OTUs which are informative (predictive of metadata)
 - **otutab_sortotus**  Sort OTU table in order of decreasing OTU size
 - **otutab_stats**  Report summary information about OTU table
 - **otutab_trim**  Trim OTU table to remove small counts, OTU and/or samples
 - **otutab_xtalk**  Estimate and filter cross-talk in OTU table
 - **qiimemap2otutab**  Convert QIIME map file to OTU table
 - **tabbed2otutab**  Convert read mapping file (read+OTU) to OTU table
 - **unbias**  Correct abundance bias in OTU table

#### Next-generation reads
 - **fastq_chars**  Report frequencies of Q score ASCII characters in FASTQ file
 - **fastq_eestats**  Report quality/e.e. per position for reads in FASTQ file
 - **fastq_eestats2**  Report number of reads retained at difference length and e.e. cutoffs
 - **fastq_filter**  Filter reads in FASTQ file by e.e. and other criteria
 - **fastq_join**  Concatenate forward (R1) and reverse (R2) paired reads
 - **fastq_mergepairs**  Assemble (merge) paired reads
 - **fastq_sra_splitpairs**  Recover paired reads from SRA interleaved or concatenated format
 - **fastx2qiime**  Convert sample labels from usearch to QIIME format
 - **fastx_demux**  Assign reads to samples (demultiplex)
 - **fastx_findorfs**  Identify ORFs in nucleotide sequences
 - **fastx_get_sample_names**  Extract sample names from FASTx file
 - **fastx_info**  Report summary information about a FASTx file
 - **fastx_learn**  Estimate error rates from amplicon reads
 - **fastx_subsample**  Extract random sub-sample from FASTx file
 - **fastx_syncpairs**  Sort forward and reverse reads into the same order
 - **fastx_trim_primer**  Remove primer-binding sequence from FASTx file
 - **fastx_truncate**  Truncate sequences in FASTx file
 - **filter_lowc**  Filter low-complexity sequences from FASTx file
 - **filter_phix**  Remove PhiX spike sequences from FASTx file
 - **search_oligodb**  Search for matches to short nucleotide sequences, e.g. primers
 - **search_pcr**  In-silico PCR, search for matches to pairs of primers in database
 - **search_pcr2**  In-silico PCR, search for matches to primer pair
 - **search_phix**  Search for matches to PhiX sequence

#### Sequence database search
 - **makeudb_sintax**  Make UDB database file for sintax
 - **makeudb_ublast**  Make UDB database file for ublast
 - **makeudb_usearch**  Make UDB database file for usearch_global
 - **search_exact**  Search for identical sequences
 - **search_global**  Search database using global alignment without speed heuristics
 - **search_local**  Search database using locbal alignment without speed heuristics
 - **search_oligodb**  Search for matches to short nucleotide sequences, e.g. primers
 - **search_pcr**  In-silico PCR, search for matches to pairs of primers in database
 - **search_pcr2**  In-silico PCR, search for matches to primer pair
 - **search_peptidedb**  Search for matches to short peptide sequences
 - **search_phix**  Search for matches to PhiX sequence
 - **ublast**  Fast database search using local alignment (much faster than BLAST)
 - **usearch_global**  Fast database search using global alignment
 - **usearch_local**  Fast database search using local alignment

#### Taxonomy commands
 - **calc_lcr_probs**  Calculate Lowest Common Rank probabilities from dist. matrix with taxonomy
 - **makeudb_sintax**  Make UDB database file for sintax
 - **nbc_tax**  Predict taxonomy using RDP Naive Bayesian Classifier algorithm
 - **sintax**  Predict taxonomy using SINTAX algorithm
 - **sintax_summary**  Generate summary report from sintax output

#### Tree commands
 - **calc_distmx**  Calculate sparse distance matrix
 - **cluster_tree**  Construct clusters from tree using distance cutoff
 - **subtree**  Extract subtree under given node
 - **tree2distmx**  Calculate distance matrix implied by tree
 - **tree_cvt**  Convert tree between tabbed and Newick formats
 - **tree_subset**  Extract tree for subset of leaves

#### Labels and annotations
 - **fastx2qiime**  Convert sample labels from usearch to QIIME format
 - **fastx_getlabels**  Extract sequence labels from FASTx file
 - **fastx_relabel**  Re-label sequences in FASTx file with prefix plus sequential number
 - **fastx_strip_annots**  Remove usearch-style annotations (name=xxx) from FASTx file
 - **otutab_otus**  Extract OTU names from OTU table
 - **otutab_samples**  Extract sample names from OTU table
  
