# USEARCH 11  
 
:no_entry: = Obsolete/Not recommended 
 
:trophy: = Selected topic 
 
 
## Sections 
 
 - [USEARCH topics](topics.md) 
  
  
## List of all commands 
 
#### Sequence alignment commands 
 
 - See also: [Output formats](https://www.drive5.com/usearch/manual/output_files.html),  [Accept hit options](https://www.drive5.com/usearch/manual/accept_options.html), [Alignment parameters](https://www.drive5.com/usearch/manual/aln_params.html) 
 
 - **allpairs_global**  Align all pairs in FASTx file using global alignment [[link](https://www.drive5.com/usearch/manual/cmd_allpairs_global.html)]
    -  _E.g._ `usearch -allpairs_global pairs.fasta -alnout results.aln` 
  - **allpairs_local**  Align all pairs in FASTx file using local alignment [[link](https://www.drive5.com/usearch/manual/cmd_allpairs_local.html)]
    - _E.g._ `usearch -allpairs_local seqs.fasta -alnout results.aln` 
  - **pairs_global**  Align pairs of sequences in FASTx file using global alignment [[link](https://www.drive5.com/usearch/manual/cmd_pairs_global.html)]
    - _E.g._ `usearch -pairs_global pairs.fasta -alnout results.aln` 
  - **pairs_local**  Align pairs of sequences in FASTx file using local alignment [[link](https://www.drive5.com/usearch/manual/cmd_pairs_local.html)]
    - _E.g._ `usearch -pairs_local pairs.fasta -alnout results.aln` 
  
 
#### Chimera detection and filtering 
- **annot**  Annotate OTU sequences as known (mock or large ref. db.), chimeric etc. [[link](https://www.drive5.com/usearch/manual/cmd_annot.html)]
    - _E.g._ `usearch -threads 8 -annot mock_reads.fq -knowndb mock_ref.fa -db silva.udb -tabbedout annot.txt -fastqout annot.fq` 
- **uchime2_ref** :no_entry:  Chimera search using UCHIME2 algorithm [[link](https://www.drive5.com/usearch/manual/cmd_uchime2_ref.html)]
  - Better using _unoise3_ or _cluster\_otus_ 
  - _E.g._ `usearch -uchime_ref reads.fasta -db 16s_ref.udb -uchimeout out.txt -strand plus -mode sensitive` 
  - [Read more](https://www.drive5.com/usearch/manual/cmd_uchime2_ref.html) 
- **uchime3_denovo**  :no_entry:  Chimera search using UCHIME3 de-novo algorithm [[link](https://www.drive5.com/usearch/manual/cmd_uchime3_denovo.html)]
  - Better using _unoise3_ or _cluster\_otus_ 
  - [Read more](https://www.drive5.com/usearch/manual/cmd_uchime3_denovo.html) 
- **unoise3** :trophy:   Denoise amplicon reads [[link](https://www.drive5.com/usearch/manual/cmd_unoise3.html)]
  - _E.g._ `usearch -unoise3 uniques.fa -zotus zotus.fa -tabbedout unoise3.txt` 
  - [Read more](https://www.drive5.com/usearch/manual/cmd_unoise3.html) 
 
#### Sequence, tree and graph-based clustering 
- **closed_ref** :no_entry:   Make OTU table using closed-reference clustering [[link](https://www.drive5.com/usearch/manual/cmd_closed_ref.html)]
   - Not recommended 
- **cluster_aggd**  Cluster distance matrix using agglomerative clustering [[link](https://www.drive5.com/usearch/manual/cmd_cluster_aggd.html)]
    - _E.g._ `usearch -cluster_aggd mx.txt -treeout clusters.tree -clusterout clusters.txt   -id 0.80 -linkage min` 
    - See _calc\_distmx_ 
- **cluster_edges**  Find connected components of graph (single-linkage clustering) [[link](https://www.drive5.com/usearch/manual/cmd_cluster_edges.html)]
- **cluster_fast** :trophy:  Cluster sequences using UCLUST [[link](https://www.drive5.com/usearch/manual/cmd_cluster_fast.html)]
  - `usearch -cluster_fast query.fasta -id 0.9 -centroids nr.fasta -uc clusters.uc` 
  - [Read more](cluster_fast.md) 
- **cluster_otus** :trophy:  Cluster sequences using UPARSE [[link](https://www.drive5.com/usearch/manual/cmd_cluster_otus.html)]
- **cluster_smallmem**  Cluster sequencees using UCLUST [[link](https://www.drive5.com/usearch/manual/cmd_cluster_smallmem.html)]
- **cluster_tree**  Construct clusters from tree using distance cutoff [[link](https://www.drive5.com/usearch/manual/cmd_cluster_tree.html)]
 
#### Distance matrices 
- **calc_distmx**  Calculate sparse distance matrix [[link](https://www.drive5.com/usearch/manual/cmd_calc_distmx.html)]
- **calc_lcr_probs**  Calculate Lowest Common Rank probabilities from dist. matrix with taxonomy [[link](https://www.drive5.com/usearch/manual/cmd_calc_lcr_probs.html)]
- **distmx_split_identity**  Split distance matrix into test/training pair for CVI [[link](https://www.drive5.com/usearch/manual/cmd_distmx_split_identity.html)]
- **tree2distmx**  Calculate distance matrix implied by tree [[link](https://www.drive5.com/usearch/manual/cmd_tree2distmx.html)]
 
#### Commands for diversity analysis 
- **alpha_div**  Calculate alpha diversity metric(s) from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_alpha_div.html)]
- **alpha_div_rare**  Calculate alpha diversity metric(s) from OTU table with rarefaction [[link](https://www.drive5.com/usearch/manual/cmd_alpha_div_rare.html)]
- **alpha_div_sig**  Statistical significance of alpha diversity correlation with metadata [[link](https://www.drive5.com/usearch/manual/cmd_alpha_div_sig.html)]
- **beta_div**  Calculate beta diversity metric(s) from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_beta_div.html)]
 
#### Commands for reads in FASTQ format 
- **fastq_chars**  Report frequencies of Q score ASCII characters in FASTQ file [[link](https://www.drive5.com/usearch/manual/cmd_fastq_chars.html)]
- **fastq_eestats**  Report quality/e.e. per position for reads in FASTQ file [[link](https://www.drive5.com/usearch/manual/cmd_fastq_eestats.html)]
- **fastq_eestats2**  Report number of reads retained at difference length and e.e. cutoffs [[link](https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html)]
- **fastq_filter**  Filter reads in FASTQ file by e.e. and other criteria [[link](https://www.drive5.com/usearch/manual/cmd_fastq_filter.html)]
- **fastq_join**  Concatenate forward (R1) and reverse (R2) paired reads [[link](https://www.drive5.com/usearch/manual/cmd_fastq_join.html)]
- **fastq_mergepairs**  Assemble (merge) paired reads [[options](https://www.drive5.com/usearch/manual/merge_options.html)]  [[description](https://www.drive5.com/usearch/manual/cmd_fastq_mergepairs.html)]
- **fastq_sra_splitpairs**  Recover paired reads from SRA interleaved or concatenated format [[link](https://www.drive5.com/usearch/manual/cmd_fastq_sra_splitpairs.html)]
 
#### Commands for sequences in FASTx format (FASTA and FASTQ) 
- **allpairs_global**  Align all pairs in FASTx file using global alignment [[link](https://www.drive5.com/usearch/manual/cmd_allpairs_global.html)]
- **allpairs_local**  Align all pairs in FASTx file using local alignment [[link](https://www.drive5.com/usearch/manual/cmd_allpairs_local.html)]
- **fasta_explode**  De-unique FASTA file with size=nnn annotations [[link](https://www.drive5.com/usearch/manual/cmd_fasta_explode.html)]
- **fasta_stripgaps**  Remove gap symbols from FASTA file [[link](https://www.drive5.com/usearch/manual/cmd_fasta_stripgaps.html)]
- **fastx2qiime**  Convert sample labels from usearch to QIIME format [[link](https://www.drive5.com/usearch/manual/cmd_fastx2qiime.html)]
- **fastx_demux**  Assign reads to samples (demultiplex) [[link](https://www.drive5.com/usearch/manual/cmd_fastx_demux.html)]
- **fastx_findorfs**  Identify ORFs in nucleotide sequences [[link](https://www.drive5.com/usearch/manual/cmd_fastx_findorfs.html)]
- **fastx_get_sample_names**  Extract sample names from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_get_sample_names.html)]
- **fastx_getlabels**  Extract sequence labels from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_getlabels.html)]
- **fastx_getseq**  Extract sequence(s) matching label [[link](https://www.drive5.com/usearch/manual/cmd_fastx_getseq.html)]
- **fastx_getseqs**  Extract sequence(s) matching labels [[link](https://www.drive5.com/usearch/manual/cmd_fastx_getseqs.html)]
- **fastx_getsubseq**  Extract subsequence given label, start, stop [[link](https://www.drive5.com/usearch/manual/cmd_fastx_getsubseq.html)]
- **fastx_info**  Report summary information about a FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_info.html)]
- **fastx_learn**  Estimate error rates from amplicon reads [[link](https://www.drive5.com/usearch/manual/cmd_fastx_learn.html)]
- **fastx_mask**  Mask low-complexity sequence [[link](https://www.drive5.com/usearch/manual/cmd_fastx_mask.html)]
- **fastx_relabel**  Re-label sequences in FASTx file with prefix plus sequential number [[link](https://www.drive5.com/usearch/manual/cmd_fastx_relabel.html)]
- **fastx_revcomp**  Reverse-complement nucleotide sequence [[link](https://www.drive5.com/usearch/manual/cmd_fastx_revcomp.html)]
- **fastx_split**  Divide sequences in FASTx file into given number of files [[link](https://www.drive5.com/usearch/manual/cmd_fastx_split.html)]
- **fastx_strip_annots**  Remove usearch-style annotations (name=xxx) from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_strip_annots.html)]
- **fastx_subsample**  Extract random sub-sample from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_subsample.html)]
- **fastx_syncpairs**  Sort forward and reverse reads into the same order [[link](https://www.drive5.com/usearch/manual/cmd_fastx_syncpairs.html)]
- **fastx_trim_primer**  Remove primer-binding sequence from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_trim_primer.html)]
- **fastx_truncate**  Truncate sequences in FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_truncate.html)]
- **fastx_uniques**  Identify unique sequences in FASTx file (dereplicate) [[link](https://www.drive5.com/usearch/manual/cmd_fastx_uniques.html)]
- **fastx_uniques_persample**  Identify unique sequences per sample in FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_uniques_persample.html)]
- **filter_lowc**  Filter low-complexity sequences from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_filter_lowc.html)]
- **filter_phix**  Remove PhiX spike sequences from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_filter_phix.html)]
- **sortbylength**  Sort sequences in FASTx file by decreasing length [[link](https://www.drive5.com/usearch/manual/cmd_sortbylength.html)]
- **sortbysize**  Sort sequences in FASTx file by decreasing size=nnn [[link](https://www.drive5.com/usearch/manual/cmd_sortbysize.html)]
 
#### Machine learning and finding informative OTUs 
 - **forest_classify**  Classify data using random forest [[link](https://www.drive5.com/usearch/manual/cmd_forest_classify.html)]
 - **forest_train**  Train random forest [[link](https://www.drive5.com/usearch/manual/cmd_forest_train.html)]
 - **otutab_core**  Identify core microbiome in OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_core.html)]
 - **otutab_forest_classify**  Classify samples using random forest [[link](https://www.drive5.com/usearch/manual/cmd_otutab_forest_classify.html)]
 - **otutab_forest_train**  Train random forest on OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_forest_train.html)]
 - **otutab_select**  Identify OTUs which are informative (predictive of metadata) [[link](https://www.drive5.com/usearch/manual/cmd_otutab_select.html)]
 
#### Miscellaneous commands 
 - **search_16s**  Identify 16S sequences in chromosomes or contigs [[link](https://www.drive5.com/usearch/manual/cmd_search_16s.html)]
   - `usearch -search_16s contigs.fa -bitvecgg97.bitvec -fastaout 16s.fa` 
   - See [search_16s](search_16s.md) to prepare the database 
 - **udb2bitvec**  Create database for search_16s command [[link](https://www.drive5.com/usearch/manual/cmd_udb2bitvec.html)]
 
#### Commands for OTU analysis and denoising 
 - **alpha_div**  Calculate alpha diversity metric(s) from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_alpha_div.html)]
 - **alpha_div_rare**  Calculate alpha diversity metric(s) from OTU table with rarefaction [[link](https://www.drive5.com/usearch/manual/cmd_alpha_div_rare.html)]
 - **alpha_div_sig**  Statistical significance of alpha diversity correlation with metadata [[link](https://www.drive5.com/usearch/manual/cmd_alpha_div_sig.html)]
 - **annot**  Annotate OTU sequences as known (mock or large ref. db.), chimeric etc. [[link](https://www.drive5.com/usearch/manual/cmd_annot.html)]
 - **beta_div**  Calculate beta diversity metric(s) from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_beta_div.html)]
 - **closed_ref**  Make OTU table using closed-reference clustering [[link](https://www.drive5.com/usearch/manual/cmd_closed_ref.html)]
 - **cluster_aggd**  Cluster distance matrix using agglomerative clustering [[link](https://www.drive5.com/usearch/manual/cmd_cluster_aggd.html)]
 - **cluster_otus**  Cluster sequences using UPARSE [[link](https://www.drive5.com/usearch/manual/cmd_cluster_otus.html)]
 - **fastx_learn**  Estimate error rates from amplicon reads [[link](https://www.drive5.com/usearch/manual/cmd_fastx_learn.html)]
 - **filter_lowc**  Filter low-complexity sequences from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_filter_lowc.html)]
 - **filter_phix**  Remove PhiX spike sequences from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_filter_phix.html)]
 - **nbc_tax**  Predict taxonomy using RDP Naive Bayesian Classifier algorithm [[link](https://www.drive5.com/usearch/manual/cmd_nbc_tax.html)]
 - **otutab**  Generate OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab.html)]
 - **otutab2biom**  Convert OTU table from tabbed to biom (json) format [[link](https://www.drive5.com/usearch/manual/cmd_otutab2biom.html)]
 - **otutab_binary**  Convert OTU table with counts to presence(1)/absence(0) [[link](https://www.drive5.com/usearch/manual/cmd_otutab_binary.html)]
 - **otutab_core**  Identify core microbiome in OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_core.html)]
 - **otutab_counts2freqs**  Convert counts to frequencies in OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_counts2freqs.html)]
 - **otutab_forest_classify**  Classify samples using random forest [[link](https://www.drive5.com/usearch/manual/cmd_otutab_forest_classify.html)]
 - **otutab_forest_train**  Train random forest on OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_forest_train.html)]
 - **otutab_group**  Sum subsets of samples in OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_group.html)]
 - **otutab_merge**  Merge two or more OTU tables [[link](https://www.drive5.com/usearch/manual/cmd_otutab_merge.html)]
 - **otutab_octave**  Generate octave plot visualizing OTU abundance distribution [[link](https://www.drive5.com/usearch/manual/cmd_otutab_octave.html)]
 - **otutab_otu_subset**  Extract subset of OTUs from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_otu_subset.html)]
 - **otutab_otus**  Extract OTU names from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_otus.html)]
 - **otutab_rare**  Rarefy OTU table so that samples have same number of reads [[link](https://www.drive5.com/usearch/manual/cmd_otutab_rare.html)]
 - **otutab_sample_subset**  Extract subset of samples from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_sample_subset.html)]
 - **otutab_samples**  Extract sample names from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_samples.html)]
 - **otutab_select**  Identify OTUs which are informative (predictive of metadata) [[link](https://www.drive5.com/usearch/manual/cmd_otutab_select.html)]
 - **otutab_sortotus**  Sort OTU table in order of decreasing OTU size [[link](https://www.drive5.com/usearch/manual/cmd_otutab_sortotus.html)]
 - **otutab_stats**  Report summary information about OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_stats.html)]
 - **otutab_trim**  Trim OTU table to remove small counts, OTU and/or samples [[link](https://www.drive5.com/usearch/manual/cmd_otutab_trim.html)]
 - **otutab_xtalk**  Estimate and filter cross-talk in OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_xtalk.html)]
 - **qiimemap2otutab**  Convert QIIME map file to OTU table [[link](https://www.drive5.com/usearch/manual/cmd_qiimemap2otutab.html)]
 - **search_oligodb**  Search for matches to short nucleotide sequences, e.g. primers [[link](https://www.drive5.com/usearch/manual/cmd_search_oligodb.html)]
 - **search_pcr**  In-silico PCR, search for matches to pairs of primers in database [[link](https://www.drive5.com/usearch/manual/cmd_search_pcr.html)]
 - **search_pcr2**  In-silico PCR, search for matches to primer pair [[link](https://www.drive5.com/usearch/manual/cmd_search_pcr2.html)]
 - **search_phix**  Search for matches to PhiX sequence [[link](https://www.drive5.com/usearch/manual/cmd_search_phix.html)]
 - **sinaps**  Predict traits [[link](https://www.drive5.com/usearch/manual/cmd_sinaps.html)]
 - **sintax**  Predict taxonomy using SINTAX algorithm [[link](https://www.drive5.com/usearch/manual/cmd_sintax.html)]
 - **sintax_summary**  Generate summary report from sintax output [[link](https://www.drive5.com/usearch/manual/cmd_sintax_summary.html)]
 - **tabbed2otutab**  Convert read mapping file (read+OTU) to OTU table [[link](https://www.drive5.com/usearch/manual/cmd_tabbed2otutab.html)]
 - **uchime2_ref**  Chimera search using UCHIME2 algorithm [[link](https://www.drive5.com/usearch/manual/cmd_uchime2_ref.html)]
 - **uchime3_denovo**  Chimera search using UCHIME3 de-novo algorithm [[link](https://www.drive5.com/usearch/manual/cmd_uchime3_denovo.html)]
 - **unbias**  Correct abundance bias in OTU table [[link](https://www.drive5.com/usearch/manual/cmd_unbias.html)]
 - **unoise3**  Denoise amplicon reads [[link](https://www.drive5.com/usearch/manual/cmd_unoise3.html)]
 - **uparse_ref**  Classify sequences derived from mock community sample [[link](https://www.drive5.com/usearch/manual/cmd_uparse_ref.html)]
 
#### OTU table commands 
 - **alpha_div**  Calculate alpha diversity metric(s) from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_alpha_div.html)]
 - **alpha_div_rare**  Calculate alpha diversity metric(s) from OTU table with rarefaction [[link](https://www.drive5.com/usearch/manual/cmd_alpha_div_rare.html)]
 - **alpha_div_sig**  Statistical significance of alpha diversity correlation with metadata [[link](https://www.drive5.com/usearch/manual/cmd_alpha_div_sig.html)]
 - **beta_div**  Calculate beta diversity metric(s) from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_beta_div.html)]
 - **otutab**  Generate OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab.html)]
 - **otutab2biom**  Convert OTU table from tabbed to biom (json) format [[link](https://www.drive5.com/usearch/manual/cmd_otutab2biom.html)]
 - **otutab_binary**  Convert OTU table with counts to presence(1)/absence(0) [[link](https://www.drive5.com/usearch/manual/cmd_otutab_binary.html)]
 - **otutab_core**  Identify core microbiome in OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_core.html)]
 - **otutab_counts2freqs**  Convert counts to frequencies in OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_counts2freqs.html)]
 - **otutab_forest_classify**  Classify samples using random forest [[link](https://www.drive5.com/usearch/manual/cmd_otutab_forest_classify.html)]
 - **otutab_forest_train**  Train random forest on OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_forest_train.html)]
 - **otutab_group**  Sum subsets of samples in OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_group.html)]
 - **otutab_merge**  Merge two or more OTU tables [[link](https://www.drive5.com/usearch/manual/cmd_otutab_merge.html)]
 - **otutab_octave**  Generate octave plot visualizing OTU abundance distribution [[link](https://www.drive5.com/usearch/manual/cmd_otutab_octave.html)]
 - **otutab_otu_subset**  Extract subset of OTUs from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_otu_subset.html)]
 - **otutab_otus**  Extract OTU names from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_otus.html)]
 - **otutab_rare**  Rarefy OTU table so that samples have same number of reads [[link](https://www.drive5.com/usearch/manual/cmd_otutab_rare.html)]
 - **otutab_sample_subset**  Extract subset of samples from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_sample_subset.html)]
 - **otutab_samples**  Extract sample names from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_samples.html)]
 - **otutab_select**  Identify OTUs which are informative (predictive of metadata) [[link](https://www.drive5.com/usearch/manual/cmd_otutab_select.html)]
 - **otutab_sortotus**  Sort OTU table in order of decreasing OTU size [[link](https://www.drive5.com/usearch/manual/cmd_otutab_sortotus.html)]
 - **otutab_stats**  Report summary information about OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_stats.html)]
 - **otutab_trim**  Trim OTU table to remove small counts, OTU and/or samples [[link](https://www.drive5.com/usearch/manual/cmd_otutab_trim.html)]
 - **otutab_xtalk**  Estimate and filter cross-talk in OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_xtalk.html)]
 - **qiimemap2otutab**  Convert QIIME map file to OTU table [[link](https://www.drive5.com/usearch/manual/cmd_qiimemap2otutab.html)]
 - **tabbed2otutab**  Convert read mapping file (read+OTU) to OTU table [[link](https://www.drive5.com/usearch/manual/cmd_tabbed2otutab.html)]
 - **unbias**  Correct abundance bias in OTU table [[link](https://www.drive5.com/usearch/manual/cmd_unbias.html)]
 
#### Next-generation reads 
 - **fastq_chars**  Report frequencies of Q score ASCII characters in FASTQ file [[link](https://www.drive5.com/usearch/manual/cmd_fastq_chars.html)]
 - **fastq_eestats**  Report quality/e.e. per position for reads in FASTQ file [[link](https://www.drive5.com/usearch/manual/cmd_fastq_eestats.html)]
 - **fastq_eestats2**  Report number of reads retained at difference length and e.e. cutoffs [[link](https://www.drive5.com/usearch/manual/cmd_fastq_eestats2.html)]
 - **fastq_filter**  Filter reads in FASTQ file by e.e. and other criteria [[link](https://www.drive5.com/usearch/manual/cmd_fastq_filter.html)]
 - **fastq_join**  Concatenate forward (R1) and reverse (R2) paired reads [[link](https://www.drive5.com/usearch/manual/cmd_fastq_join.html)]
 - **fastq_mergepairs**  Assemble (merge) paired reads [[link](https://www.drive5.com/usearch/manual/cmd_fastq_mergepairs.html)]
 - **fastq_sra_splitpairs**  Recover paired reads from SRA interleaved or concatenated format [[link](https://www.drive5.com/usearch/manual/cmd_fastq_sra_splitpairs.html)]
 - **fastx2qiime**  Convert sample labels from usearch to QIIME format [[link](https://www.drive5.com/usearch/manual/cmd_fastx2qiime.html)]
 - **fastx_demux**  Assign reads to samples (demultiplex) [[link](https://www.drive5.com/usearch/manual/cmd_fastx_demux.html)]
 - **fastx_findorfs**  Identify ORFs in nucleotide sequences [[link](https://www.drive5.com/usearch/manual/cmd_fastx_findorfs.html)]
 - **fastx_get_sample_names**  Extract sample names from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_get_sample_names.html)]
 - **fastx_info**  Report summary information about a FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_info.html)]
 - **fastx_learn**  Estimate error rates from amplicon reads [[link](https://www.drive5.com/usearch/manual/cmd_fastx_learn.html)]
 - **fastx_subsample**  Extract random sub-sample from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_subsample.html)]
 - **fastx_syncpairs**  Sort forward and reverse reads into the same order [[link](https://www.drive5.com/usearch/manual/cmd_fastx_syncpairs.html)]
 - **fastx_trim_primer**  Remove primer-binding sequence from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_trim_primer.html)]
 - **fastx_truncate**  Truncate sequences in FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_truncate.html)]
 - **filter_lowc**  Filter low-complexity sequences from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_filter_lowc.html)]
 - **filter_phix**  Remove PhiX spike sequences from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_filter_phix.html)]
 - **search_oligodb**  Search for matches to short nucleotide sequences, e.g. primers [[link](https://www.drive5.com/usearch/manual/cmd_search_oligodb.html)]
 - **search_pcr**  In-silico PCR, search for matches to pairs of primers in database [[link](https://www.drive5.com/usearch/manual/cmd_search_pcr.html)]
 - **search_pcr2**  In-silico PCR, search for matches to primer pair [[link](https://www.drive5.com/usearch/manual/cmd_search_pcr2.html)]
 - **search_phix**  Search for matches to PhiX sequence [[link](https://www.drive5.com/usearch/manual/cmd_search_phix.html)]
 
#### Sequence database search 
 - :no_entry: **makeudb_sintax**  Make UDB database file for sintax [[link](https://www.drive5.com/usearch/manual/cmd_makeudb_sintax.html)]. See makeudb_usearch instead.
 - **makeudb_ublast**  Make UDB database file for ublast [[link](https://www.drive5.com/usearch/manual/cmd_makeudb_ublast.html)]
 - **makeudb_usearch**  Make UDB database file for usearch_global [[link](https://www.drive5.com/usearch/manual/cmd_makeudb_usearch.html)]
 - **search_exact**  Search for identical sequences [[link](https://www.drive5.com/usearch/manual/cmd_search_exact.html)]
 - **search_global**  Search database using global alignment without speed heuristics [[link](https://www.drive5.com/usearch/manual/cmd_search_global.html)]
 - **search_local**  Search database using locbal alignment without speed heuristics [[link](https://www.drive5.com/usearch/manual/cmd_search_local.html)]
 - **search_oligodb**  Search for matches to short nucleotide sequences, e.g. primers [[link](https://www.drive5.com/usearch/manual/cmd_search_oligodb.html)]
 - **search_pcr**  In-silico PCR, search for matches to pairs of primers in database [[link](https://www.drive5.com/usearch/manual/cmd_search_pcr.html)]
 - **search_pcr2**  In-silico PCR, search for matches to primer pair [[link](https://www.drive5.com/usearch/manual/cmd_search_pcr2.html)]
 - **search_peptidedb**  Search for matches to short peptide sequences [[link](https://www.drive5.com/usearch/manual/cmd_search_peptidedb.html)]
 - **search_phix**  Search for matches to PhiX sequence [[link](https://www.drive5.com/usearch/manual/cmd_search_phix.html)]
 - **ublast**  Fast database search using local alignment (much faster than BLAST) [[link](https://www.drive5.com/usearch/manual/cmd_ublast.html)]
 - **usearch_global**  Fast database search using global alignment [[link](https://www.drive5.com/usearch/manual/cmd_usearch_global.html)]
 - **usearch_local**  Fast database search using local alignment [[link](https://www.drive5.com/usearch/manual/cmd_usearch_local.html)]
 
#### Taxonomy commands 
 - **calc_lcr_probs**  Calculate Lowest Common Rank probabilities from dist. matrix with taxonomy [[link](https://www.drive5.com/usearch/manual/cmd_calc_lcr_probs.html)]
 - **makeudb_sintax**  Make UDB database file for sintax [[link](https://www.drive5.com/usearch/manual/cmd_makeudb_sintax.html)]
 - **nbc_tax**  Predict taxonomy using RDP Naive Bayesian Classifier algorithm [[link](https://www.drive5.com/usearch/manual/cmd_nbc_tax.html)]
 - **sintax**  Predict taxonomy using SINTAX algorithm [[link](https://www.drive5.com/usearch/manual/cmd_sintax.html)]
 - **sintax_summary**  Generate summary report from sintax output [[link](https://www.drive5.com/usearch/manual/cmd_sintax_summary.html)]
 
#### Tree commands 
 - **calc_distmx**  Calculate sparse distance matrix [[link](https://www.drive5.com/usearch/manual/cmd_calc_distmx.html)]
 - **cluster_tree**  Construct clusters from tree using distance cutoff [[link](https://www.drive5.com/usearch/manual/cmd_cluster_tree.html)]
 - **subtree**  Extract subtree under given node [[link](https://www.drive5.com/usearch/manual/cmd_subtree.html)]
 - **tree2distmx**  Calculate distance matrix implied by tree [[link](https://www.drive5.com/usearch/manual/cmd_tree2distmx.html)]
 - **tree_cvt**  Convert tree between tabbed and Newick formats [[link](https://www.drive5.com/usearch/manual/cmd_tree_cvt.html)]
 - **tree_subset**  Extract tree for subset of leaves [[link](https://www.drive5.com/usearch/manual/cmd_tree_subset.html)]
 
#### Labels and annotations 
 - **fastx2qiime**  Convert sample labels from usearch to QIIME format [[link](https://www.drive5.com/usearch/manual/cmd_fastx2qiime.html)]
 - **fastx_getlabels**  Extract sequence labels from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_getlabels.html)]
 - **fastx_relabel**  Re-label sequences in FASTx file with prefix plus sequential number [[link](https://www.drive5.com/usearch/manual/cmd_fastx_relabel.html)]
 - **fastx_strip_annots**  Remove usearch-style annotations (name=xxx) from FASTx file [[link](https://www.drive5.com/usearch/manual/cmd_fastx_strip_annots.html)]
 - **otutab_otus**  Extract OTU names from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_otus.html)]
 - **otutab_samples**  Extract sample names from OTU table [[link](https://www.drive5.com/usearch/manual/cmd_otutab_samples.html)]
   
 
