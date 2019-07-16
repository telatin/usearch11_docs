# Cluster_fast 

Clusters sequences in a FASTA or FASTQ file using a variant of the UCLUST algorithm designed to maximize speed.

See: [Official documentation](https://www.drive5.com/usearch/manual/cmd_cluster_fast.html)

### Example
```
usearch -cluster_fast query.fasta -id 0.9 -centroids nr.fasta -uc clusters.uc
```

### Parameters
 - **-evalue** _real_: Maximum E-value. Required for most commands that use local alignments.
 - **-id** _real_: Minimum identity. Required for most commands that use global alignments.
 - **-query_cov** _real_: Fraction of the query sequence that is aligned, in the range 0.0 to 1.0.*
 - **-target_cov** _real_: Fraction of the target sequence that is aligned, in the range 0.0 to 1.0.*
 - **-maxdiffs** _int_: Maximum number of differences between the sequences, i.e. the maximum edit distance.
 - **-{min,max}sizeratio** _real_: Desired query\_size / target\_size ratio
 
 -  **-maxaccepts** and **-maxrejects**: Termination options
 
 See also:
  - [Indexing](indexing.md)
  - [Masking low complexity](masking.md)
 
 (*) With local alignments, this test is applied AFTER a local alignment is already created, so the effect is to reject local alignments that are too short, NOT to extend them further. With global alignments, columns containing terminal gaps are discarded before the test is applied.
