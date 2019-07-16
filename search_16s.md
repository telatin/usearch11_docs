# USEARCH 11

## Prediction of 16S sequences from contigs

Bit vector files are created by the udb2bitvec command. The recommended sequence database is the Greengenes 97% representative set. You should use the `-wordlength 13` option of **makeudb_usearch** to create the udb file, as shown below.

#### Prepare Database
```bash
wget "https://github.com/biocore/qiime-default-reference/raw/master/qiime_default_reference/gg_13_8_otus/rep_set/97_otus.fasta.gz"

gunzip 97_otus.fasta.gz
mv 97_otus.fasta gg97.fa

usearch -makeudb_usearch gg97.fa -wordlength 13 -output gg97.udb
usearch -udb2bitvec gg97.udb -output gg97.bitvec
```

#### Predict 16S from contigs
```bash
DB=/path/to/gg97.bitvec
OUT=$(echo $CONTIGS | rev | cut -f 2- -d . | rev )

usearch -search_16s $CONTIGS -bitvec $DB -fastaout ${BASENAME}.16s.fa -tabbedout ${BASENAME}.16Sresults.txt -fragout ${BASENAME}.part16s.fa -hitsout ${BASENAME}.16Shits.fa
```

### Output
 - **FASTA** file with the predicted genes (-fastaout, -hitsout, -fragout)
 - **Tabular file** with the predictions (-tabbedout):
```
ctg1     query   length=4584796  wins=7  genes=7 frags=0
ctg1     win     strand=+        lo=3633812      hi=3636839      un=947956       len=3028        genes=1 starts=1(458)/0 ends=1(1929)/0
ctg1     win     strand=+        lo=3727139      hi=3730189      un=854606       len=3051        genes=1 starts=1(475)/0 ends=1(1946)/0
ctg1     win     strand=+        lo=3860285      hi=3863331      un=721464       len=3047        genes=1 starts=1(479)/0 ends=1(1950)/0
ctg1     win     strand=+        lo=3901582      hi=3904664      un=680131       len=3083        genes=1 starts=1(507)/0 ends=1(1978)/0
ctg1     win     strand=+        lo=4558872      hi=4561911      un=22884        len=3040        genes=1 starts=1(467)/0 ends=1(1941)/0
ctg1     win     strand=-        lo=1461040      hi=1464066      un=3120729      len=3027        genes=1 starts=1(456)/0 ends=1(1928)/0
ctg1     win     strand=-        lo=2156711      hi=2159775      un=2425020      len=3065        genes=1 starts=1(484)/0 ends=1(1956)/0
ctg1     gene    strand=+        lo=3634271      hi=3635762      len=1492        start=GTTTGATCATGGC/0   end=AGTCGTAACAAGGTAACCGTA/0
ctg1     gene    strand=+        lo=3727615      hi=3729106      len=1492        start=GTTTGATCATGGC/0   end=AGTCGTAACAAGGTAACCGTA/0
ctg1     gene    strand=+        lo=3860765      hi=3862256      len=1492        start=GTTTGATCATGGC/0   end=AGTCGTAACAAGGTAACCGTA/0
ctg1     gene    strand=+        lo=3902090      hi=3903581      len=1492        start=GTTTGATCATGGC/0   end=AGTCGTAACAAGGTAACCGTA/0
ctg1     gene    strand=+        lo=4559340      hi=4560834      len=1495        start=GTTTGATCATGGC/0   end=AGTCGTAACAAGGTAACCGTA/0
ctg1     gene    strand=-        lo=1461497      hi=1462989      len=1493        start=GTTTGATCATGGC/0   end=AGTCGTAACAAGGTAACCGTA/0
ctg1     gene    strand=-        lo=2157196      hi=2158688      len=1493        start=GTTTGATCATGGC/0   end=AGTCGTAACAAGGTAACCGTA/0
```

## Sample script
```bash
#!/bin/bash
DB=$PATH_TO/gg97.bitvec
usearch=usearch_11

if [ -e "$1" ]; then
  CONTIGS=$1
  if [ "NO$2" == "NO" ]; then
    OUT=$(readlink -f "$CONTIGS" | rev | cut -f 2-100 -d . | rev)
  else
    OUT="$2"
  fi

  echo "Output fasta: $OUT.16S.fa"
  echo "Fragments fa: $OUT.16Spart.fa"
  echo "Tabular file: $OUT.16S.txt"
  echo "Fasta hits:   $OUT.16S.hits"

  $usearch -search_16s "$CONTIGS" -bitvec "$DB" -fastaout "$OUT.16S.fa" \
     -fragout "$OUT.16Spart.fa" -tabbedout "$OUT.16S.txt" -hitsout "$OUT.16S.hits"

else
  echo "USAGE: Contigs.fa OutPrefix"
  echo "Specify contigs file"

fi
```
