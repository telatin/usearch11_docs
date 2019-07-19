## Databases

Taxonomy and reference (FASTA) databases


### Check FASTA format

To be used as taxonomy reference databases, the FASTA files supplied to USEARCH should be formatted as:

`>Identifier;tax=d__Bacteria,p__Animalcules,o__Pilicoccales`

### Download USEARCH-formatted files

Providers:
 - [GreenGenes license terms](http://greengenes.secondgenome.com/)
 - [SILVA licence terms](https://www.arb-silva.de/silva-license-information/)
 - [UNITE licence terms](https://unite.ut.ee/repository.php)
 
```
DB_DIR='ref/'
mkdir -p "$DB_DIR"
cd "$DB_DIR"
wget -O "$DB_DIR/dp_16s_v16_sp.fa.gz" "https://www.drive5.com/sintax/rdp_16s_v16_sp.fa.gz"
wget -O "$DB_DIR/silva_16s_v123.fa.gz" "https://www.drive5.com/sintax/silva_16s_v123.fa.gz"
wget -O "$DB_DIR/gg_16s_13.5.fa.gz" "https://www.drive5.com/sintax/gg_16s_13.5.fa.gz"
wget -O "$DB_DIR/uchime_reference_dataset_28.06.2017.zip" "https://unite.ut.ee/sh_files/uchime_reference_dataset_28.06.2017.zip"
wget -O "$DB_DIR/unite.zip" "https://files.plutof.ut.ee/doi/BD/0A/BD0A35FEEF2DDF018E67970C7F7D56D02AB0C82A6396C14327FA0B5B1DBAF981.zip"
gunzip "$DB_DIR"/*.gz
unzip "$DB_DIR/unite.zip"
unzip "$DB_DIR/uchime_reference_dataset_28.06.2017.zip"
```
