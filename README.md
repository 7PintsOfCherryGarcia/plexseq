# plexseq
Demultiplex fastq files (accepts gzip compression) acording to barcode sequences present on reads.

plexseq takes fastq formated reads and a barcode file to separate sequencing
reads according to said barcodes.

Data may be provided in several ways:

	1) single fastqfile (single end reads, -r option)
	
	2) two fastq files (R1 R2 paired reads, -r and --r2 options)
	
	3) interleaved fastqfile (R1 followed by R2, --ri option)
	
	4) fastq data from STDIN (single or interleaved paired end, --s option)

A barcode file must be provided and data should be structured as follows:

tab separated file with at least two columns

BARCODE_ID|BARCODE_SEQUENCE|BARCODE_SEQUENCE2(optional)
----------|----------------|---------------------------

Example:

\#ID|Seq|R2(optional)
----|----|----
A1F1|CAGTCGAAT|TCACAGATG
A1F2|TCAGTCGAA|TCACAGATG
A1F3|CTCAGTCGA|TGCAGATGA
A1F4|CAGAAAAGT|TTGTTCTCC
A1F5|TCAGAAAAG|TGTTCTCCA
A1F6|CTCAGAAAA|TGTTCTCCA
A1F7|TCAGAAAAG|CTTGTTCTC
A1F8|TCAGTCGAA|CTCACAGAT

lines starting with # are ignored

plexseq works by looking the barcode sequence in specific parts of the sequence
read. For now plexseq looks for barcode sequence 1 in the first 9 bases of the
5' end of read1. If a second barcode sequence is provided, plexseq will look
for this sequence in the first 9 bases of the 5' end of read2 in the case of
paired reads. 

plexseq looks for exact matches between reads and barcodes.

This limited functionality is just for the moment. Future versions will automatically
set the number of bases to look in each read based on the length of the barcode,
and lookups with mismatches will be introduced.

Usage: 

python demultiplex.py -r1 "read_file" -i "indexfile" -o "outdir" [--r2 "read2_file" --ri "interleaved_read_file"] 

Or: 

cat seqfile | python demultiplex.py -i "indexfile" -o "outdir" --s   


# Issues

For the moment plexseq can only separate sequences which exactly match barcode sequences

All barcode sequences have to be of the same length

Can not write gzip compressed files

Barcode sequences have to be located at the 5' end of reads

Only fastq files can be demultiplexed
