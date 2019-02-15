# Allele Counting Algorithm Proposal
## 12 Feb 18

### Objective

To remove all reads that could be affected by bisulfite conversion that would result in incorrect SNP reference/alt counting.


### Read examples

Four total cases-
  1. read is reverse, bisulfite conversion on reverse strand (C/T)
  2. read is reverse, bisulfite conversion on forward strand (A/G)
  3. read is forward, bisulfite conversion on forward strand (C/T)
  4. read is forward, bisulfite conversion on reverse strand (A/G)


In the following example reads from the data (from four separate SNPs) the observed read from the individual sequenced and the reference read it was mapped to are reported.

Differences between the reference sequence and the observed read are indicated in the reference sequence by a lower case letter.


```bash
# forward read example, A/G conversion
# chr4:154520218
obs-read: GAATTCAATAAACTAAAAATCAATCGTTTAAAAACAATAATCTTTCAAACTAACTAACTAATAATAAATACAACAAATAAAAATAAACAAACATAATAAAAAAAACAAAAATAATATTAAAAAAACTCAAA
ref-read: GgATTCAATgAgCTAAggATCAATCGTTTAAgAACAATAATCTTTCAAACcAACTAgCTAATAATAgATgCAACAAATAAAAATAAgCAAACATgATgAAAAAAgCAAAAgTAATATTAAAAgAgCTCAAA
```


```bash
# forward read example, C/T conversion
# chr4:16954354
obs-read: AATTTATATATTATTGTTAAGTATAGTGATTATATATTTTATATTTTATGTTAATTTTATAAAAATGTTAGTTATGTTTA GTGGAAATGTTATTTTTTTAGGGGAAGTGTTTATTTTATAGATAGGTTTTTA
ref-read: AATcTAcAcATTAcTGccAAGcAcAGTGAcTAcAcATTTcAcATTTcATGTTAATTTTAcAAAAATGTTAGccATGTcTAGGTGGAAATGcTATTcTTTcAGGGGAAGTGccTAccTTAcAGATAGGTTcccA
```


```bash
# reverse read example, A/G conversion
# chr4:154520218
obs-read: ATCGTTTAAAAACAATAATCTTTCAAACTAACTAACTAATAATAAATACAACAAATAAAAAATAAACAAACATAATAAAAAAACAAAAATAATATTAAAAAAACTCAAAAAATAACTAAAAAAATAAAATTTT
ref-read: ATCGTTTAAgAACAATAATCTTTCAAACcAACTAgCTAATAATAgATgCAACAAATAAAAAATAAgCAAACATgATgAAAAAAgCAAAAgTAATATTAAAAgAgCTCAAAgAATAACTgAAAgAATAAAATTTT
```


```bash
# reverse read example, C/T conversion
# chr4:57909378
obs-read: TTTATAATTTAAATGTTTTTTTTTATAAATGAGGTGATAAAATACGTTTTTTGGAAATTTTATTTGTTAATTTGAGTTATTGTATGTTTGGTGTTAGAAGAAAGGATTTAAATAGTTGTTATGAAAGT
ref-read: TTcAcAAccTAAATGTTTTcTTcTATAAATGAGGTGAcAAAAcACGTTTTTTGGAAAcTTcATcTGTTAAccTGAGccATTGcATGcTTGGTGcTAGAAGAAAGGAccTAAATAGcTGTcATGAAAGT
```

### SNP examples

```bash
# example 1 - would count as other
is_reverse: True
SNP: chr4:34393136
SNP ref: C
SNP alt: A
obs allele: T
mapped reference allele: c #lower case indicates that reference should be C, but observed read is different. In this case, the mapped reference allele is the same as the SNP ref.
```

```bash
# example 2 - would count as alt
is_reverse: True  
SNP: chr4 141036796
SNP ref: G
SNP alt: A
obs allele: A
mapped reference allele: g
```

```bash
# example 3 - would count as other
is reverse: False
SNP: chr4:74049709
SNP ref: T
SNP alt G
obs allele: A
mapped reference allele: t
```

``` bash
# example 4 - would count as other
is reverse: True
SNP: chr4:154366382
SNP ref: T
SNP alt: A
obs allele: C
mapped reference allele: t
```

```bash
# example 5 - would count as other
is reverse: False
SNP: chr4:169071642
SNP ref: T
SNP alt: G
obs allele: A
mapped reference allele: t
```

```bash
# example 6 - would count as other
is reverse: True
SNP: chr4:80897520
SNP ref: T
SNP alt: A
obs allele: C
mapped reference allele: t
```

```bash
#example 7 - count as other 
is reverse: False
SNP: chr5 158451151
SNP ref: A
SNP alt: C
obs allele: G
mapped reference: a
obs read: TGATTGTTGGAATAATAAAGATGAAAATAGTTATTTTTTTTTTTGAAGAAGTTTATGGTGGAGTGGGGGTATAGTAATTAGTTAATTAAGATGATAGTGTTTGTTAGATATTAGAGTTTTTAAGATAATTGATA
reference: TGATTGTTGGAATAAcAAAGATGAAAAcAGTTAcTccccTcccTGAAGAAGTTTATGGTGGAGTGGGGGTATAaTAATTAGcTAATcAAGATGATAGTGTTTGTcAGAcAccAGAGTcTcTAAGAcAATTGATA
```
### Proposal

Exclude all SNPs in (C:T, G:A, indel)

Let the snp of interest be SNP = (ref, alt).

Get the observed and reference read base pair at that snp: observed_read[SNP_POS], reference_read[SNP_POS]

Exclude SNPs where the reference_read[SNP_POS:, SNP_POS] == "CG" or "GC" (for reverse)

Exclude SNPs where reference_read[SNP_POS].upper() does not equal ref. This checks that the reference expected by the SNP is the same in the mapping genome I used.
