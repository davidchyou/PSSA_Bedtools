# PSSA_Bedtools

**What it does?**

Detecting adverse drug reactions (ADR) is very important in minimising drug-related harms especially in older individuals aged 65 years and above. In the Big-Data era, many computational algorithms had been developed to predict ADR of the various drugs. The prescription sequence-symmetry analysis algorithm, [PSSA](https://www.jstor.org/stable/pdf/3702146.pdf), is a very powerful and computationally simple algorithm to predict ADR, and is scalable to big, national level healthcare dataset with more than 1 billion rows of data. It requires only minimal amount of data: prescription dates and drugs prescribed, grouped by individuals, and we only need prescription dataset. The idea is that if drug A causes medical condition X, and X can be treated with drug B, then, the number of prescription sequence A -> B is expected to be greater than the number of prescription sequence B -> A. The statistical quantity "sequence ratio" (SR), n(AB)/n(BA), is the main ADR predictor in PSSA, although log(SR) is preferred because the standard error of it is easy to compute. 

To avoid confounding, the principle of "new user design" applies and we only consider first-time prescriptions of all drugs of interests. The prescription dataset is therefore processed into a "drug-initiation dataset". For each individual, first-time prescription of every drug the individual was used was extracted. To mitigate confounding further, sequence A -> B are two-drug sequences, where drug B is the next drug initiated after initiating drug A, and no additional drugs were initiated in between. Furthermore, we only count drug initiation events that are less than 100 days apart.

**How the genomic application Bedtools can be applied to solve this pharmacoepidemiological problem?** 

In the first step of the PSSA algorithm, lists of two-drug initiation sequences were extracted from all individuals. Because we want initiation sequence A -> B to be genuinely A -> B, not either A -> B or A -> C -> B (and so on), B is the closest event from A forwards in time. Therefore, the computation is analogous to the following bioinformatic data-extraction problem: for each gene in the genome, within the same genomic contig, find the nearest gene downstream. The genomic analysis application [Bedtools](https://bedtools.readthedocs.io/en/latest/) has a sub-application ["closest"](https://bedtools.readthedocs.io/en/latest/content/tools/closest.html) that does exactly this computation, and therefore we can use it to perform PSSA. All we need to do is to generate a BED6-formatted drug-initiation dataset, which is bioinformatician's day-to-day practice.

In genomics, each genomic features (genes, intron, exon, repeat-region etc.) are recorded in a BED6-formatted file, that has 6 fields separated by tabs: contig-accession, feature-start, feature-end, feature-name, score (e.g. statistical confidence if the feature predicted computationally), and the genomic strand (+ or -). The score field can also be used for other attributes. In pharmacoepidemiology, we can use BED6-formatted file to record prescription or medical-event data by recording individual-ID, event-start, event-end, event-name or drug-name, and ICD-code or chemical-ID. Since we are dealing with a timeline not genomic sequences, the last field is always "+" as time always goes forwards.

In PSSA, the BED6-formatted drug-initiation dataset looks like:

    Person-1  1324  1324  Drug-A  1641  +
    Person-2  2113  2113  Drug-B  3877  +
    Person-3  4218  4218  Drug-C  8594  +
    
Because drug-initiations are one-day event, so start = end, and Bedtools allows this.

**Usage**

The Shell script reads in a BED6-formatted drug-initiation dataset as above-mentioned. The drug-initiation dataset need to be sorted by individual-ID, and event-start. On a Linux console, this can be done by the command:

    sort -k1,1 -k2,2n <input.bed> > <input.sorted.bed>
    
Then, to run the script:

    sh PSSABed6.sh <input.sorted.bed>
    sh PSSABed6.sh <input.sorted.bed> > <seq_ratios.txt>

Log sequence ratios of all two-drug initiation sequences of all drugs in the dataset will be tabulated with 95% confidence intervals and supporting information. By default, all statistical data will be displayed on the console, but can be directed into a file.  

To avoid dividing-by-zero or computing log(0), observed n(AB) and n(BA) were both incremented by 1 (i.e. a pseudocount) when computing log sequence ratios and the standard errors, but the original n(AB) and n(BA) is included in the output.

**Output**

The output is a tab-delimited table showing (from left to right) the following:

1. Chemical ID of Drug A (The first drug initiated).
2. Chemical ID of Drug B (The drug initiated after Drug A).
3. Chemical name of Drug A.
4. Chemical name of Drug B.
5. Frequence of sequence A -> B (n(AB))
6. Frequence of sequence B -> A (n(BA))
7. Frequence of sequence A -> B with a pseudocount (n(AB) + 1)
8. Frequence of sequence B -> A with a pseudocount (n(BA) + 1)
9. log(sequence ratio with pseudocounts) = log((n(AB) + 1) / (n(BA) + 1))
10. Lower-bound of the 95% confidence interval.
11. Upper-bound of the 95% confidence interval.
12. Is it an ADR signal? (1 = yes,0 = no)

A log(SR) is mathematically equivalent to the relative risk calculated by univariate conditional logistic regression, which can be expressed in closed form as log(n(BA)/n(BA)), and the standard error is therefore sqrt(1/n(AB) + 1/n(BA)), or sqrt(1/(n(AB)+1) + 1/(n(BA)+1)) with pseudocounts. 

We observed an ADR signal if log(SR) is significantly above zero statistically, in which case the lower-bound of the 95% confidence interval is greater than zero.
