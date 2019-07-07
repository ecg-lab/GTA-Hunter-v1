# GTA_Hunter
  GTA_Hunter is a 'support-vector machine' (SVM) classifier that distinguishes gene transfer agent (GTA) genes from their viral homologs. 
  
  The classifier was developed to detect GTA genes that are homologous to 11 of the genes that encode the GTA in *Rhodobacter capsulatus* (RcGTA) (see Kogay et al, 2019). However, the classifier can be trained to detect other RcGTA genes, other GTAs, or even other virus-like elements, as long as there is training data available and the two classes have distinct amino acid composition. 

## System Requirements
 `Python v3.5.1` 

## Python packages
- NumPy v. 1.13.1
- CVXOPT v. 1.1.8
- matplotlib v. 2.1.0
- pandas v. 0.21.0
- Biopython v. 1.69

## Third party requirements
- BLAST v. 2.2.31+

**GTA_Hunter.py is a command-line program with the following options:**

```
usage: GTA_Hunter_blast.py [-h] [-b] [-g GTA] [-v VIRUS] [-q QUERIES]
                           [-o OUTDIR] [-k [KMER]] [-p [PSEAAC]] [-y] [-m]
                           [-O] [-f FOLDER] [-W] [-w WEIGHT WEIGHT]
                           [-t CLUSTER_TYPE] [-d DIST] [-c C] [-x [XVAL]]
                           [-e KERNEL KERNEL] [-s]

Gene Classification Using SVM.

optional arguments:
  -h, --help            Show this help message and exit.
  -b, --blast           Run BLASTP search to first identify GTA homologs for the classification.
  -g GTA, --GTA GTA     The FASTA-formatted (.faa or .fna) "true GTA" sequences used for training.
  -v VIRUS, --virus VIRUS
                        The FASTA-formatted (.faa or .fna) "true viruses" sequences used for training.
  -q QUERIES, --queries QUERIES
                        The FASTA-formatted (.faa or .fna) sequences to be classified.
  -o OUTDIR, --outdir OUTDIR
                        The folder path in which to store the output.
  -k [KMER], --kmer [KMER]
                        The size of k-mer (default=4).
  -p [PSEAAC], --pseaac [PSEAAC]
                        Expand feature set to include pseudo amino acid
                        composition (default=3). As the parameter, specify value of lambda. Weight = 0.05 (after Chou 2001). 
  -y, --physico         Expand feature set to include physicochemical
                        properties of amino acids.
  -m, --min             Print bare minimum of the results.
  -O, --optimal         Use the optimal parameters for the RcGTA gene homolog
                        classification as listed in Table 2 in Kogay et al 2019.
  -f FOLDER, --folder FOLDER
                        Provide a folder with one or multiple proteomes (*.faa files).
  -W                    Weigh the training set. Distance files will be
                        supplied automatically.
  -w WEIGHT WEIGHT, --weight WEIGHT WEIGHT
                        Weigh the training set. Specify the two pairwise distance files needed for
                        weighting (first file for GTAs, second file for viruses).
  -t CLUSTER_TYPE, --cluster_type CLUSTER_TYPE
                        Specify 'farthest' or 'nearest' neighbors clustering
                        (default='farthest').
  -d DIST, --dist DIST  Specify the cutoff distance for clustering in the
                        weighting scheme (default=0.01).
  -c C, --soft_margin C
                        The soft margin for the SVM (default=1.0).
  -x [XVAL], --xval [XVAL]
                        To perform cross validation of the training set. Specify
                        folds over 10 repetitions (default=5).
  -e KERNEL KERNEL, --kernel KERNEL KERNEL
                        Specify kernel to be used and sigma if applicable
                        (i.e. Gaussian) (default='linear', 0).
  -s, --svs             Show support vectors.

```

## How to use it - Example 1: Classify homologs of RcGTA-like gene using the provided training data
In this example, we will classify five homologs of the RcGTA portal protein (g3) using k-mer of size 3 as a feature and setting raw softness of the SVM margin to 100. The sequences to classify are in FASTA format in one file, and are located in the "example_run" folder; the training data is in "data/training" folder. 

```
python GTA_Hunter.py -g data/training/gta/3_gta.faa -v data/training/viral/3_viral.faa -c 10000 -w data/training/gta/3_gta.dist data/training/viral/3_viral.dist -t 0.02 -k 3 -q example_run/g3_example.faa
```

In the output, the first two sequences will be classified as a "GTA" and the remaining three as a "virus".

```
Gene                                                                                           Score          Classification    GTA Gene
>WP_100283140.1 phage portal protein [Sphingomonas sp. Cra20]                                   -0.824313      GTA               3
>WP_100367008.1 phage portal protein [Yoonia maricola]                                          -0.947548      GTA               3
>WP_102855116.1 phage portal protein [Phaeobacter inhibens]                                     0.873095       virus             3
>WP_121219731.1 phage portal protein [Oceanibaculum indicum]                                    0.965013       virus             3
>WP_128262196.1 phage portal protein [Bradyrhizobium yuanmingense]                              0.702121       virus             3
time to run: 10.004424810409546
```

## How to use it - Example 2: Cross-validate classifier for a gene

In this example, we will identify the accuracy of the classifier using five-fold cross-validation.
```
python GTA_Hunter.py -g data/training/gta/4_gta.faa -v data/training/viral/4_viral.faa -c 10000 -w data/training/gta/4_gta.dist data/training/viral/4_viral.dist -t 0.02 -k 3 -p -x 5
```
In the output, on average all sequences will be classified correctly, indicating that tested classifer has 100% of accuracy.

```
We correctly classified (on average) 62.00/62 GTA and 465.00/465 Viral genes.
time to run: 61.76529502868652
```

## How to use it - Example 3: Find homologs of RcGTA-like gene in a genome and classify them using the provided training data and the most optimal SVM parameters
In this example, we will first identify homologs of the RcGTA genes among the encoded proteins in a *Paracoccus marcussi* genome using a BLASTP search (Altschul et al. 1997), and then classify them using training data provided in data/taining folder and the optimal parameters that were detected via cross validation (see Table 2 in Kogay et al. 2019). The input file is FASTA-formatted file with amino acid sequences of encoded proteins (*.faa format provided by GenBank@NCBI) and is located in the "example_blast" folder. 

Note: please create 'example_outdir' folder before running the program.

```
python GTA_Hunter.py -b -f example_blast/ -o example_outdir/ -O
```
The results will be written to the file 'result_(file_name).out' in the output directory.
The BLASTP search will identify 15 RcGTA homologs and all of them will be classifier as RcGTA-like genes (see 'example_output' for details)

## Are the detected GTA-like genes located in the same neighborhood in a genome?

`GTA_Hunter.py` only evaluates one gene at a time, even if queries come from the same genome or even when the whole genome was scanned. To see if the detected "GTA" genes are found in the same neighborhood in a genome, one can run an additional script (`Clustering_genes.py`; located in the "Clustering genes" folder). The script uses the DBSCAN algorithm (Ester et al. 1996).

**Clustering_genes.py is a command-line program with the following options:**  

```  
usage: Clustering_genes.py [-h] --data DATA --feature FEATURE [-s SIZE]
                           [-e DIST]

Clustering of genes via DBSCAN to get RcGTA-like clusters

optional arguments:
  -h, --help            Show this help message and exit
  --data DATA           File that containes one line with space or tab-separated entries. First entry is the genome or contig ID, 
						the rest are RefSeq IDs of the encoded proteins to cluster
  
  --feature FEATURE     Genome feature table in the NCBI feature table file format.
  -s SIZE, --Cluster SIZE
                        The minimum cluster size to consider as RcGTA-like
                        cluster (default=6).
  -e DIST, --Max DIST   The maximum spatial distance to cluster genes
                        (default=8000)

```

  For example, let's examine 12 RcGTA-like homologs detected in the *Rhodobacter sphaeroides strain HJ* genome from the example 3 above. To call a region an "RcGTA-like cluster", we will require that at least 6 of the inputted genes have no more than 8,000 base pairs between the adjacent genes. 
  
```
python Clustering_genes.py --data example_cl/data.txt --feature example_cl/Rhodobacter_sphae_HJ.txt -s 6 -e 8000
```
  The program will report that 11 out of the 17 provided RcGTA-like genes are in the same neighborhood in the genome, and therefore can be designated as an "RcGTA-like cluster".

```
NZ_CP036419.1 has RcGTA-like cluster; The cluster size is 11 genes
WP_137457888.1  WP_137457889.1  WP_137457890.1  WP_115474024.1  WP_009564481.1  WP_137457892.1  WP_137457893.1  WP_101328016.1  WP_115474016.1  WP_137457895.1  WP_137457896.1
```
  
  All example outputs are stored at the folder "example_output".

## References

- Altschul SF, et al. 1997. Gapped BLAST and PSI-BLAST: a new generation of protein database search programs. Nucleic Acids Res 25: 3389-3402
- Chou KC. 2001. Prediction of protein cellular attributes using pseudo‐amino acid composition. Proteins 43: 246-255.
- Ester M, Kriegel H-P, Sander J, Xu X editors. Kdd. 1996
- Kogay et al 2019
