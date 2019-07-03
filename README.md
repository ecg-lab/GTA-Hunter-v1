# GTA_Hunter
  It is a bioinformatics tool that classifies homologs from RcGTA structural gene as true GTAs or viruses. 

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

**GTA_Hunter.py is the main file, and is used through command line using the following commands:**
```
usage: GTA_Hunter_blast.py [-h] [-b] [-g GTA] [-v VIRUS] [-q QUERIES]
                           [-o OUTDIR] [-k [KMER]] [-p [PSEAAC]] [-y] [-m]
                           [-O] [-f FOLDER] [-W] [-w WEIGHT WEIGHT]
                           [-t CLUSTER_TYPE] [-d DIST] [-c C] [-x [XVAL]]
                           [-e KERNEL KERNEL] [-s]

Gene Classification Using SVM.

optional arguments:
  -h, --help            show this help message and exit
  -b, --blast           Run blast to identify gta homologs.
  -g GTA, --GTA GTA     The .faa or .fna training file for GTA genes.
  -v VIRUS, --virus VIRUS
                        The .faa or .fna training file for viral genes.
  -q QUERIES, --queries QUERIES
                        The .faa or .fna query file to be classified.
  -o OUTDIR, --outdir OUTDIR
                        The folder path in which to store output.
  -k [KMER], --kmer [KMER]
                        The kmer size needed for feature generation
                        (default=4).
  -p [PSEAAC], --pseaac [PSEAAC]
                        Expand feature set to include pseudo amino acid
                        composition. Specify lambda (default= None). Weight =
                        0.05.
  -y, --physico         Expand feature set to include physicochemical
                        composition.
  -m, --min             Print bare minimum results.
  -O, --optimal         Pick the optimal parameters for the GTA gene homolog
                        classification
  -f FOLDER, --folder FOLDER
                        Feed GTA_Hunter a folder from which to run the program
                        on all
  -W                    Weight training set if desired. Distance files will be
                        supplied automatically.
  -w WEIGHT WEIGHT, --weight WEIGHT WEIGHT
                        Allows for weighting of training set. Will need to
                        specify the two pairwise distance files needed for
                        weighting (GTA first, then virus).
  -t CLUSTER_TYPE, --cluster_type CLUSTER_TYPE
                        Specify 'farthest' or 'nearest' neighbors clustering
                        (default='farthest').
  -d DIST, --dist DIST  Specify the cutoff distance for clustering in the
                        weighting scheme (default=0.01).
  -c C, --soft_margin C
                        The soft margin for the SVM (default=1.0).
  -x [XVAL], --xval [XVAL]
                        Performs cross validation of training set. Specify
                        folds over 10 repetitions (default=5).
  -e KERNEL KERNEL, --kernel KERNEL KERNEL
                        Specify kernel to be used and sigma if applicable
                        (i.e. gaussian) (default='linear', 0).
  -s, --svs             Show support vectors.

```
## Example run
  Try to make a test run by classifying homologs of RcGTA phage portal protein using k-mer of size 3
```
python GTA_Hunter.py -g data/training/gta/3_gta.faa -v data/training/viral/3_viral.faa -c 100 -k 3 -q example_run/g3_example.faa
```
  As an output, the first two genes will be classified as "GTAs" and other three as "viruses".
  If you want to identify RcGTA-like genes from the whole proteome, the 'blast' function can be used. Try the following example run for Paracoccus marcussi (note: create 'example_outdir' folder):
```
python GTA_Hunter.py -b -f example_blast/ -o example_outdir/ -O
```
  The genes will be classified using optimal parameters that were detected via cross validation. The results will be written to the file 'result_(file_name).out'.

## From genes to cluster
  Detected RcGTA-like genes can be clustered using DBSCAN via additional script, which is located at the "Clustering genes" folder.
**Clustering_genes.py is used through command line via the following commands:**  
```  
usage: Clustering_genes.py [-h] --data DATA --feature FEATURE [-s SIZE]
                           [-e DIST]

Clustering of genes via DBSCAN to get RcGTA-like clusters

optional arguments:
  -h, --help            show this help message and exit
  --data DATA           The text files with the chromosome ID first and genes
                        to cluster.
  --feature FEATURE     Supplied feature table from the NCBI.
  -s SIZE, --Cluster SIZE
                        The minimum cluster size to consider as RcGTA-like
                        cluster (default=6).
  -e DIST, --Max DIST   The maximum spatial distance to cluster genes
                        (default=8000)
```
  The test run can be done using the following example:
```
python Clustering_genes.py --data example_cl/data.txt --feature example_cl/GCF_006151785.1_feature.txt
```
  As an output, test run will show that RcGTA-like cluster contains 11 genes.
