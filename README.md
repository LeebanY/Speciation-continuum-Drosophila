# Speciation-continuum-Drosophila
Inferring demographic history and effective migration genome-wide in species pairs across Drosophila.

#### The analysis pipeline is as follows:

01. Sample Metdata: Random selection (if possible) of publically-available sequences to use from NCBI. This contains filtering of metadata to ensure only WGS datasets containing a single individual were used, as well as data that 

2. Genome annotation via BRAKER2 with D. melangoaster proteins.

3. Read filtering, mapping and variant calling for all pairs.

4. Preprocessing of vcfs.

5. Summary statistic calculation (coverage, hetA, hetB, hetAB, absolute genetic divergence).

6. Filtering intronic intervals using BEDTOOLS.

7. Extraction of distribution of pairwise divergence.

8. Modelling via Mathematica notebooks.

9. Analysis of modelling results in R.
