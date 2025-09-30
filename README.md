# MIXPRS
MIXPRS is a data fission-based multi-population PRS integration framework designed to effectively combine PRS derived from multiple populations and methods. This framework is described in detail in:

- **Xu, L., Dong, Y., Zeng, X., Bian, Z., Zhou, G., Guan, L., & Zhao, H. (2025). Almost Free Enhancement of Multi-Population PRS: From Data-Fission to Pseudo-GWAS Subsampling. bioRxiv, 2025-06.**

The MIXPRS pipeline requires **only GWAS summary statistics and LD reference panels**, and involves three main steps (Figure 1):

* **Step1: GWAS subsampling in MIXPRS**  
  Generate independent subsampled training and tuning GWAS datasets.

* **Step2: Estimation of MIXPRS combination weights for different PRS methods**  
  Estimate combining weights using two component methods, **JointPRS-auto** and **SDPRX**, each requiring only GWAS summary statistics and LD reference panels.

* **Step3: Derivation of MIXPRS**  
  Calculate the integrated MIXPRS using weights estimated in Step 2.

Detailed implementations of these component methods are available in their respective repositories:

* [JointPRS-auto](https://github.com/LeqiXu/JointPRS)
* [SDPRX](https://github.com/eldronzhou/SDPRX)

<p align="center">
  <img src="https://github.com/user-attachments/files/20983662/Figure1.pdf" alt="MIXPRS Workflow"/>
</p>

For full details of codes, see the [MIXPRS_analysis repository](https://github.com/LeqiXu/MIXPRS_analysis).

## Getting Started
In this section, we provide detailed, step-by-step instructions for implementing MIXPRS. Please replace all placeholders with the appropriate paths and filenames specific to your computing environment.

### 1. MIXPRS Installation
For the first time, you need to use the following code to install MIXPRS (using [conda](https://docs.conda.io/) for environment management):
```
git clone https://github.com/LeqiXu/MIXPRS.git
cd MIXPRS
conda env create -f environment.yml
conda activate MIXPRS
```
> **Note:** If installing the packages from `environment.yml` is too slow on your cluster, you can alternatively run:
>
> ```bash
> git clone https://github.com/LeqiXu/MIXPRS.git
> cd MIXPRS
> conda create -n MIXPRS python=3.8
> conda activate MIXPRS
> pip install -r requirements.txt
> ```
After this, you only need to use
```
conda activate MIXPRS
```

### 2. LD Reference Panel Download
We use reference panels from [PRS-CSx](https://github.com/getian107/PRScsx#getting-started) and you can follow their instructions to download them. It is strongly recommended to create two subfolders within your reference directory:

* **1KG**: This subfolder should contain LD reference panels constructed using the 1000 Genomes Project phase 3 samples and the corresponding SNP information file. Specifically, include:

  * `snpinfo_mult_1kg_hm3`
  * Population-specific directories (`ldblk_1kg_afr`, `ldblk_1kg_amr`, `ldblk_1kg_eas`, `ldblk_1kg_eur`, `ldblk_1kg_sas`), each containing:

    * `ldblk_1kg_chr1.hdf5` through `ldblk_1kg_chr22.hdf5`
    * `snpinfo_1kg_hm3`

* **UKBB**: This subfolder should contain LD reference panels constructed using the UK Biobank data and the corresponding SNP information file. Specifically, include:

  * `snpinfo_mult_ukbb_hm3`
  * Population-specific directories (`ldblk_ukbb_afr`, `ldblk_ukbb_amr`, `ldblk_ukbb_eas`, `ldblk_ukbb_eur`, `ldblk_ukbb_sas`), each containing:

    * `ldblk_ukbb_chr1.hdf5` through `ldblk_ukbb_chr22.hdf5`
    * `snpinfo_ukbb_hm3`

Then place the downloaded LD reference panels and the SNP information file into their corresponding subfolders.

### 3. Pruned SNP List Preparation
In the `MIXPRS` repository, we provide precomputed pruned SNP lists for five populations (EUR, EAS, AFR, SAS, and AMR), derived from genotype data of the [Phase-3 1000 Genomes Project](https://www.internationalgenome.org/data). These SNP lists were generated using PLINK2 with parameters: window size = 250, step size = 5, and correlation threshold $r^2 = 0.5$.

To replicate our results or generate customized pruned SNP lists, use the following PLINK2 commands:
```bash
module load PLINK/2

i=1
pval=1
r2=0.5
wc=250

for pop in EUR EAS AFR SAS AMR; do

plink2 --bfile ./1000g_phase3_data/geno_data/${pop} \
       --indep-pairwise ${wc} 5 ${r2} \
       --out ./snplist/${pop}_prune_pval${pval}_r2${r2}_wc${wc}_${i}
done
```
You can adjust these parameters or apply the commands to other genotype datasets according to your analytical needs.

### 4. Summary Statistics Preparation
We require GWAS summary statistics provided as input to MIXPRS to follow the format below (including the header line):

```
SNP         A1  A2  BETA        SE          Z                P      N
rs3934834   T   C   0.00492590  0.00761858  0.646564057869   0.518  100061
rs3766192   C   T   0.00144597  0.00580197  0.249220523374   0.803  100953
rs9442372   A   G   0.00261076  0.00586576  0.445084694907   0.656  98059
...
```

Here
* `SNP`: SNP rsID.
*  `A1`: Effect allele.
*  `A2`: Alternative allele.
*  `BETA`: Effect size of allele A1, which is only used to determine the direction of an association.
*  `SE`: Standard error of the effect size.
*  `Z`: Z-score (BETA divided by SE).
*  `P`: P-value of the effect size, which is used to calculate the standardized effect size.
*  `N`: GWAS sample size.

### 5. MIXPRS Implementation
#### Step 0: MIXPRS data preparation
Before running MIXPRS, ensure these datasets are ready (see previous sections):
* **LD Reference Panel (Section 2).**
* **Pruned SNP List (Section 3).**
* **GWAS Summary Statistics (Section 4).**

Example assumptions:
* 1KG LD panel used.
* Precomputed pruned SNP lists from MIXPRS repository (`snplist`).
* Formatted GWAS summary statistics `${trait}_${pop}_MIXPRS_sumstat.txt`. We provide an example summary statistics dataset for EAS HDL (500 SNPs on chromosome 1) within the cloned MIXPRS repository (`example_data` folder), obtained from [GLGC](https://csg.sph.umich.edu/willer/public/glgc-lipids2021/).

#### Step1: GWAS subsampling in MIXPRS
This step partitions a single original GWAS dataset into independent subsampled training and tuning GWAS datasets using data fission principles for the target population (`pop`).

The default parameters for subsampling are:
* Requires a pruned SNP list (`--prune_snplist`) provided in the MIXPRS repository (`snplist` folder).
* `indep_approx=TRUE`: uses an identity covariance matrix for the pruned SNPs instead of LD reference panels (the `--ref_dir` flag is still required).
* `train_tune_ratio=3`: divides the original GWAS dataset into subsampled training GWAS (3/4 of samples) and subsampled tuning GWAS (1/4 of samples).
* `repeat=4`: generates four independent subsampled training-tuning GWAS pairs for robust evaluation.

Run the following command to subsample training and tuning GWAS:
```bash
conda activate MIXPRS

python ${MIXPRS_path}/MIX_subsample2.py \
  --ref_dir=${ref_data_path}/1KG \
  --sst_file=${summary_stat_path}/${trait}_${pop}_MIXPRS_sumstat.txt \
  --pop=${pop} \
  --prune_snplist=${MIXPRS_path}/snplist/${pop}_prune_pval1_r20.5_wc250_1.snplist \
  --indep_approx=TRUE \
  --train_tune_ratio=3 \
  --repeat=4 \
  --out_dir=${output_path}/subsample/clean \
  --out_name=${trait}_prune_snplist_1
```
This command generates four pairs of independent subsampled training and tuning GWAS summary statistics files, named as follows (`repeat=1,2,3,4`):  
* **Subsampled training GWAS**:
  `${output_path}/subsample/clean/${trait}_prune_snplist_1_${pop}_train_GWAS_approxTRUE_ratio3.00_repeat${repeat}.txt`

* **Subsampled tuning GWAS**:
  `${output_path}/subsample/clean/${trait}_prune_snplist_1_${pop}_tune_GWAS_approxTRUE_ratio3.00_repeat${repeat}.txt`

Example format of these generated subsampled GWAS files:  
* **Subsampled training GWAS example**:
  ```
  SNP         CHR     BP      A1  A2  A1_Frq      BETA           SE           Z         P         N
  rs3131969   1       754182  A   G   2.669e-01  -3.482222e-03  7.831456e-03 -4.446456e-01  6.565759e-01  74586
  rs1048488   1       760912  C   T   1.161e-01  -3.537603e-03  7.797450e-03 -4.536872e-01  6.500540e-01  85875
  rs1806509   1       853954  C   A   4.137e-01   4.359470e-03  7.256207e-03  6.007918e-01  5.479786e-01  60434
  ...
  ```

* **Subsampled tuning GWAS example**:
  ```
  SNP         CHR     BP      A1  A2  A1_Frq      BETA           SE           Z         P         N
  rs3131969   1       754182  A   G   2.669e-01   2.987062e-02  1.356448e-02  2.202121e+00  2.765678e-02  24862
  rs1048488   1       760912  C   T   1.161e-01   1.037583e-02  1.350558e-02  7.682627e-01  4.423312e-01  28625
  rs1806509   1       853954  C   A   4.137e-01   1.163264e-03  1.256812e-02  9.255673e-02  9.262557e-01  20144
  ...
  ```

#### Step2: Estimation of MIXPRS combination weights for different PRS methods
This step includes two substeps:

##### Step 2.1: Obtain LD-pruned PRS beta files
* Use the subsampled training GWAS for the target population obtained from **Step 1**.
* For GWAS from other populations, apply the LD-pruned SNP list from the **target population** (instead of their original populations) to filter and obtain LD-pruned GWAS summary statistics from their original GWAS summary statistics.
* Format these aligned, pruned GWAS summary statistics from all populations for each method and implement each method according to their respective repositories:
  * [JointPRS-auto](https://github.com/LeqiXu/JointPRS)
  * [SDPRX](https://github.com/eldronzhou/SDPRX)

Repeat this procedure for each of the four subsampled training GWAS sets generated in **Step 1** (`repeat=1,2,3,4`). After completing this step, you should obtain beta files similar to the following for each repeat (`repeat=1,2,3,4`):
* `${JointPRS_EUR_beta_file}`, `${JointPRS_EAS_beta_file}`, `${JointPRS_AFR_beta_file}`, `${JointPRS_SAS_beta_file}`, `${JointPRS_AMR_beta_file}`
* `${SDPRX_EUR_beta_file}`, `${SDPRX_EAS_beta_file}`, `${SDPRX_AFR_beta_file}`, `${SDPRX_SAS_beta_file}`, `${SDPRX_AMR_beta_file}`

**Note**: The number of population-specific beta files obtained depends on the availability of GWAS summary statistics for each population.

##### Step 2.2: Obtain PRS combining weights
* Use the subsampled tuning GWAS and format them according to the required MIXPRS summary statistics format detailed in **Section 4**.
  We assume the formatted subsampled tuning GWAS file is named as follows:
  `${output_path}/subsample/MIXPRS/${trait}_prune_snplist_1_${pop}_tune_MIXPRS_approxTRUE_ratio3.00_repeat${repeat}.txt`
* Use the LD-pruned PRS (beta files) obtained from **Step 2.1**.
* Default parameters for this step are:
  * `indep_approx=TRUE`: uses an identity covariance matrix for the pruned SNPs instead of LD reference panels (the `--ref_dir` flag is still required).
  * `selection_criterion=NNLS`: employs the Lawson–Hanson Non-Negative Least Squares (NNLS) algorithm to estimate non-negative PRS combining weights.

Run the following command to calculate the optimal PRS combining weights:

```bash
conda activate MIXPRS

python ${MIXPRS_path}/MIX_linear_weight.py \
  --ref_dir=${ref_data_path}/1KG \
  --sst_file=${output_path}/subsample/MIXPRS/${trait}_prune_snplist_1_${pop}_tune_MIXPRS_approxTRUE_ratio3.00_repeat${repeat}.txt \
  --pop=${pop} \
  --prs_beta_file=${JointPRS_EUR_beta_file},${JointPRS_EAS_beta_file},${JointPRS_AFR_beta_file},${JointPRS_SAS_beta_file},${JointPRS_AMR_beta_file},${SDPRX_EUR_beta_file},${SDPRX_EAS_beta_file},${SDPRX_AFR_beta_file},${SDPRX_SAS_beta_file},${SDPRX_AMR_beta_file} \
  --indep_approx=TRUE \
  --selection_criterion=NNLS \
  --out_dir=${output_path}/Final_weight/no_val \
  --out_name=${trait}_prune_snplist_1_JointPRS_SDPRX_EUR_EAS_AFR_SAS_AMR_repeat${repeat}
```

Repeat this procedure for each of the four subsampled tuning GWAS sets generated in **Step 1** (`repeat=1,2,3,4`). After completing this step, you should obtain four PRS combining weights files, named as follows for each repeat (`repeat=1,2,3,4`):
* **PRS combining weights**:
  `${output_path}/Final_weight/no_val/${trait}_prune_snplist_1_JointPRS_SDPRX_EUR_EAS_AFR_SAS_AMR_repeat${repeat}_${pop}_non_negative_linear_weights_approxTRUE.txt`

**Note**: For each repeat (`repeat=1,2,3,4`), ensure you use the corresponding subsampled tuning GWAS dataset from **Step 1** along with its matching LD-pruned PRS beta file obtained from **Step 2.1**.

#### Step3: Derivation of MIXPRS
This step includes two substeps:

##### Step 3.1: Obtain full SNPs PRS beta files
* Use the original GWAS summary statistics from all populations.
* Format these original GWAS summary statistics according to the requirements of each method and implement each method following their respective repositories:
  * [JointPRS-auto](https://github.com/LeqiXu/JointPRS)
  * [SDPRX](https://github.com/eldronzhou/SDPRX)

After completing this step, you should obtain full SNPs beta files for each method and population, such as:

* `${JointPRS_EUR_beta_file}`, `${JointPRS_EAS_beta_file}`, `${JointPRS_AFR_beta_file}`, `${JointPRS_SAS_beta_file}`, `${JointPRS_AMR_beta_file}`
* `${SDPRX_EUR_beta_file}`, `${SDPRX_EAS_beta_file}`, `${SDPRX_AFR_beta_file}`, `${SDPRX_SAS_beta_file}`, `${SDPRX_AMR_beta_file}`

**Note**: The number of population-specific beta files depends on the availability of GWAS summary statistics for each population.

##### Step 3.2: Obtain MIXPRS
* Use the original GWAS summary statistics for the target population formatted according to **Section 4: Summary Statistics Preparation**. This should match the file used in Step 1:
  ```
  ${summary_stat_path}/${trait}_${pop}_MIXPRS_sumstat.txt
  ```
* Use the calculated PRS combining weights obtained in **Step 2.2**, named as follows:
  ```bash
  weight_file1="${output_path}/Final_weight/no_val/${trait}_prune_snplist_1_JointPRS_SDPRX_EUR_EAS_AFR_SAS_AMR_repeat1_${pop}_non_negative_linear_weights_approxTRUE.txt"
  weight_file2="${output_path}/Final_weight/no_val/${trait}_prune_snplist_1_JointPRS_SDPRX_EUR_EAS_AFR_SAS_AMR_repeat2_${pop}_non_negative_linear_weights_approxTRUE.txt"
  weight_file3="${output_path}/Final_weight/no_val/${trait}_prune_snplist_1_JointPRS_SDPRX_EUR_EAS_AFR_SAS_AMR_repeat3_${pop}_non_negative_linear_weights_approxTRUE.txt"
  weight_file4="${output_path}/Final_weight/no_val/${trait}_prune_snplist_1_JointPRS_SDPRX_EUR_EAS_AFR_SAS_AMR_repeat4_${pop}_non_negative_linear_weights_approxTRUE.txt"
  ```
* Use the full SNPs PRS beta files obtained from **Step 3.1**.

Run the following command to obtain the final MIXPRS:

```bash
conda activate MIXPRS

python ${MIXPRS_path}/MIX_final_combine.py \
  --ref_dir=${ref_data_path}/1KG \
  --sst_file=${summary_stat_path}/${trait}_${pop}_MIXPRS_sumstat.txt \
  --pop=${pop} \
  --prs_beta_file=${JointPRS_EUR_beta_file},${JointPRS_EAS_beta_file},${JointPRS_AFR_beta_file},${JointPRS_SAS_beta_file},${JointPRS_AMR_beta_file},${SDPRX_EUR_beta_file},${SDPRX_EAS_beta_file},${SDPRX_AFR_beta_file},${SDPRX_SAS_beta_file},${SDPRX_AMR_beta_file} \
  --weight_file=${weight_file1},${weight_file2},${weight_file3},${weight_file4} \
  --out_dir=${output_path}/MIXPRS \
  --out_name=${trait}
```

This command produces two results files in the specified output directory (`${output_path}/MIXPRS`):
* **Final MIXPRS beta**:
  ```
  ${output_path}/MIXPRS/${trait}_${pop}_MIXPRS.txt
  ```
* **Weighted PRS beta for each method and population**:
  ```
  ${output_path}/MIXPRS/${trait}_${pop}_MIXPRS_separate.txt
  ```

**Note**: Ensure the order and naming of `--prs_beta_file` are consistent between Step 2.2 and Step 3.2.

## Acknowledgment
Part of the code is adapted from [PRS-CSx](https://github.com/getian107/PRScsx/tree/master). We thank Dr. Tian Ge for sharing his code and LD reference panels.

## Support
Please direct any problems or questions to Leqi Xu (leqi.xu@yale.edu).
