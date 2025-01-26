# **Breast Cancer Prognostic Gene Signature**

## **Project Overview**

This project aims to develop a **breast tumor-based gene signature** that can be used to predict clinical outcomes for breast cancer patients. By leveraging **genome-wide gene expression data**, we analyze multiple patient cohorts to identify significant genes associated with survival outcomes.

## **Objectives**

- Identify a **gene signature with <100 genes** to predict patient survival.
- Develop the signature using **one or more discovery cohorts**.
- Validate the performance of the gene signature in independent **validation cohorts**.
- Generate **Kaplan-Meier survival curves** to visualize the prognostic impact.

---

## **Study Cohorts**

We analyze genome-wide gene expression data obtained from **Affymetrix GeneChip™ Human Genome U133 Array** across six different breast cancer patient cohorts:

| Cohort | Expression Data | Phenotype Data   |
| ------ | --------------- | ---------------- |
| FRA    | `FRA_expr.txt`  | `FRA_pheno.txt`  |
| GER    | `GER_expr.txt`  | `GER_pheno.txt`  |
| SWE    | `SWE_expr.txt`  | `SWE_pheno.txt`  |
| TWN    | `TWN_expr.txt`  | `TWN_pheno.txt`  |
| USA1   | `USA1_expr.txt` | `USA1_pheno.txt` |
| USA2   | `USA2_expr.txt` | `USA2_pheno.txt` |

---

## **Methodology**

### **1. Data Preparation**

- Load gene expression and phenotype data for each cohort.
- Match phenotype data with gene expression matrices.

### **2. Gene Signature Selection**

- Perform **Cox proportional hazards regression** to identify genes significantly associated with patient survival.
- Adjust **p-values** using the **Bonferroni** or **Benjamini-Hochberg** method.
- Retain genes with significant hazard ratios (`HR`).

### **3. Signature Validation**

- Apply the identified gene signature to **validation cohorts**.
- Compute **risk scores** and classify patients into **high-risk** or **low-risk** groups.
- Validate the prognostic performance using:
  - **Cox proportional hazards model**
  - **Kaplan-Meier survival analysis**

---

## **Process Description**

1. **Data Loading:**

   - Read gene expression and phenotype data files for all cohorts.
   - Match sample identifiers in expression and phenotype datasets.

2. **Signature Development:**

   - Perform survival analysis using **Cox proportional hazards regression** on each discovery cohort.
   - Adjust p-values to control for multiple comparisons.
   - Select significant genes with p-values below the adjusted threshold.
   - Create a **gene signature file** containing gene names and associated risk weights.

3. **Signature Validation:**

   - Apply the gene signature to independent validation cohorts.
   - Compute risk scores for each patient based on gene expression levels.
   - Classify patients into **high-risk** or **low-risk** groups.
   - Use **Kaplan-Meier survival analysis** to compare survival probabilities between risk groups.
   - Plot survival curves and perform **log-rank tests** to assess statistical significance.

---

## **Results**

The gene signature was validated across multiple independent datasets, and **Kaplan-Meier survival curves** were generated for each cohort:

### **Validation Cohort: France**



### **Validation Cohort: Germany**



### **Validation Cohort: Sweden**



### **Validation Cohort: Taiwan**



### **Validation Cohort: USA1**



### **Validation Cohort: USA2**



---

## **Key Findings**

- **High-risk** patients (red line) consistently showed lower survival probabilities compared to **low-risk** patients (blue line).
- The gene signature effectively stratified patients into two distinct risk groups across all validation cohorts.
- Significant differences in survival curves were observed using the **Cox proportional hazards model**.

---

## **How to Use This Repository**

### **1. Clone the Repository**

```sh
git clone https://github.com/your-repo/breast-cancer-signature.git
cd breast-cancer-signature
```

### **2. Install Dependencies in R**

```r
install.packages("survival")
install.packages("dplyr")
```

### **3. Run the Analysis**

Execute the **breast\_cancer\_signature.R** script in R:

```r
source("breast_cancer_signature.R")
```

---

## **Future Work**

- Further optimize the gene signature using **machine learning models**.
- Expand validation to **larger multi-omics datasets**.
- Investigate the **biological pathways** associated with selected genes.

---

## **References**

1. Cox, D. R. (1972). Regression models and life tables. *Journal of the Royal Statistical Society*.
2. Affymetrix GeneChip™ Human Genome U133 Array.
3. Kaplan, E. L., & Meier, P. (1958). Nonparametric estimation from incomplete observations.

---

## **License**

This project is licensed under the **MIT License**. Feel free to contribute and modify the code.

---

## **Contributions**

We welcome contributions! Please create an **issue** or submit a **pull request** for improvements.

