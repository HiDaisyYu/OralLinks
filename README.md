# Oral vs Systemic Disease Cross-Trait Analysis

This repository documents a reproducible pipeline for exploring the shared genetic architecture between **oral/dental diseases** and **systemic diseases** using GWAS datasets.  
The framework is adapted from recent large-scale cross-trait studies (e.g., lung vs gastrointestinal diseases) and extended to oral-systemic interactions.

---

## Methods Framework

We apply a five-stage cross-trait genetic analysis workflow:

### 1. Data Collection
- Integrate GWAS datasets across multiple populations (European, East Asian, African).
- Oral/dental traits: periodontitis, dental caries, tooth loss, oral cancer, gingivitis, etc.
- Systemic traits: cardiovascular disease, metabolic disease, autoimmune disease, cognitive disorders, digestive disease, addiction, head & neck cancers.

### 2. Genetic Correlation Analysis
- **Global level:** Estimate genome-wide genetic correlations using LDSC (Linkage Disequilibrium Score Regression).
- **Local level:** Identify shared chromosomal regions using SUPERGNOVA.
- **Functional regions:** Examine overlap in regulatory elements (e.g., DNase I hypersensitivity sites, transcribed regions).

### 3. Identification of Pleiotropic Variants
- **Cross-trait meta-analysis:** Detect SNPs influencing both oral and systemic traits.
- **Fine-mapping:** Construct 99% credible sets to prioritize causal variants.
- **Colocalization analysis:** Confirm shared signals across traits.
- **Gene-level analysis:** Identify pleiotropic genes associated with multiple disease categories.

### 4. Functional Annotation
- **Gene Set Enrichment Analysis (GSEA):** Highlight immune and inflammatory pathways.
- **Tissue-Specific Enrichment (TSEA):** Assess enrichment in oral tissues (gingiva, oral mucosa) and systemic tissues (heart, colon, pancreas).
- **Cell-Type Specific Enrichment (CSEA):** Identify immune cell involvement (e.g., T cells, neutrophils).
- **TWAS/PWAS:** Integrate transcriptome and proteome data to link gene/protein expression with disease risk.

### 5. Causal Inference and Microbiome Analysis
- **Bi-directional Mendelian Randomization (MR):** Test causal relationships between oral and systemic diseases.
- **Mediation analysis:** Introduce the oral microbiome as a mediator of disease associations.
  - Example pathways:
    - *Porphyromonas gingivalis → systemic inflammation → cardiovascular disease*
    - *Oral microbiome dysbiosis → diabetes → periodontal disease*

---

## Keyword Dictionary

We define keyword groups to classify GWAS traits into oral/dental vs systemic disease categories.

```python
KEYWORDS = {
    'Oral': [
        'periodont', 'gingiv', 'gum', 'tooth', 'teeth', 'dental', 'caries', 'decay',
        'mouth', 'oral', 'tongue', 'saliv', 'xerostomia', 'denture', 'edentul',
        'leukoplakia', 'canker', 'stomatitis'
    ],
    'Cognitive': [
        'alzheimer', 'dementia', 'parkinson', 'cognitive', 'memory', 'schizophrenia',
        'depress', 'anxiety', 'bipolar', 'insomnia', 'neuroticism'
    ],
    'Addiction': [
        'smoking', 'tobacco', 'cigarette', 'nicotine', 'alcohol', 'drink', 'liquor',
        'cannabis', 'marijuana', 'weed', 'dependence'
    ],
    'Autoimmune': [
        'sjogren', 'lupus', 'sle', 'rheumatoid', 'arthritis', 'psoriasis',
        'autoimmune', 'type 1 diabet', 'multiple sclerosis', 'celiac'
    ],
    'Cardio': [
        'coronary', 'myocardial', 'stroke', 'heart failure', 'atrial fib',
        'hypertension', 'blood pressure'
    ],
    'Metabolic': [
        'type 2 diabet', 'bmi', 'obesity',
        'lipid', 'cholesterol', 'hdl', 'ldl', 'triglyceride'
    ],
    'Digestive': [
        'bowel', 'colitis', 'crohn', 'ibd', 'gastric', 'stomach', 'esophag',
        'reflux', 'gerd', 'peptic', 'pancrea', 'liver', 'gut'
    ],
    'Head_Neck_Cancer': [
        'head and neck', 'oral cancer', 'mouth cancer', 'tongue cancer',
        'esophageal cancer', 'nasopharyngeal', 'larynx', 'pharynx'
    ]
}
