#!/usr/bin/env python3
import pandas as pd

# === 配置文件路径 ===
INPUT_FILE = "summary_statistics_table_export.tsv"  # 上传的 GWAS Catalog 导出文件, this file was exported from: https://www.ebi.ac.uk/gwas/downloads/summary-statistics #List of published studies with summary statistics
OUTPUT_FILE = "gwas_targets.tsv"  # 输出筛选结果

# === 关键词字典 (拖网式搜索) === Keywords below were selected from this file "summary_statistics_table_export.tsv"
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

def run_filtering():
    print(f"Loading {INPUT_FILE}...")
    try:
        df = pd.read_csv(INPUT_FILE, sep='\t', low_memory=False)
    except Exception as e:
        print(f"Error reading file: {e}")
        return

    print(f"Total studies in file: {len(df)}")

    tasks = []

    print("Filtering traits...")
    for idx, row in df.iterrows():
        text_content = str(row.get('reportedTrait', '')).lower() + " " + str(row.get('efoTraits', '')).lower()

        matched_cat = None
        for cat, roots in KEYWORDS.items():
            if any(root in text_content for root in roots):
                matched_cat = cat
                break

        if matched_cat:
            tasks.append({
                'Category': matched_cat,
                'Trait': row.get('reportedTrait'),
                'Accession': row.get('accessionId'),
                'Sample_Size': row.get('initialSampleDescription', 'NA')
            })

    if tasks:
        result_df = pd.DataFrame(tasks)
        result_df = result_df.drop_duplicates(subset=['Accession'])
        result_df.to_csv(OUTPUT_FILE, sep='\t', index=False)
        print(f"\n✅ 筛选完成！")
        print(f"共找到 {len(result_df)} 个相关研究。")
        print("分布如下:")
        print(result_df.groupby('Category').size())
        print(f"任务列表已保存至: {OUTPUT_FILE}")
    else:
        print("❌ 未找到匹配的研究。")

if __name__ == "__main__":
    run_filtering()
