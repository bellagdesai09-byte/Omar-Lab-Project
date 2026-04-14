import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings

class ProteomicsDiscoveryPipeline:
    def __init__(self, file_path, config):
        print(f"Loading data from {file_path}...")
        # Use sep='\t' because your file is actually Tab-Separated
        self.df = pd.read_csv(file_path, sep='\t')
        self.config = config
        self.results = {}

    def clean_and_format(self):
        print("Cleaning data (removing Decoys)...")
        decoy_mask = self.df['PG.ProteinGroups'].str.contains('Decoy', na=False, case=False)
        self.df = self.df[~decoy_mask].copy()

        self.group_cols = {}
        for logical_name, substrings in self.config['group_mapping'].items():
            matched_cols = []
            for sub in substrings:
                match = [col for col in self.df.columns if sub in col and '.PG.Quantity' in col]
                if match:
                    matched_cols.append(match[0])
            self.group_cols[logical_name] = matched_cols
            
            for col in matched_cols:
                self.df[col] = pd.to_numeric(self.df[col], errors='coerce')

    def impute_and_normalize(self):
        print("\n" + "="*50)
        print(" PHASE 1: DATA PRE-PROCESSING")
        print("="*50)
        print(f"[*] Normalizing and Imputing missing values...")
        
        all_quant_cols = [col for cols in self.group_cols.values() for col in cols]
        sample_medians = self.df[all_quant_cols].median()
        global_median = sample_medians.median()
        
        for col in all_quant_cols:
            self.df[col] = self.df[col] * (global_median / sample_medians[col])

        self.df[all_quant_cols] = self.df[all_quant_cols].fillna(0.1)
        self.df[all_quant_cols] = np.log2(self.df[all_quant_cols] + 1)
        print(f"[OK] Normalization complete. Base-line floor set to 0.1.")

    def run_statistics(self):
        print("\n" + "="*50)
        print(" PHASE 2: STATISTICAL TIERING")
        print("="*50)
        warnings.filterwarnings('ignore', category=RuntimeWarning)

        for test_group, control_group in self.config['comparisons']:
            comp_label = f"{test_group} vs {control_group}"
            print(f"\n▶ ANALYZING: {comp_label}")
            print("-" * 30)
            
            cols_test = self.group_cols[test_group]
            cols_ctrl = self.group_cols[control_group]

            comp_df = self.df.copy()
            p_values_raw, log2_fcs = [], []

            for index, row in comp_df.iterrows():
                val_test = row[cols_test].values.astype(float)
                val_ctrl = row[cols_ctrl].values.astype(float)
                try:
                    t_stat, p_val = stats.ttest_ind(val_test, val_ctrl, equal_var=False)
                    lfc = np.mean(val_test) - np.mean(val_ctrl)
                except:
                    p_val, lfc = np.nan, np.nan
                
                p_values_raw.append(p_val)
                log2_fcs.append(lfc)

            comp_df['Log2FC'] = log2_fcs
            comp_df['Raw_P_Value'] = p_values_raw
            comp_df = comp_df.dropna(subset=['Raw_P_Value', 'Log2FC'])

            _, pvals_corrected, _, _ = multipletests(comp_df['Raw_P_Value'], alpha=self.config['p_value_threshold'], method='fdr_bh')
            comp_df['FDR_Adjusted_P'] = pvals_corrected
            
            pass_raw = (comp_df['Raw_P_Value'] < self.config['p_value_threshold']) & (comp_df['Log2FC'].abs() > self.config['log2fc_threshold'])
            pass_fdr = (comp_df['FDR_Adjusted_P'] < self.config['p_value_threshold']) & (comp_df['Log2FC'].abs() > self.config['log2fc_threshold'])
            
            comp_df['Passed_Raw_Threshold'] = pass_raw
            comp_df['Passed_FDR_Threshold'] = pass_fdr

            conditions = [pass_fdr, pass_raw & ~pass_fdr]
            choices = ['High Confidence (FDR < 0.05)', 'Exploratory (Raw P < 0.05)']
            comp_df['Significance_Tier'] = np.select(conditions, choices, default='Not Significant')

            self.results[f"{test_group}_vs_{control_group}"] = comp_df
            
            print(f"  Tier 1 (High Confidence):   {pass_fdr.sum():>4} proteins")
            print(f"  Tier 2 (Exploratory):       {(pass_raw & ~pass_fdr).sum():>4} proteins")
            print(f"  [√] Comparison Successful")

    def export_and_plot(self):
        print("\n" + "="*50)
        print(" PHASE 3: EXPORTING VISUALS")
        print("="*50)
        if not os.path.exists("Output"):
            os.makedirs("Output")

        for comp_name, df_res in self.results.items():
            output_file = f"Output/{comp_name}_Results.csv"
            export_cols = ['PG.ProteinGroups', 'PG.Genes', 'Log2FC', 'Raw_P_Value', 'FDR_Adjusted_P', 'Passed_Raw_Threshold', 'Passed_FDR_Threshold', 'Significance_Tier']
            df_res[export_cols].to_csv(output_file, index=False)
            
            plt.figure(figsize=(10, 6))
            df_res['neg_log_raw_p'] = -np.log10(df_res['Raw_P_Value'])
            palette = {'High Confidence (FDR < 0.05)': 'red', 'Exploratory (Raw P < 0.05)': 'royalblue', 'Not Significant': 'lightgrey'}

            sns.scatterplot(data=df_res, x='Log2FC', y='neg_log_raw_p', hue='Significance_Tier', palette=palette, alpha=0.7, edgecolor=None)
            plt.axhline(-np.log10(self.config['p_value_threshold']), color='black', linestyle='--', linewidth=1)
            plt.axvline(self.config['log2fc_threshold'], color='black', linestyle=':', linewidth=1)
            plt.axvline(-self.config['log2fc_threshold'], color='black', linestyle=':', linewidth=1)
            plt.title(f'Discovery Volcano Plot: {comp_name}')
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            plt.savefig(f"Output/{comp_name}_Volcano.png", dpi=300, bbox_inches='tight')
            plt.close()
        print("[OK] All plots saved to 'Output' folder.\n")

if __name__ == "__main__":
    experiment_config = {
        'group_mapping': {
            'Treatment_C': ['HEX_C2', 'HEX_C3', 'HEX_C4', 'HEX_C5'],
            'Treatment_N': ['HEX_N2', 'HEX_N3', 'HEX_N4', 'HEX_N5'],
            'Baseline_Control':  ['HEX_NES2', 'HEX_NES3', 'HEX_NES4', 'HEX_NES5']
        },
        'comparisons': [('Treatment_C', 'Baseline_Control'), ('Treatment_N', 'Baseline_Control'), ('Treatment_C', 'Treatment_N')],
        'p_value_threshold': 0.05,    
        'log2fc_threshold': 1.0       
    }

    # UPDATE FILE NAME HERE
    data_file = '20260313_093105_Omar_AKAP11_ServiceAward_DIA_TT_20250701_Report_pivot for Bella (1).csv'
    
    pipeline = ProteomicsDiscoveryPipeline(data_file, experiment_config)
    pipeline.clean_and_format()
    pipeline.impute_and_normalize()
    pipeline.run_statistics()
    pipeline.export_and_plot()