import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import os
import warnings
from adjustText import adjust_text
import plotly.express as px

class ProteomicsDiscoveryPipeline:
    def __init__(self, file_path, config):
        print(f"Loading data from {file_path}...")
        
        # Smart Separator Detector: Handles both your tab-files and colleague's comma-files
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
            first_line = f.readline()
            file_sep = '\t' if '\t' in first_line else ','
            
        self.df = pd.read_csv(file_path, sep=file_sep)
        self.config = config
        self.results = {}

    def clean_and_format(self):
        print("Cleaning data (removing Decoys and background)...")
        if 'PG.ProteinGroups' in self.df.columns:
            decoy_mask = self.df['PG.ProteinGroups'].str.contains('Decoy', na=False, case=False)
            self.df = self.df[~decoy_mask].copy()

        # Wide-Format U2OS Assassin: Purges background noise
        u2os_cols = [col for col in self.df.columns if 'U2OS' in col.upper()]
        self.df = self.df.drop(columns=u2os_cols)
        
        if len(u2os_cols) > 0:
            print(f"  [-] Successfully purged {len(u2os_cols)} U2OS columns from the dataset.")

        self.group_cols = {}
        for logical_name, substrings in self.config['group_mapping'].items():
            matched_cols = []
            for sub in substrings:
                # Case-insensitive search for all matching replicates
                match = [col for col in self.df.columns if sub.lower() in col.lower() and '.pg.quantity' in col.lower()]
                if match:
                    matched_cols.extend(match)
            
            self.group_cols[logical_name] = list(set(matched_cols))
            print(f"  [+] Mapped '{logical_name}' to {len(self.group_cols[logical_name])} columns.")
            
            for col in self.group_cols[logical_name]:
                self.df[col] = pd.to_numeric(self.df[col], errors='coerce')

    def impute_and_normalize(self):
        print("\n" + "="*50)
        print(" PHASE 1: DATA PRE-PROCESSING")
        print("="*50)
        
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
            
            cols_test = self.group_cols[test_group]
            cols_ctrl = self.group_cols[control_group]

            if len(cols_test) < 2 or len(cols_ctrl) < 2:
                print(f"  [!] ERROR: Not enough replicates to run T-Test for {comp_label}!")
                continue

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
            
            # --- UPDATED 3-TIER LOGIC ---
            pass_fdr = (comp_df['FDR_Adjusted_P'] < self.config['p_value_threshold']) & (comp_df['Log2FC'].abs() > self.config['log2fc_threshold'])
            pass_raw_fc = (comp_df['Raw_P_Value'] < self.config['p_value_threshold']) & (comp_df['Log2FC'].abs() > self.config['log2fc_threshold'])
            pass_raw_only = (comp_df['Raw_P_Value'] < self.config['p_value_threshold']) & (comp_df['Log2FC'].abs() <= self.config['log2fc_threshold'])

            conditions = [pass_fdr, pass_raw_fc & ~pass_fdr, pass_raw_only]
            choices = ['Tier 1: High Confidence (FDR)', 'Tier 2: Exploratory (Raw P + FC)', 'Tier 3: Statistical Trend (Raw P)']
            comp_df['Significance_Tier'] = np.select(conditions, choices, default='Not Significant')

            self.results[f"{test_group}_vs_{control_group}"] = comp_df
            
            print(f"  Tier 1 (FDR):        {pass_fdr.sum():>4}")
            print(f"  Tier 2 (P + FC):     {(pass_raw_fc & ~pass_fdr).sum():>4}")
            print(f"  Tier 3 (P only):     {pass_raw_only.sum():>4}")

    def export_and_plot(self):
        print("\n" + "="*50)
        print(" PHASE 3: EXPORTING VISUALS (PNG & HTML)")
        print("="*50)
        if not os.path.exists("Output"):
            os.makedirs("Output")

        for comp_name, df_res in self.results.items():
            print(f"  [+] Processing: {comp_name}...")
            
            output_file = f"Output/{comp_name}_Results.csv"
            
            # These are the columns we want at the very front of the CSV for the GUI
            priority_cols = [
                'PG.ProteinGroups', 
                'PG.Genes', 
                'PG.ProteinDescriptions', 
                'PG.Coverage', 
                'PG.Qvalue', 
                'Log2FC', 
                'Raw_P_Value', 
                'FDR_Adjusted_P', 
                'Significance_Tier'
            ]
            
            # Organize columns: priority first, then everything else
            stat_cols = [c for c in priority_cols if c in df_res.columns]
            other_cols = [c for c in df_res.columns if c not in stat_cols and c != 'neg_log_raw_p']
            
            # Save the final Results CSV
            df_res[stat_cols + other_cols].to_csv(output_file, index=False)
            
            # --- Visualization Math & Plots ---
            df_res['neg_log_raw_p'] = -np.log10(df_res['Raw_P_Value'])
            palette = {
                'Tier 1: High Confidence (FDR)': 'red', 
                'Tier 2: Exploratory (Raw P + FC)': 'royalblue', 
                'Tier 3: Statistical Trend (Raw P)': 'green',
                'Not Significant': 'lightgrey'
            }

            # 1. Create Static PNG Volcano Plot
            plt.figure(figsize=(10, 6))
            sns.scatterplot(data=df_res, x='Log2FC', y='neg_log_raw_p', hue='Significance_Tier', palette=palette, alpha=0.7, edgecolor=None)
            plt.axhline(-np.log10(self.config['p_value_threshold']), color='black', linestyle='--', linewidth=1)
            plt.axvline(self.config['log2fc_threshold'], color='black', linestyle=':', linewidth=1)
            plt.axvline(-self.config['log2fc_threshold'], color='black', linestyle=':', linewidth=1)
            
            texts = []
            top_hits = df_res[df_res['Significance_Tier'] != 'Not Significant'].nsmallest(15, 'Raw_P_Value')
            for i, row in top_hits.iterrows():
                texts.append(plt.text(row['Log2FC'], row['neg_log_raw_p'], str(row['PG.Genes']), fontsize=8))
            if texts:
                adjust_text(texts, arrowprops=dict(arrowstyle='-', color='black', lw=0.5))

            plt.title(f'Discovery Volcano Plot: {comp_name}')
            plt.savefig(f"Output/{comp_name}_Volcano.png", dpi=300, bbox_inches='tight')
            plt.close()

            # 2. Create Interactive HTML Volcano Plot
            hover_dict = {'Log2FC': ':.2f', 'Raw_P_Value': ':.4f'}
            if 'PG.ProteinDescriptions' in df_res.columns:
                hover_dict['PG.ProteinDescriptions'] = True

            fig = px.scatter(
                df_res, x='Log2FC', y='neg_log_raw_p', color='Significance_Tier',
                hover_name='PG.Genes', hover_data=hover_dict,
                color_discrete_map=palette,
                title=f'Interactive Volcano: {comp_name}'
            )
            fig.write_html(f"Output/{comp_name}_Interactive_Volcano.html")
            
        print("[OK] Exports saved to 'Output' folder with full metadata.\n")

if __name__ == "__main__":
#     # --- CONFIGURATION SWITCHBOARD ---
#     # Toggle between your project and the colleague's data by commenting/uncommenting
    
#     # --- OPTION A: Colleague's Alpha/Beta Data ---
#     experiment_config = {
#         'group_mapping': {
#             'Alpha': ['alpha_'], 
#             'Beta':  ['beta_'],
#             'QC':    ['QC_']    
#         },
#         'comparisons': [('Alpha', 'Beta')],
#         'p_value_threshold': 0.05,    
#         'log2fc_threshold': 1.0       
#     }
#     data_file = 'Translated_Alpha_Beta_Data.csv'

# --- OPTION B: Original AKAP11 Data ---
    experiment_config = {
        'group_mapping': {
            'Treatment_C': ['HEX_C2', 'HEX_C3', 'HEX_C4', 'HEX_C5'],
            'Treatment_N': ['HEX_N2', 'HEX_N3', 'HEX_N4', 'HEX_N5'],
            'Baseline_Control': ['HEX_NES2', 'HEX_NES3', 'HEX_NES4', 'HEX_NES5']
        },
        'comparisons': [('Treatment_C', 'Baseline_Control'), ('Treatment_N', 'Baseline_Control')],
        'p_value_threshold': 0.05,    
        'log2fc_threshold': 1.0       
    }
    
    data_file = '20260313_093105_Omar_AKAP11_ServiceAward_DIA_TT_20250701_Report_pivot for Bella (1).csv'
    
    # These must be aligned with the experiment_config above!
    pipeline = ProteomicsDiscoveryPipeline(data_file, experiment_config)
    pipeline.clean_and_format()
    pipeline.impute_and_normalize()
    pipeline.run_statistics()
    pipeline.export_and_plot()