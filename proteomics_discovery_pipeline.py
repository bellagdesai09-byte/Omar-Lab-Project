import pandas as pd
import numpy as np
from scipy import stats
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import seaborn as sns
import os

class ProteomicsDiscoveryPipeline:
    def __init__(self, file_path, config):
        print(f"Loading data from {file_path}...")
        self.df = pd.read_csv(file_path)
        self.config = config
        self.results = {}

    def clean_and_format(self):
        print("Cleaning data (removing Decoys)...")
        # 1. Remove ONLY the fake Decoy sequences. 
        decoy_mask = self.df['PG.ProteinGroups'].str.contains('Decoy', na=False, case=False)
        self.df = self.df[~decoy_mask].copy()

        # 2. Map config substrings to actual column names
        self.group_cols = {}
        for logical_name, substrings in self.config['group_mapping'].items():
            matched_cols = []
            for sub in substrings:
                match = [col for col in self.df.columns if sub in col and '.PG.Quantity' in col]
                if match:
                    matched_cols.append(match[0])
            self.group_cols[logical_name] = matched_cols
            
            # Force numeric types
            for col in matched_cols:
                self.df[col] = pd.to_numeric(self.df[col], errors='coerce')

    def impute_and_normalize(self):
        print("Normalizing and Imputing missing values to 0.1...")
        all_quant_cols = [col for cols in self.group_cols.values() for col in cols]
        
        # 1. Median Centering (Calculate medians ignoring NaNs to get accurate scaling)
        sample_medians = self.df[all_quant_cols].median()
        global_median = sample_medians.median()
        
        for col in all_quant_cols:
            self.df[col] = self.df[col] * (global_median / sample_medians[col])

        # 2. Aggressive Imputation (Rescue False Negatives)
        # Replace all remaining NaNs with 0.1 so no data is dropped
        self.df[all_quant_cols] = self.df[all_quant_cols].fillna(0.1)

        # 3. Log2 Transformation
        self.df[all_quant_cols] = np.log2(self.df[all_quant_cols] + 1)

    def run_statistics(self):
        for test_group, control_group in self.config['comparisons']:
            print(f"Running statistics for {test_group} vs {control_group}...")
            
            cols_test = self.group_cols[test_group]
            cols_ctrl = self.group_cols[control_group]

            comp_df = self.df.copy()
            p_values_raw = []
            log2_fcs = []

            # Perform Welch's T-Test
            for index, row in comp_df.iterrows():
                val_test = row[cols_test].values.astype(float)
                val_ctrl = row[cols_ctrl].values.astype(float)

                # Because of imputation, there are no NaNs, so variance won't crash
                # Add a tiny try/except block just in case standard deviation is perfectly 0
                try:
                    t_stat, p_val = stats.ttest_ind(val_test, val_ctrl, equal_var=False)
                    lfc = np.mean(val_test) - np.mean(val_ctrl)
                except:
                    p_val = np.nan
                    lfc = np.nan
                
                p_values_raw.append(p_val)
                log2_fcs.append(lfc)

            comp_df['Log2FC'] = log2_fcs
            comp_df['Raw_P_Value'] = p_values_raw

            # Drop math errors (if any) before FDR
            comp_df = comp_df.dropna(subset=['Raw_P_Value', 'Log2FC'])

            # Apply FDR correction
            reject, pvals_corrected, _, _ = multipletests(comp_df['Raw_P_Value'], alpha=self.config['p_value_threshold'], method='fdr_bh')
            comp_df['FDR_Adjusted_P'] = pvals_corrected
            
            # Create Tiered Flags
            pass_raw = (comp_df['Raw_P_Value'] < self.config['p_value_threshold']) & (comp_df['Log2FC'].abs() > self.config['log2fc_threshold'])
            pass_fdr = (comp_df['FDR_Adjusted_P'] < self.config['p_value_threshold']) & (comp_df['Log2FC'].abs() > self.config['log2fc_threshold'])
            
            comp_df['Passed_Raw_Threshold'] = pass_raw
            comp_df['Passed_FDR_Threshold'] = pass_fdr

            # Assign categories for the plot
            conditions = [pass_fdr, pass_raw & ~pass_fdr]
            choices = ['High Confidence (FDR < 0.05)', 'Exploratory (Raw P < 0.05)']
            comp_df['Significance_Tier'] = np.select(conditions, choices, default='Not Significant')

            self.results[f"{test_group}_vs_{control_group}"] = comp_df
            print(f"  -> {pass_fdr.sum()} High Confidence Hits | {(pass_raw & ~pass_fdr).sum()} Exploratory Hits")

    def export_and_plot(self):
        if not os.path.exists("Output"):
            os.makedirs("Output")

        for comp_name, df_res in self.results.items():
            # 1. Export CSV
            output_file = f"Output/{comp_name}_Results.csv"
            export_cols = ['PG.ProteinGroups', 'PG.Genes', 'Log2FC', 'Raw_P_Value', 'FDR_Adjusted_P', 'Passed_Raw_Threshold', 'Passed_FDR_Threshold', 'Significance_Tier']
            df_res[export_cols].to_csv(output_file, index=False)
            
            # 2. Multi-Tier Volcano Plot
            plt.figure(figsize=(10, 6))
            df_res['neg_log_raw_p'] = -np.log10(df_res['Raw_P_Value'])
            
            palette = {
                'High Confidence (FDR < 0.05)': 'red',
                'Exploratory (Raw P < 0.05)': 'royalblue',
                'Not Significant': 'lightgrey'
            }

            sns.scatterplot(
                data=df_res, x='Log2FC', y='neg_log_raw_p',
                hue='Significance_Tier', palette=palette,
                alpha=0.7, edgecolor=None
            )

            # Draw threshold lines
            plt.axhline(-np.log10(self.config['p_value_threshold']), color='black', linestyle='--', linewidth=1, label='Raw p = 0.05')
            plt.axvline(self.config['log2fc_threshold'], color='black', linestyle=':', linewidth=1)
            plt.axvline(-self.config['log2fc_threshold'], color='black', linestyle=':', linewidth=1)

            plt.title(f'Discovery Volcano Plot: {comp_name}')
            plt.xlabel('Log2 Fold Change')
            plt.ylabel('-Log10(Raw P-Value)')
            
            # Move legend outside the plot
            plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            
            plot_file = f"Output/{comp_name}_Volcano.png"
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()


# ==========================================
# CONFIGURATION BLOCK (Abstracted for any experiment)
# ==========================================
if __name__ == "__main__":
    
    experiment_config = {
        # Notice U2oS is completely ignored here. Only the mapped groups are processed.
        # You can name the keys whatever you want (e.g., 'Treatment_1', 'Baseline').
        'group_mapping': {
            'Treatment_C': ['HEX_C2', 'HEX_C3', 'HEX_C4', 'HEX_C5'],
            'Treatment_N': ['HEX_N2', 'HEX_N3', 'HEX_N4', 'HEX_N5'],
            'Baseline_Control':  ['HEX_NES2', 'HEX_NES3', 'HEX_NES4', 'HEX_NES5']
        },
        
        # Define the exact comparisons you want. 
        'comparisons': [
            ('Treatment_C', 'Baseline_Control'),  # C vs Control
            ('Treatment_N', 'Baseline_Control'),  # N vs Control
            ('Treatment_C', 'Treatment_N')        # C vs N
        ],
        
        # Statistical Cutoffs
        'p_value_threshold': 0.05,    
        'log2fc_threshold': 1.0       
    }

    # Run the Pipeline
    pipeline = ProteomicsDiscoveryPipeline('BioLab_Report_for_Excel.csv', experiment_config)
    pipeline.clean_and_format()
    pipeline.impute_and_normalize()
    pipeline.run_statistics()
    pipeline.export_and_plot()
    
    print("Analysis Complete! Check the 'Output' folder.")