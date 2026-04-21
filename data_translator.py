import pandas as pd

SHEET_NAME = 'MSdataprotein3-4-26' 
print(f"Loading sheet '{SHEET_NAME}'...")
df = pd.read_excel('MSdataprotein3-4-26.xlsx', sheet_name=SHEET_NAME) 
df.columns = df.columns.str.strip()

# 1. Split out the quantities for the math
df['New_Column_Name'] = df['R.Condition'].astype(str) + "_" + df['R.FileName'].astype(str) + ".PG.Quantity"
df_quant = df.pivot_table(
    index='PG.ProteinNames', 
    columns='New_Column_Name', 
    values='PG.Quantity', 
    aggfunc='mean'
).reset_index()

# 2. Split out the metadata for the GUI (dropping duplicates so we have 1 row per protein)
# Add any columns here you want to see in your side-panel
gui_columns = ['PG.ProteinNames', 'PG.ProteinAccessions', 'PG.ProteinDescriptions', 'PG.Coverage', 'PG.Qvalue']
df_meta = df[gui_columns].drop_duplicates(subset=['PG.ProteinNames'])

# 3. Merge them back together safely
df_final = pd.merge(df_meta, df_quant, on='PG.ProteinNames', how='inner')

# Rename for the pipeline
df_final.rename(columns={'PG.ProteinNames': 'PG.Genes', 'PG.ProteinAccessions': 'PG.ProteinGroups'}, inplace=True)

df_final.to_csv('Translated_Alpha_Beta_Data.csv', index=False)
print(f"[OK] Ready! {len(df_final)} proteins matched with full GUI metadata.")