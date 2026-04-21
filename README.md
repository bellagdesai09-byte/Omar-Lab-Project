# Omar-Lab-Project
**Version 1.0 **

This tool is designed to automate the statistical analysis of Mass Spec (DIA) data. It handles the "dirty work" of cleaning data, fixing missing values, and generating plots.

---


## Quick Start Guide (No Coding Required)

### 1. Setup the Data
* Save your proteomics report as a **.csv** file.
* Make sure your file has these columns: `PG.ProteinGroups`, `PG.Genes`, and your sample columns (ending in `.PG.Quantity`).
* Place the file inside the `Omar-Lab-Project` folder.

### 2. Configure the Experiment
Open `proteomics_discovery_pipeline.py` and scroll to the bottom. Update the **Experiment Config** block to match your current run:
```python
experiment_config = {
    'group_mapping': {
        'Treatment_Name': ['Sample_ID_1', 'Sample_ID_2'], 
        'Control_Name':   ['Control_ID_1', 'Control_ID_2']
    },
    'p_value_threshold': 0.05,
    'log2fc_threshold': 1.0
}

### 3. Run the Analysis
Open your terminal (WSL/Ubuntu) and type these two lines:

Bash
source venv/bin/activate
python3 proteomics_discovery_pipeline.py

### 4.  Understanding the Output
After the program finishes, look in the new /Output folder:
The Results Table (.csv)

High Confidence: Hits that passed the strict FDR (False Discovery Rate) correction. 

Use these for main figures.Exploratory: Hits that passed the raw p-value threshold but not FDR. 

Use these for hypothesis generation and "rescuing" lower-abundance proteins.

The Volcano Plot (.png)

A visual map of your data:
🔴 Red Dots: High Confidence hits (P < 0.05, Fold Change > 2).
🔵 Blue Dots: Exploratory hits (P < 0.05).
🔘 Grey Dots: Non-significant background.