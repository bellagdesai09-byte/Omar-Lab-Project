import pandas as pd
xl = pd.ExcelFile('MSdataprotein3-4-26.xlsx')
print("Available sheets:", xl.sheet_names)