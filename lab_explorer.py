import pandas as pd
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

class LabExplorer:
    def __init__(self, root):
        self.root = root
        self.root.title("Omar Lab: Proteomics Explorer")
        self.root.geometry("800x500")
        
        # Split the window into Left and Right panes
        self.paned = ttk.PanedWindow(root, orient=tk.HORIZONTAL)
        self.paned.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)
        
        # ============================================
        # LEFT SIDE: The File Loader & Gene List
        # ============================================
        self.frame_left = ttk.Frame(self.paned)
        self.paned.add(self.frame_left, weight=1)
        
        # The button to pick your CSV
        self.btn_load = ttk.Button(self.frame_left, text="Load Results CSV", command=self.load_data)
        self.btn_load.pack(fill=tk.X, pady=(0,5))
        
        #The interactive list
        self.tree = ttk.Treeview(self.frame_left, columns=("Gene"), show='headings')
        self.tree.heading("Gene", text="Significant Genes")
        self.tree.pack(fill=tk.BOTH, expand=True)
        self.tree.bind("<<TreeviewSelect>>", self.show_details)
        
        # ============================================
        # RIGHT SIDE: The Details Panel 
        # ============================================
        self.frame_right = ttk.Frame(self.paned)
        self.paned.add(self.frame_right, weight=2)
        
        self.lbl_title = ttk.Label(self.frame_right, text="Protein Details", font=("Helvetica", 14, "bold"))
        self.lbl_title.pack(anchor=tk.W)
        
        # The text box where all the "fluff" goes
        self.details = tk.Text(self.frame_right, wrap=tk.WORD, state=tk.DISABLED, font=("Consolas", 10))
        self.details.pack(fill=tk.BOTH, expand=True, pady=5)
        
        self.df = None
        
    def load_data(self):
        #Opens the standard Windows file picker
        filepath = filedialog.askopenfilename(filetypes=[("CSV Files", "*.csv")])
        if not filepath: return
        
        try: 
            self.df = pd.read_csv(filepath)
            self.tree.delete(*self.tree.get_children()) #Clear old data
            
            # Check for the correct column names
            gene_col = 'PG.Genes' if 'PG.Genes' in self.df.columns else self.df.columns[0]
            
            #Filter to show ONLY significant hits in the list
            if 'Significance_Tier' in self.df.columns:
                sig_df = self.df[self.df['Significance_Tier'].str.contains('High|Exploratory', na=False)]
                genes = sig_df[gene_col].dropna().unique()
            else:
                genes = self.df[gene_col].dropna().unique()
                
            # Add genes to the left-side list
            for gene in genes:
                self.tree.insert("", tk.END, values=(gene,))
                
            messagebox.showinfo("Success", f"Loaded {len(genes)} significant hits!")
        except Exception as e:
            messagebox.showerror("Error", f"Could not load file:\n{e}")
            
    def show_details(self, event):
        #Triggered when you click a gene
        if not self.tree.selection(): return
        selected_item = self.tree.selection()[0]
        gene_name = self.tree.item(selected_item)['values'][0]
        
        gene_col = 'PG.Genes' if 'PG.Genes' in self.df.columns else self.df.columns[0]
        
        gene_data = self.df[self.df[gene_col] == gene_name].iloc[0]
        
        #Unlock the text box, clear it, and write the new data
        self.details.config(state=tk.NORMAL)
        self.details.delete(1.0, tk.END)
        
        #Dynamically loop through Every column and display it
        for col in self.df.columns:
            val = gene_data[col]
            self.details.insert(tk.END, f"{col}:\n", "header")
            self.details.insert(tk.END, f"  {val}\n\n")
            
        self.details.tag_config("header", font=("Consolas", 10, "bold"))
        self.details.config(state=tk.DISABLED) #Lock it so user can't type in it
        
if __name__ == "__main__":
    root = tk.Tk()
    app = LabExplorer(root)
    root.mainloop()