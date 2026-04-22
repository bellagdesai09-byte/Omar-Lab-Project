import customtkinter as ctk
import pandas as pd
import numpy as np
import requests
import threading
from tkinter import filedialog
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg

# --- THEME SETTINGS ---
ctk.set_appearance_mode("Dark")
ctk.set_default_color_theme("blue")

class LabDiscoveryDatabase(ctk.CTk):
    def __init__(self):
        super().__init__()

        self.title("Omar Lab | Protein Discovery Hub v3.1")
        self.geometry("1400x900")

        # --- GRID SYSTEM ---
        self.grid_columnconfigure(1, weight=1)
        self.grid_rowconfigure(0, weight=1)

        # --- SIDEBAR (Left) ---
        self.sidebar = ctk.CTkFrame(self, width=240, corner_radius=0)
        self.sidebar.grid(row=0, column=0, sticky="nsew")
        
        ctk.CTkLabel(self.sidebar, text="DISCOVERY HUB", font=("Inter", 22, "bold")).pack(pady=30)
        
        self.load_btn = ctk.CTkButton(self.sidebar, text="📂 Load Results CSV", command=self.load_file, height=45)
        self.load_btn.pack(padx=20, pady=10)

        ctk.CTkLabel(self.sidebar, text="Filter by Significance:").pack(pady=(20,0))
        self.tier_menu = ctk.CTkOptionMenu(self.sidebar, 
                                          values=["All Tiers", "Tier 1: High Confidence", "Tier 2: Exploratory", "Tier 3: Statistical Trend"],
                                          command=lambda x: self.apply_filter())
        self.tier_menu.pack(padx=10, pady=10)

        ctk.CTkLabel(self.sidebar, text="Switch View:").pack(pady=(20,0))
        self.view_selector = ctk.CTkSegmentedButton(self.sidebar, 
                                                     values=["Table View", "Volcano Plot", "PCA Dashboard"], 
                                                     command=self.toggle_view)
        self.view_selector.set("Table View")
        self.view_selector.pack(padx=10, pady=10)

        self.stats_label = ctk.CTkLabel(self.sidebar, text="", font=("Inter", 12))
        self.stats_label.pack(side="bottom", pady=20)

        # --- MAIN CONTAINER (Center) ---
        self.container = ctk.CTkFrame(self, fg_color="transparent")
        self.container.grid(row=0, column=1, padx=20, pady=20, sticky="nsew")
        self.container.grid_rowconfigure(0, weight=1)
        self.container.grid_columnconfigure(0, weight=1)

        # View Frames
        self.table_view = ctk.CTkFrame(self.container, fg_color="transparent")
        self.table_view.grid(row=0, column=0, sticky="nsew")
        
        self.search_bar = ctk.CTkEntry(self.table_view, placeholder_text="🔍 Search Gene Symbols or Tags...", height=45)
        self.search_bar.pack(fill="x", pady=(0, 15))
        self.search_bar.bind("<KeyRelease>", lambda e: self.apply_filter())

        self.scroll_frame = ctk.CTkScrollableFrame(self.table_view, label_text="Experimental Discoveries")
        self.scroll_frame.pack(fill="both", expand=True)

        self.pca_view = ctk.CTkFrame(self.container, corner_radius=15)
        self.volcano_view = ctk.CTkFrame(self.container, corner_radius=15)
        
        # --- RIGHT INTEL PANEL ---
        self.intel_panel = ctk.CTkFrame(self, width=420, corner_radius=15)
        self.intel_panel.grid(row=0, column=2, padx=20, pady=20, sticky="nsew")
        
        ctk.CTkLabel(self.intel_panel, text="Intelligence Panel", font=("Inter", 18, "bold")).pack(pady=20)

        self.info_box = ctk.CTkTextbox(self.intel_panel, width=380, height=700, font=("Inter", 13), wrap="word")
        self.info_box.pack(padx=15, pady=10)
        self.info_box.configure(state="disabled")

        # Memory Anchors
        self.data = None
        self.buttons = []
        self.pca_samples = []
        self.v_canvas = None
        self.p_canvas = None

    def load_file(self):
        path = filedialog.askopenfilename(filetypes=[("CSV files", "*.csv")])
        if path:
            self.data = pd.read_csv(path)
            self.tag_biology() 
            self.apply_filter()
            self.render_pca_plot()
            self.render_volcano_plot()
            self.stats_label.configure(text=f"Loaded {len(self.data)} Proteins")

    def tag_biology(self):
        def search_tags(desc):
            desc = str(desc).lower()
            tags = []
            if 'kinase' in desc: tags.append("[KINASE]")
            if 'zinc finger' in desc: tags.append("[DNA-BINDING]")
            if 'mitochondr' in desc: tags.append("[MITO]")
            if 'transport' in desc: tags.append("[TRANSPORT]")
            if 'receptor' in desc: tags.append("[RECEPTOR]")
            return " ".join(tags)

        desc_col = next((col for col in self.data.columns if 'description' in col.lower()), None)
        if desc_col:
            self.data['Bio_Tags'] = self.data[desc_col].apply(search_tags)
        else:
            self.data['Bio_Tags'] = ""

    def toggle_view(self, view):
        self.table_view.grid_forget()
        self.pca_view.grid_forget()
        self.volcano_view.grid_forget()
        if view == "Table View": self.table_view.grid(row=0, column=0, sticky="nsew")
        elif view == "PCA Dashboard": self.pca_view.grid(row=0, column=0, sticky="nsew")
        elif view == "Volcano Plot": self.volcano_view.grid(row=0, column=0, sticky="nsew")

    def apply_filter(self):
        if self.data is None: return
        for b in self.buttons: b.destroy()
        self.buttons = []

        q = self.search_bar.get().strip().lower()
        t_choice = self.tier_menu.get()
        df = self.data.copy()
        
        if q:
            df = df[df['PG.Genes'].str.lower().str.contains(q, na=False) | 
                    df['Bio_Tags'].str.lower().str.contains(q, na=False)]
                    
        if "Tier 1" in t_choice: df = df[df['Significance_Tier'].str.contains("Tier 1", na=False)]
        elif "Tier 2" in t_choice: df = df[df['Significance_Tier'].str.contains("Tier 2", na=False)]
        elif "Tier 3" in t_choice: df = df[df['Significance_Tier'].str.contains("Tier 3", na=False)]

        for _, row in df.head(50).iterrows():
            bg = "#2b2b2b"
            if "Tier 1" in str(row['Significance_Tier']): bg = "#5c1a1a"
            elif "Tier 2" in str(row['Significance_Tier']): bg = "#1a315c"

            btn = ctk.CTkButton(self.scroll_frame, text=f"{row['PG.Genes']}  {row['Bio_Tags']} | FC: {row['Log2FC']:.2f}", 
                                fg_color=bg, anchor="w", command=lambda r=row: self.get_external_intel(r))
            btn.pack(fill="x", padx=5, pady=3)
            self.buttons.append(btn)

    def render_pca_plot(self):
        for w in self.pca_view.winfo_children(): w.destroy()
        ctk.CTkLabel(self.pca_view, text="Interactive PCA (Click dots for Sample Info)", font=("Inter", 18, "bold")).pack(pady=20)

        self.pca_samples = [c for c in self.data.columns if '.PG.Quantity' in c and 'QC' not in c]
        if len(self.pca_samples) < 3: return

        pca_engine = PCA(n_components=2)
        coords = pca_engine.fit_transform(self.data[self.pca_samples].T)

        fig, ax = plt.subplots(figsize=(7, 5), facecolor='#2b2b2b')
        ax.set_facecolor('#2b2b2b')
        
        colors = ['#ff4b4b' if 'alpha' in s.lower() else '#4b88ff' for s in self.pca_samples]
        ax.scatter(coords[:, 0], coords[:, 1], c=colors, s=150, edgecolors='white', alpha=0.8, picker=7)
        
        for i, sample in enumerate(self.pca_samples):
            ax.text(coords[i, 0], coords[i, 1] + 0.1, sample.split('_')[0], color='white', ha='center', fontsize=9)

        ax.tick_params(colors='white')
        self.p_canvas = FigureCanvasTkAgg(fig, master=self.pca_view)
        self.p_canvas.draw()
        self.p_canvas.get_tk_widget().pack(fill="both", expand=True, padx=40, pady=20)
        fig.canvas.mpl_connect('pick_event', self.on_pca_click)

    def on_pca_click(self, event):
        idx = event.ind[0]
        col = self.pca_samples[idx]
        det = (self.data[col] > 0.11).sum()
        self.update_info(f"--- SAMPLE INFO ---\nName: {col.split('.PG.Quantity')[0]}\nProteins Detected: {det:,}")

    def render_volcano_plot(self):
        for w in self.volcano_view.winfo_children(): w.destroy()
        ctk.CTkLabel(self.volcano_view, text="Interactive Volcano (Click dots for Gene Intel)", font=("Inter", 18, "bold")).pack(pady=20)

        fig, ax = plt.subplots(figsize=(7, 5), facecolor='#2b2b2b')
        ax.set_facecolor('#2b2b2b')

        y_vals = -np.log10(self.data['Raw_P_Value'].clip(lower=1e-300))
        x_vals = self.data['Log2FC']

        colors = []
        for tier in self.data['Significance_Tier']:
            if "Tier 1" in str(tier): colors.append('#ff4b4b')
            elif "Tier 2" in str(tier): colors.append('#4b88ff')
            elif "Tier 3" in str(tier): colors.append('#50c878')
            else: colors.append('#555555')

        ax.scatter(x_vals, y_vals, c=colors, s=40, alpha=0.7, picker=7, edgecolors='none')
        ax.tick_params(colors='white')

        self.v_canvas = FigureCanvasTkAgg(fig, master=self.volcano_view)
        self.v_canvas.draw()
        self.v_canvas.get_tk_widget().pack(fill="both", expand=True, padx=40, pady=20)
        fig.canvas.mpl_connect('pick_event', self.on_volcano_click)

    def on_volcano_click(self, event):
        idx = event.ind[0]
        row = self.data.iloc[idx]
        self.get_external_intel(row)

    def get_external_intel(self, row):
        raw_gene = str(row['PG.Genes'])
        gene_id = raw_gene.split('_')[0] if '_' in raw_gene else raw_gene
        uniprot_id = str(row.get('PG.ProteinGroups', '')).split(';')[0]
        
        self.update_info(f"Connecting to databases for {gene_id}...")
        threading.Thread(target=self._fetch_intel_thread, args=(row, gene_id, uniprot_id)).start()

    def _fetch_intel_thread(self, row, gene_id, uniprot_id):
        intel = ""
        # 1. Broad Search by Symbol FIRST
        try:
            url = f"https://mygene.info/v3/query?q={gene_id}&fields=summary,name,go&species=human"
            r = requests.get(url, timeout=5).json()
            
            # 2. Backup Search by Accession SECOND
            if not r.get('hits') and uniprot_id != 'nan':
                r = requests.get(f"https://mygene.info/v3/query?q={uniprot_id}&fields=summary,name,go&species=human", timeout=5).json()

            summary, full_name, location = "No summary found.", "Unknown Protein.", "Unknown"
            if r.get('hits') and len(r['hits']) > 0:
                hit = r['hits'][0]
                summary = hit.get('summary', summary)
                full_name = hit.get('name', full_name)
                
                # Bulldozer Location Logic
                go_dict = hit.get('go', {})
                cc_data = go_dict.get('CC', go_dict.get('cc', []))
                if isinstance(cc_data, dict): cc_data = [cc_data]
                loc_list = [str(item.get('term')).title() for item in cc_data if isinstance(item, dict) and 'term' in item]
                if loc_list:
                    location = ", ".join(list(dict.fromkeys(loc_list))[:2])

            intel += f"NAME: {full_name}\nLOCATION: {location}\n\nSUMMARY:\n{summary}\n\n"
        except Exception as e:
            intel += f"⚠️ DB Error: {e}\n\n"

        try:
            # 12-second Timeout for slow Literature servers
            r_pmc = requests.get(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={gene_id}&format=json&pageSize=2", timeout=12)
            if r_pmc.status_code == 200:
                data = r_pmc.json()
                papers = "".join([f"• {a.get('title')} ({a.get('pubYear')})\n\n" for a in data.get('resultList', {}).get('result', [])])
                intel += f"--- RECENT RESEARCH ---\n{papers or 'None found.'}"
        except:
            intel += "--- RECENT RESEARCH ---\n⚠️ Literature server timeout.\n"

        self.after(0, self.update_info, intel + f"\n--- STATS ---\nLFC: {row['Log2FC']:.2f}\nP-Val: {row['Raw_P_Value']:.2e}")

    def update_info(self, txt):
        self.info_box.configure(state="normal")
        self.info_box.delete("1.0", "end")
        self.info_box.insert("1.0", txt)
        self.info_box.configure(state="disabled")

if __name__ == "__main__":
    app = LabDiscoveryDatabase()
    app.mainloop()