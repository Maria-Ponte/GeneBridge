import customtkinter
import tkinter
from tkinter import ttk
import euclidean_distance as euc
import bray_curtis_distance as bc
import cossine_similarity as cos
import manhattan_distance as man
import canberra_distance as can
import shannons_entropy as shan
import lempel_ziv as lempel_ziv
import kullback_leibler as kl
import chaos_game_representation as cgr

def parse_fasta(fasta_string):
    names = []
    sequences = []
    
    lines = fasta_string.split('\n')
    current_name = ''
    current_sequence = ''
    
    for line in lines:
        line = line.strip()
        
        if line.startswith('>'):
            # Header line, extract the name/identifier
            current_name = line[1:]
            
            # Store the previous sequence, if any
            if current_sequence:
                sequences.append(current_sequence)
                current_sequence = ''
            
            names.append(current_name)
        else:
            # Sequence line, append to the current sequence
            current_sequence += line
    
    # Store the last sequence
    if current_sequence:
        sequences.append(current_sequence)
    
    return names, sequences



class MyScrollableCheckboxFrame(customtkinter.CTkScrollableFrame):
    def __init__(self, master, title, values, bg, height):
        super().__init__(master, label_text=title, label_fg_color= bg, label_text_color="whitesmoke")
        self.grid_columnconfigure(0, weight=1)
        self.values = values
        self.checkboxes = []
        

        for i, value in enumerate(self.values):
            checkbox = customtkinter.CTkCheckBox(self, text=value)
            checkbox.grid(row=i+1, column=0, padx=10, pady=(10, 0), sticky="w")
            self.checkboxes.append(checkbox)

        self.configure(height = height)

    def get(self):
        checked_checkboxes = []
        for checkbox in self.checkboxes:
            if checkbox.get() == 1:
                checked_checkboxes.append(checkbox.cget("text"))
        return checked_checkboxes


class MyScrollableCheckboxFrame2(customtkinter.CTkScrollableFrame):
    def __init__(self, master, title2, values2, bg):
        super().__init__(master, label_text=title2, label_fg_color= bg, label_text_color="whitesmoke")
        self.grid_columnconfigure(0, weight=1)
        self.values2 = values2
        self.checkboxes2 = []

        for i, value in enumerate(self.values2):
            checkbox = customtkinter.CTkCheckBox(self, text=value)
            checkbox.grid(row=i+1, column=0, padx=10, pady=(10, 0), sticky="w")
            self.checkboxes2.append(checkbox)
            
   
    def get(self):
        checked_checkboxes_2 = []
        for checkbox in self.checkboxes2:
            if checkbox.get() == 1:
                checked_checkboxes_2.append(checkbox.cget("text"))
        return checked_checkboxes_2
    
class MyScrollableCheckboxFrame3(customtkinter.CTkScrollableFrame):
    def __init__(self, master, title3, values3, bg):
        super().__init__(master, label_text=title3, label_fg_color= bg, label_text_color="whitesmoke")
        self.grid_columnconfigure(0, weight=1)
        self.values3 = values3
        self.checkboxes3 = []

        for i, value in enumerate(self.values3):
            checkbox = customtkinter.CTkCheckBox(self, text=value)
            checkbox.grid(row=i+1, column=0, padx=10, pady=(10, 0), sticky="w")
            self.checkboxes3.append(checkbox)
            
   
    def get(self):
        checked_checkboxes_3 = []
        for checkbox in self.checkboxes3:
            if checkbox.get() == 1:
                checked_checkboxes_3.append(checkbox.cget("text"))
        return checked_checkboxes_3

class App(customtkinter.CTk):
    def __init__(self):
        super().__init__()

        self.title("Gene Bridge")
        self.geometry("700x630")
        self.grid_columnconfigure(0, weight=1)
        self.grid_rowconfigure(1, weight=1)

        self.label = customtkinter.CTkLabel(self, text=" Sequences" ,font=("Arial",15) , text_color= '#388583')
        self.label.grid(row=0, column=0, padx=10, pady=(5, 0), sticky="w")

        self.label = customtkinter.CTkLabel(self, text=" Sequences must be written in Fasta Format")
        self.label.grid(row=1, column=0, columnspan= 3, padx=10, pady=(0, 0), sticky="nw")

        self.textbox = customtkinter.CTkTextbox(self)
        self.textbox.grid(row=2, column=0, columnspan = 3, padx=10, pady=10, sticky="nsew")

        values = ["[WB] Euclidean Distance" , "[WB] Normalized Euclidean Distance", "[WB] Bray-Curtis Distance", "[WB] Normalized Bray-Curtis Distance", "[WB] Cosine Similarity", "[WB] Manhattan Distance" , "[WB] Normalized Manhattan Distance" , "[WB] Canberra distance", "[WB] Normalized Canberra distance", "[IB] Shannon's Entropy", "[IB] Lempel-Ziv Complexity", "[IB] Kullback-Leibler Divergence", "[GR] Chaos Game Representation"]
        self.scrollable_checkbox_frame = MyScrollableCheckboxFrame(self, title="Methods", values=values, bg='#388583', height = 10)
        self.scrollable_checkbox_frame.grid(row=3, column=0, padx=10, pady=(10, 0), sticky="ew")
        
        values3 = ["single", "complete", "average", "weighted", "centroid", "median", "ward"]
        self.scrollable_checkbox_frame3 = MyScrollableCheckboxFrame3(self, title3="Hierarchical Cluster", values3=values3, bg='#388583')
        self.scrollable_checkbox_frame3.grid(row=3, column=2, padx=10, pady=(10, 0), sticky="e")
        
        self.label = customtkinter.CTkLabel(self, text=" [WB] word-based methods    ||    [IB] info-based methods    ||    [GR] graphical representation")
        self.label.grid(row=5, column=0, columnspan= 3, padx=10, pady=(5, 0), sticky="w")


        self.button = customtkinter.CTkButton(self, text="Start", command=self.button_callback, fg_color='#388583')
        self.button.grid(row=6, column=0, columnspan = 3, padx=10, pady=10, sticky="ewn")

        

    def button_callback(self):
        names,seq = parse_fasta(self.textbox.get("1.0","end-1c"))
        clusters = self.scrollable_checkbox_frame3.get()
        methods = self.scrollable_checkbox_frame.get()
        n_methods = len(methods)
        
        
        print("Labels:", names)
        print("Sequences:" , seq)
        print("Hierarquical Clusters:", clusters)
        print("Methods:", methods)
        
    

        for i in range(1, n_methods + 1):
            if methods[i - 1] == "[WB] Euclidean Distance":
                print(euc.plot_eucl(seq, names,clusters))

            elif methods[i - 1] == "[WB] Normalized Euclidean Distance":
                print(euc.plot_eucl_norm(seq, names, clusters))

            elif methods[i - 1] == "[WB] Bray-Curtis Distance":
                print(bc.plot_bc(seq, names,clusters))

            elif methods[i - 1] == "[WB] Normalized Bray-Curtis Distance":
                print(bc.plot_bc_norm(seq, names,clusters))

            elif methods[i - 1] == "[WB] Cosine Similarity":
                print(cos.plot_cos(seq, names,clusters))

            elif methods[i - 1] == "[WB] Manhattan Distance":
                print(man.plot_man(seq, names,clusters))

            elif methods[i - 1] == "[WB] Normalized Manhattan Distance":
                print(man.plot_man_norm(seq, names,clusters))

            elif methods[i - 1] == "[WB] Canberra distance":
                print(can.plot_can(seq, names,clusters))

            elif methods[i - 1] == "[WB] Normalized Canberra distance":
                print(can.plot_can_norm(seq, names,clusters))

            elif methods[i - 1] == "[IB] Shannon's Entropy":
                print(shan.plot_shan(seq, names, clusters))

            elif methods[i - 1] == "[IB] Lempel-Ziv Complexity":
                print(lempel_ziv.plot_lempel_ziv(seq, names, clusters))

            elif methods[i - 1] == "[IB] Kullback-Leibler Divergence":
                print(kl.plot_kl(seq, names, clusters))

            elif methods[i - 1] == "[GR] Chaos Game Representation":
                print(cgr.plot_cgr(seq, names))


app = App()
app.mainloop()




