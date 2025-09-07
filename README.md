# Supporting R Code for Manuscript: "Fungal diversity drives soil health and agricultural sustainability in black soils"

This repository contains six R scripts that support the analyses presented in **Yang et al.**'s manuscript titled *"Fungal diversity drives soil health and agricultural sustainability in black soils"*. These scripts cover data preprocessing, network construction, motif analysis, and taxonomy assignment.

## Repository Contents

### 5. `Code1_taxonomy_assignment.r`
- **Purpose**: Assigns taxonomy to OTUs using multiple datasets.
- **Key Functions**:
  - Merges taxonomy assignments from different sources (e.g., ITS, 18SV9).
  - Selects the most informative assignment based on scoring.
- **Output**: 
  - Final merged taxonomy assignment (`final_merged_data.csv`).

### 1. `Code2_network_construction.r`
- **Purpose**: Constructs microbial networks from OTU tables using the SpiecEasi method.
- **Key Functions**:
  - Filters OTUs based on abundance and prevalence thresholds.
  - Rarefies OTU tables to an even depth.
  - Constructs microbial co-occurrence networks.
- **Output**: 
  - Rarefied OTU table (`rarefied_otu_table.csv`).
  - Network edge list and graph structure.

### 2. `Code3_qcmi.r`
- **Purpose**: Analyzes microbial network cohesion and ecological associations based on environmental factors.
- **Key Functions**:
  - Splits environmental factors into soil, climate, and plant categories.
  - Evaluates network properties influenced by ecological factors.
  - Calculates cohesion metrics using `qcmi`.
- **Output**: 
  - CSV files for soil, climate, plant, and biotic network edges.
  - Network cohesion results (`re_qcmi1.csv`).

### 3. `Code4_link_count_z.r`
- **Purpose**: Calculates interaction counts between different node types and performs Z-score and P-value calculations.
- **Key Functions**:
  - Constructs subgraphs for each sample.
  - Calculates interaction counts between node types (e.g., Fungi, Protists, Metazoa, Plants).
  - Computes Z-scores and P-values against random network distributions.
- **Output**: 
  - Interaction count matrix (`df_total.csv`).
  - Z-score and P-value results (`z_scores.csv`, `p_values.csv`).

### 4. `Code5_motifs.r`
- **Purpose**: Performs motif analysis on microbial networks.
- **Key Functions**:
  - Identifies network motifs using the `netmotif` function.
  - Compares observed motifs with random networks to calculate Z-scores and P-values.
- **Output**: 
  - Motif analysis results (`result_netmotif.csv`).
  - Z-scores and P-values for motifs (`motif_z_scores.csv`, `motif_p_values.csv`).

### 6. `network_motif_functions.r`
- **Purpose**: Defines functions for network motif analysis.
- **Key Functions**:
  - Implements various motif calculations, including cycle facilitation, cycle competition, and transitive associations.
- **Usage**: Functions are sourced in `Code4_motifs.r` to calculate motif distributions.

## Usage Instructions

1. Clone this repository:
   ```bash
   git clone https://github.com/joshualiuxu/Diversity-supports-SOC-and-productivity.git




For questions or issues regarding the code, please contact:

Xu Liu
Email: xliu@issas.ac.cn

Teng Yang
Email: tyang@issas.ac.cn   
