# EC Number Prediction 

This repository contains the code and data used to evaluate computational tools for enzyme function prediction. We benchmarked multiple EC (Enzyme Commission) number prediction algorithms — spanning both **similarity-based** and **machine/deep learning** approaches — using reaction SMILES as input. The pipeline includes scripts for dataset preprocessing, tool evaluation under various conditions, performance assessment across EC hierarchy levels and classes, and visualization of results. Overall, this repository provides a reproducible and extensible framework for benchmarking EC number prediction methods and helps users identify the most suitable tool for their metabolic modeling applications.

<p align="center">
  <img src="FIG1.png" alt="General pipeline of EC number prediction methodologies. A) Pipeline of similarity-based methods. B) Pipeline of machine/deep learning methods." width="600">
  <br>
  <em>General pipeline of EC number prediction methodologies. A) Similarity-based methods. B) Machine/deep learning methods.</em>
</p>


Specifically, we assessed the tools under three conditions: 
1. We evaluated all selected methods — [E-zyme](https://www.genome.jp/tools/e-zyme/), [E-zyme2](https://www.genome.jp/tools/e-zyme2/), [BridgIT](https://lcsb-databases.epfl.ch/Bridgit), [SelenzymeRF](https://github.com/synbiochem/selenzyme/tree/SelenzymeRF), [SIMMER](https://github.com/aebustion/SIMMER), [Theia](https://github.com/daenuprobst/theia), [BEC-Pred](https://github.com/KeeliaQWJ/BEC-Pred) and [CLAIRE](https://github.com/zishuozeng/CLAIRE) — Using 20% of the KEGG 2025 database (1866 reactions). 
2. For all of the methods with avaiable source code — [SelenzymeRF](https://github.com/synbiochem/selenzyme/tree/SelenzymeRF), [SIMMER](https://github.com/aebustion/SIMMER), [Theia](https://github.com/daenuprobst/theia), [BEC-Pred](https://github.com/KeeliaQWJ/BEC-Pred) and [CLAIRE](https://github.com/zishuozeng/CLAIRE) — we trained or used as prior knowledge 80% of the MetaNetX database (34.046 reactions), and then queried the methods with the remaining 20% (3.783 reactions). 
3. We did a case study on 28 drugs and their associated enzyme-annotated degradation reactions, and used them to query against all selected methods. Additionally, we applied a Top1 and Top5 **majority voting strategy** using [SelenzymeRF](https://github.com/synbiochem/selenzyme/tree/SelenzymeRF), [SIMMER](https://github.com/aebustion/SIMMER), [Theia](https://github.com/daenuprobst/theia) and [BEC-Pred](https://github.com/KeeliaQWJ/BEC-Pred), to show the potential of combining multiple algorithms to correctly predict EC number. 

For more information, please refer to:  

- Josefina Arcagni: jarcagniriv@unav.es 
- Telmo Blasco: tblasco@tecnun.es

### Table of contents: 

- [Project Structure](#project-structure)
- [Installation](#installation)
- [Included Tools](#included-tools)
  - [E-zyme / E-zyme2](#e-zyme--e-zyme2)
  - [BridgIT](#bridgit)
  - [SelenzymeRF](#selenzymerf)
  - [SIMMER](#simmer)
  - [Theia](#theia)
  - [BEC-Pred](#bec-pred)
  - [CLAIRE](#claire)
- [Results](#results)
  - [Case 1](#case-1)
  - [Case 2](#case-2)
  - [CaseStudy](#casestudy)
  - [MajorityVote](#majorityvote)
  
## Cite

## Project Structure
The code has the following structure: 

```
ECNumberPrediction/
├── data/            # Input datasets for analysis
├── methods/         # Implemented EC number prediction methods
│   ├── BEC-Pred
│   ├── BridgIT
│   ├── CLAIRE
│   ├── E-zyme
│   ├── SelenzymeRF
│   ├── SIMMER
│   └── theia
├── results/         # Output results from the analyses
│   ├── Case1
│   ├── Case2
│   ├── CaseStudy
│   └── MajorityVote
└── README.md
```

### Notes:
- Put input files under `data/`.
- Each method in `methods/` has its own **bash file** with implementation steps.
- Outputs for runs and evaluations go to `results/<method_or_case>/`.
- Example workflow:
    1. Prepare inputs in `data/`.
    2. Run a method from `methods/run_<MethodName>.sh/`.
    3. Check results in `results/`.


## Installation

Clone the repository:

```bash
git clone https://github.com/yourusername/ECNumberPrediction.git
cd ECNumberPrediction
```
Each method in the `methods/` folder may have its own installation requirements. Refer to the individual method documentation for setup instructions.


## Included Tools



| **Tool**      | **Year** | **Type** | **Database**     | **Features**                                                                                          | **Open-source code** |
|----------------|----------|----------|------------------|--------------------------------------------------------------------------------------------------------|----------------------|
| **E-zyme**     | 2009     | SB       | KEGG             | RDM patterns, substrate-product pairs, Tanimoto score                                                  | [No](https://www.genome.jp/tools/e-zyme/)                   |
| **E-zyme2**    | 2016     | SB       | KEGG             | RDM patterns, substrate-product pairs, graph-based substructures                                       | [No](https://www.genome.jp/tools/e-zyme2/)                   |
| **BridgIT**    | 2019     | SB       | KEGG             | Daylight fingerprints, reactive site identification, BNICE.ch rules                                    | [No](https://lcsb-databases.epfl.ch/Bridgit)                   |
| **SelenzymeRF**| 2023     | SB       | MetaNetX         | Morgan fingerprints, RXNMapper reactive sites, fragment analysis                                       | [Yes (GitHub)](https://github.com/synbiochem/selenzyme/tree/SelenzymeRF)         |
| **SIMMER**     | 2023     | SB       | MetaCyc          | Atom-Pair fingerprints, Tanimoto score, enrichment analysis                                            | [Yes (GitHub)](https://github.com/aebustion/SIMMER)         |
| **Theia**      | 2023     | ML       | ECREACT / Rhea   | MLP, differential reaction fingerprints                                     | [Yes (GitHub)](https://github.com/daenuprobst/theia)         |
| **BEC-Pred**   | 2024     | ML       | USPTO-ECREACT    | BERT, transfer learning                                                                                | [Yes (GitHub)](https://github.com/KeeliaQWJ/BEC-Pred)         |
| **CLAIRE**     | 2025     | ML       | ECREACT          | Contrastive learning, rxnfp embeddings, differential reaction fingerprints                             | [Yes (GitHub)](https://github.com/zishuozeng/CLAIRE)         |

*The table above summarizes the tools used, detailing their year of release, type (SB: similarity-based or ML: machine learning), associated databases, key features, and availability of open-source code.*

Follow the steps below for implementing each of the tools included in the `methods/` folder:

### E-zyme / E-zyme2

1. Create & activate venv:
```bash
python3 -m venv e-zyme_venv
source e-zyme_venv/bin/activate 
```

2. Install libs:
```bash
python -m pip install --upgrade pip
pip install requests beautifulsoup4 pandas
```

3. Make runner executable and run:
```bash
chmod +x methods/E-zyme/run_ezyme.sh
bash methods/E-zyme/run_ezyme.sh
```

This will run the script for Case 1 and Case Study and the *result* csv files fr E-zyme and E-zyme 2 will both be saved in the `results/` folder under the corresponding case folders with the tool name: `E-zyme.csv`. 

### BridgIT

1. Create & activate venv:
```bash
python3 -m venv bridgit_venv
source bridgit_venv/bin/activate 
```
2. Install libs:
```bash
python -m pip install --upgrade pip
pip install requests pandas rdkit
```

3. Make runner executable and run:
```bash
chmod +x methods/BridgIT/run_BridgIT.sh
bash methods/BridgIT/run_BridgIT.sh
```

*Notes:* BridgIT can only be accessed through its own server in https://lcsb-databases.epfl.ch/Bridgit, a user account needs to be created to access it. The first part of the bash file will process reaction SMILES and create the necessary files in `methods/BridgIT/input/` to input in their web server. Once the results ready, download them from the server and place them in the `methods/BridgIT/output/` folder. The steps for extracting the results are in the second part of the bash file. 

This will run the script for Case 1 and Case Study and the *result* csv files will be saved in the `results/` folder under the corresponding case folders with the tool name: `BridgIT.csv`. 

### SelenzymeRF

1. Create & activate venv:
```bash
python3 -m venv selenzyme_venv
source selenzyme_venv/bin/activate 
```

2. Install libs:
```bash
python -m pip install --upgrade pip
pip install requests beautifulsoup4 pandas
```

3. Clone repository & download data:
```bash
git clone -b SelenzymeRF https://github.com/synbiochem/selenzyme.git methods/SelenzymeRF/SelenzymeRF_code
cd methods/SelenzymeRF/SelenzymeRF_code
```
*Notes:* unzip the required reference datasets as indicated in the upstream repository. The script already included in `Selenzyme_code/`: `/methods/SelenzymeRF/SelenzymeRF_code/start_server.sh` was specifically modified to run the cases and shouldn´t be deleted. 


4. Make runner executable and run:
```bash
chmod +x /methods/SelenzymeRF/run_selenzymerf.sh
bash /methods/SelenzymeRF/run_selenzymerf.sh
```

This will run the script for all cases and the *result* csv files will be saved in the `results/` folder under the corresponding case folders with the tool name: `SelenzymeRF.csv`.

### SIMMER

1. Create conda environment:
```bash
conda env create -f methods/SIMMER/simmer_env.yml
conda activate simmer_env
```

2. Clone repository & download data:
```bash
git clone https://github.com/aebustion/SIMMER.git methods/SIMMER/SIMMER_code
cd methods/SIMMER/SIMMER_code
```

*Notes:* import the required reference datasets as indicated in the upstream repository. The scripts already included in `SIMMER_code/`: `methods/SIMMER/SIMMER_code/SIMMER/SIMMER.py` and `methods/SIMMER/SIMMER_code/SIMMER/SIMMER2.py` were specifically modified to run the cases and shouldn´t be deleted. 

3. Make runner executable and run:
```bash
chmod +x methods/SIMMER/run_SIMMER.sh
bash methods/SIMMER/run_SIMMER.sh
```

This will run the script for all cases and the *result* csv files will be saved in the `results/` folder under the corresponding case folders with the tool name: `SIMMER.csv`.

### Theia 

*Requirements*: 1 GPU and 50G memory for training and evaluating the model.

1. Create conda environment:
```bash
conda env create -f methods/theia/theia_env.yml
conda activate theia_env
```

2. Clone repository & download data:
```bash
git clone https://github.com/daenuprobst/theia.git methods/theia/theia_code
cd methods/theia/theia_code
```

*Notes:*  All the scripts already included in `Theia_code/`  were specifically modified to run the cases and shouldn´t be deleted. 


3. Make runner executable and run:
```bash
chmod +x methods/theia/run_theia.sh
bash methods/theia/run_theia.sh
```

This will run the script for all cases and the *result* csv files will be saved in the `results/` folder under the corresponding case folders with the tool name: `theia.csv`.

### BEC-Pred 

*Requirements*: 1 GPU and 128G memory for training and evaluating the model.

1. Create conda environment:
```bash
conda env create -f methods/BEC-Pred/becpred_gpu.yml
conda activate becpred_env
```

2. Make runner executable and run:
```bash
chmod +x methods/BEC-Pred/run_becpred.sh
bash methods/BEC-Pred/run_becpred.sh
```

This will run the script for all cases and the *result* csv files will be saved in the `results/` folder under the corresponding case folders with the tool name: `BEC-Pred.csv`.

### CLAIRE

*Requirements:* 1 GPU and 128G memory for training and evaluating the model.

1. Create conda environments:
```bash
conda env create -f methods/CLAIRE/claire_env.yml
# activate when using CLAIRE tools
conda activate claire_env

conda env create -f methods/CLAIRE/rxnfp_env.yml
# activate when using rxnfp utilities
conda activate rxnfp_env
```

2. Clone repository & download data:
```bash
git clone https://github.com/zishuozeng/CLAIRE.git methods/CLAIRE/CLAIRE_code
cd methods/CLAIRE/CLAIRE_code
```
Notes: Import or unzip any reference datasets required by CLAIRE as indicated in the upstream repository. Do not remove or modify the modified scripts included in `methods/CLAIRE/CLAIRE_code/` that are needed to run the cases.

3. Make runner executable and run:
```bash
chmod +x methods/CLAIRE/run_claire.sh
bash methods/CLAIRE/run_claire.sh
```

This will run CLAIRE for the configured cases and save result CSV files under the corresponding case folders in `results/` with the tool name: `CLAIRE.csv`.

## Results

Results are organized in subfolders inside `results/`:

- `results/Case1/` – Results for the first evaluation case (queried all methods with their original dataset tested with KEGG reaction queries).

- `results/Case2/` – Results for the second evaluation case (five open-code methods trained on 80% of MetaNetX dataset and queried on the other 20%).

- `results/CaseStudy/` – Results for the Case Study (queries all methods with their original dataset with 28).

- `results/MajorityVote/` – Top1 and Top5 majority voting strategies using [SelenzymeRF](https://github.com/synbiochem/selenzyme/tree/SelenzymeRF), [SIMMER](https://github.com/aebustion/SIMMER), [Theia](https://github.com/daenuprobst/theia) and [BEC-Pred](https://github.com/KeeliaQWJ/BEC-Pred). 

To reproduce the metrics and figures used in the paper for each case, follow the steps below:

### Case 1

1. Merge all tool CSV outputs into merged_output.csv:
```bash
python3 results/Case1/join_results.py
```
Output: results/Case1/merged_output.csv

2. Compute evaluation metrics:
```bash
python3 results/Case1/get_metrics.py
```
Output: results/Case1/evaluation_summary.csv

3. Generate plots:
```bash
Rscript results/Case1/metrics-case1.R
```
The plot with be saved in ´results/Case1/case1_plot.png´.

### Case 2

1. Merge all tool CSV outputs into merged_output.csv:
```bash
python3 results/Case2/join_results.py
```
Output: results/Case2/merged_output.csv

2. Compute evaluation metrics:
```bash
python3 results/Case2/get_metrics.py
```
Output: results/Case2/evaluation_summary.csv

3. Generate plots:
```bash
Rscript results/Case2/metrics-case2.R
```
The plot with be saved in ´results/Case2/case2_plot.png´.

### CaseStudy
1. Merge all tool CSV outputs into merged_output.csv:
```bash
python3 results/CaseStudy/join_results.py
```
3. Generate plots:
```bash
Rscript results/CaseStudy/case_study_heatmap.R
```
Output: results/CaseStudy/merged_output.csv

The plot with be saved in ´results/CaseStudy/casestudyplot.png´.

### MajorityVote 

A tool to combine multiple EC number prediction outputs using two Majority Voting strategies. 
- **Top1**: Takes the first (top) EC prediction from each method and picks the EC number that appears most frequently across all methods.
- **Top5**: Takes the Top 5 predictions from each method. Each prediction contributes a weighted vote (5 points for rank 1, 4 for rank 2, etc.), the EC number with the highest cumulative weighted score across methods is selected as the consensus prediction.

#### Input: 
A CSV file containing:
- An identifier column (default: entity, e.g., enzyme, compound, or reaction ID)
- One or more method columns, each containing EC predictions in this format:

```text
1.1.1.1|1.1.1.2;2.7.1.1|2.7.1.2
```
Where:
- ; separates ranked prediction groups (rank 1 → rank 2 → rank 3 …). 
- | separates tied ECs within the same rank.

#### Usage

Command-Line Examples: 

1. Compute for All Entities 

```bash
python3 majority_ec_vote.py \
  --input_csv results/merged_predictions.csv \
  --entity all \
  --use_all \
  --output_csv results/majority_votes.csv
```

2. Compute for One Specific Reaction
```bash
python3 majority_ec_vote.py \
  --input_csv results/merged_predictions.csv \
  --entity R001 \
  --methods MethodA MethodB MethodC
```

3. Custom Identifier Column
```bash
python3 majority_ec_vote.py \
  --input_csv results/merged_predictions.csv \
  --entity all \
  --use_all \
  --id_col reaction_id
```

#### Output: 

For each entity (or a single entity if specified), the tool reports:

| **Column Name**       | **Description** |
|:----------------------|:----------------|
| `majority_vote_top1`  | Most frequent Top-1 EC prediction across methods |
| `majority_vote_top5`  | Weighted Top-5 consensus EC prediction |

(Optional) saved to a CSV file if --output_csv is provided	
