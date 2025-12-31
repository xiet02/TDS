# AI-Powered Antibody Design and Developability Optimization

This repository provides an integrated, AI-driven framework for in silico antibody design and multi-objective optimization. It combines sequence diversification, structural prediction, docking confidence assessment, and developability scoring to accelerate the development of Antibody-Drug Conjugates (ADCs) and T-cell engagers (TCEs).

The framework uses tezepelumab (anti-TSLP), an FDA-approved antibody, as a benchmark to demonstrate an unbiased ranking system where clinical antibodies serve as internal controls.

## 📂 Repository Structure

The pipeline is organized into two primary task modules, reflecting the workflow of structural prediction followed by functional docking and prioritization:

### 1. Task 1-Developability Score

Focuses on variant generation and the assessment of biophysical fitness.

* **`step0_generate_library.py`**: Generates a library of 100 CDR-focused variants while maintaining framework conservation to ensure structural integrity.

* **`step1.1_batch_af2_library.py`**: Executes batch structure predictions using AlphaFold2 via ColabFold.

* **`step1.2_parse_AF2_output.py`**: Parses AlphaFold2 output to extract confidence metrics such as pLDDT (local confidence) and pTM (global fold confidence).

* **`step2_solubility_proxy.py`**: Calculates solubility scores based on charge distribution and hydrophobic surface area.

* **`step3_liability_scoring_cdr.py`**: Screens CDR sequences for chemical liabilities, including N-glycosylation sites, deamidation hotspots, and isomerization sites.

* **`step4_generate_CDS.py`**: Integrates metrics into a Composite Developability Score (CDS) using a weighted linear combination of solubility, sequence liabilities, and immunogenicity.

### 2. Task 2-AF2 Docking

Focuses on structural binding validation and final candidate prioritization. Note this will be replaced with an AF3 based pipeline.

* **`step0_convert_fasta_to_csv.py`**: Prepares sequence data for high-throughput processing by splitting scFv sequences into VH and VL domains.

* **`step1_make_multimer_fastas.py`**: Formats antibody and antigen sequences for complex prediction.

* **`step2_run_colabfold_loop.py`**: Automates the docking of antibody variants with target antigens using AlphaFold-Multimer or AlphaFold3.

* **`step3_parse_af2_and_interface.py`**: Extracts interface quality metrics, primarily iPTM, which specifically evaluates the multi-chain interface quality.

* **`step4_merge_dev_docking_rank.py`**: Produces the final decision-ready candidate list by balancing docking interface confidence with developability scores.

## 🚀 Key Features

* **Multi-Objective Optimization**: Simultaneously optimizes competing objectives, such as binding interface quality and developability parameters including solubility and stability.

* **Advanced Structural Modeling**: Leverages AlphaFold2/3 for structure prediction and docking, incorporating diffusion-based generation for improved antibody-antigen complex modeling.

* **Scalable Pipeline**: Accelerates early-stage optimization from traditional 2–4 year cycles to approximately 1–2 weeks, prioritizing the most promising candidates for experimental testing.

* **Unbiased Ranking**: Validates the system using clinical benchmarks; for example, tezepelumab (var0000) ranked #38, proving the framework discovers novel candidates based on computational profiles rather than sequence familiarity.

## 🛠️ Requirements

* Python 3.x
* **ColabFold/AlphaFold2**: For structural and multimer predictions.
* **Biopython & Pandas**: For sequence manipulation and data aggregation.
* **NetMHCIIpan**: For predicting T-cell epitopes and immunogenicity.

## 📖 Citation

If you use this framework in your research, please cite:

*AI-powered integration of multi-source data for TAA discovery to accelerate ADC and TCE drug development (II): in silico Antibody Design and Developability Optimization.*