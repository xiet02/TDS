# TADA: Target and Antibody Discovery via AI

**TADA** is an AI-powered pipeline for accelerating **ADC (Antibody-Drug Conjugate)** and **TCE (T-cell Engager)** drug development, covering both **target discovery** and **antibody engineering**.

---

## Overview

TADA integrates multi-source data and AI to streamline the discovery and optimization of therapeutic antibodies. The repository is organized into two parts:

- **Part I:** Antigen target identification and prioritization.
- **Part II:** In silico antibody design and developability optimization.

---

## Part I: TAA Discovery System

### Overview
The TAA Discovery System is an innovative framework for identifying and prioritizing tumor-associated antigens for TCE or ADC development. It integrates diverse datasets—including multi-omics repositories and scientific publications—using a graph retrieval-augmented generation (RAG)-enhanced language model to extract insights from biological/clinical literature and curated oncology-related omics databases (TCGA, GTEx, single-cell atlases, etc.).

### Features
- **Graph RAG-enhanced Language Model:** Efficiently extracts relevant information from scientific literature.
- **Omics Database Integration:** Combines data from TCGA, GTEx, and other sources for comprehensive analysis.
- **Safety Score:** Prioritizes TAAs with high tumor selectivity and low on-target/off-tumor risk.

### Reference
For methodology and findings, see:
[Xie & Hung,bioRxiv 2025](https://www.biorxiv.org/content/10.1101/2025.05.06.652559v1)

---

## Part II: In Silico Antibody Design

### Overview
After target nomination, antibody engineering remains a major bottleneck. TADA provides an integrated, AI-driven framework for in silico antibody design and multi-objective optimization, combining sequence diversification, structure prediction, docking confidence assessment, and developability scoring.

### Features
- **Sequence Diversification:** Generates CDR-focused antibody variants.
- **Structure Prediction:** Uses AlphaFold-Multimer and AlphaFold3 for accurate modeling.
- **Developability Scoring:** Evaluates solubility, CDR liability risks, and interface confidence.
- **Unbiased Ranking:** Produces a unified score for candidate prioritization.

### Reference
For methodology and findings, see:
[Xie, bioRxiv 2025](https://www.biorxiv.org)

---

## Reference & Contact
For questions or collaboration, contact: [txie@neoomics.com](mailto:txie@neoomics.com)
