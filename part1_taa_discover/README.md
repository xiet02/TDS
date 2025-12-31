# Setup Guide

> Reference: https://microsoft.github.io/graphrag/get_started/

## Step 0: Install GraphRAG
```bash
pip install graphrag
```

## Step 1: Create Directory Structure
```bash
mkdir -p ./TDS/input
```

## Step 2: Prepare Reference Files

Retrieve biomedical literature and add omics analysis outputs to the input folder.

### Example: Fetching PubMed Central Articles
```bash
lynx -dump https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/PMC10290806/unicode
lynx -dump https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/PMC10251928/unicode
```

### Adding Omics Analysis Output

Depending on available omics analysis output, add them to the input folder.

**Example:**
- `Safety_Score.txt` - Core 0-28 safety metrics for target prioritization
- Convert `Safety_Scores.txt` to JSON format for any 'text' fields, skipping noise sections

## Step 3: Initialize the Workspace
```bash
graphrag init --root ./TDS
```

**Note:** This script will create two files: `.env` and `settings.yaml`.

Please refer to the [Get Started guide](https://microsoft.github.io/graphrag/get_started/) to set up your own custom settings.

### Customizing Settings

Refer to the Get Started guide to configure `settings.yaml`. For the TAA project, ensure the following prompt files are customized in the `prompts/` directory:

- **`local_search_system_prompt.txt`** - Tuned for Oncology Target Discovery
- **`extract_claims.txt`** - Configured to capture "Quantitative Assertions" like p-values and safety scores
- **`question_gen_system_prompt.txt`** - Updated to generate technical research questions for target triage

## Step 4: Run GraphRAG Indexing
```bash
graphrag index --root ./TDS
```

After completion, your indexing artifacts will be stored in parquet files within `./TDS/output`.

The duration of this process depends on your input data size, model choice, and text chunk size (configurable in `settings.yaml`).

## Step 5: Query the GraphRAG

Use natural language queries to evaluate targets based on validity, safety profiles, and therapeutic suitability.

### Local Search (Comparative Target Evaluation)
```bash
graphrag query \
  --root ./TDS \
  --method local \
  --query "Evaluate TACSTD2 and MUC1 as Tumor-Associated Antigens in NSCLC; which one is a better ADC/TCE target?"
```

### Global Search (Thematic Dataset Synthesis)
```bash
graphrag query \
  --root ./TDS \
  --method global \
  --query "Summarize the overarching safety risks for major TAAs across the clinical literature in this dataset."
```
