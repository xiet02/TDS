##ref: https://microsoft.github.io/graphrag/get_started/

###step 0: install graphRAG
pip install graphrag

###step 1: make dir
mkdir -p ./TDS/input

###step 2: prepare reference file(s), example:
lynx -dump https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/PMC10290806/unicode
lynx -dump https://www.ncbi.nlm.nih.gov/research/bionlp/RESTful/pmcoa.cgi/BioC_json/PMC10251928/unicode 
...
###Depending on available omics analysis output, add them to the input folder. Example: 
### Safety_Score.txt #This data provides the core 0-28 safety metrics used for target prioritization.
### Convert the Safety_Scores.txt to JSON format for any 'text' fields, skipping noise sections

###step 3: initialize the workspace 
graphrag init --root ./TDS
#note:# Note: This script will create two files: .env and settings.yaml.
#Please refer to https://microsoft.github.io/graphrag/get_started/ to set up your own custom settings.

### Customizing Settings: Refer to the Get Started guide to configure the settings.yaml. For the TAA project, 
ensure the following prompt files are customized in the prompts/ directory:
  #local_search_system_prompt.txt: Tuned for Oncology Target Discovery.
  #extract_claims.txt: Configured to capture "Quantitative Assertions" like p-values and safety scores.
  #question_gen_system_prompt.txt: Updated to generate technical research questions for target triage.

###step 4: Run graphRAG
graphrag index --root ./TDS

#After completion, your indexing artifacts will be stored in parquet files within ./TDS/output.
#The duration of this process depends on your input data size, model choice, and text chunk size (configurable in settings.yaml). 

###step 5: Query the GraphRAG
Use natural language queries to evaluate targets based on validity, safety profiles, and therapeutic suitability.
##Local Search (Comparative Target Evaluation):
Bash
graphrag query \
--root ./TDS \
--method local \
--query "Evaluate TACSTD2 and MUC1 as Tumor-Associated Antigens in NSCLC; which one is a better ADC/TCE target?"

##Global Search (Thematic Dataset Synthesis):
Bash
graphrag query \
--root ./TDS \
--method global \
--query "Summarize the overarching safety risks for major TAAs across the clinical literature in this dataset."
