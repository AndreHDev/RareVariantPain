# Code for: Rare variants in the general population as a bioinformatic gateway to pain genetics and analgesic drug discovery

This contains R and Python code for finding pain relevant papers in rare variant analysis and identifying clinically approved drugs targeting pain-related genes through OMIM-DrugBank integration.

## Project Structure

- **Pubmed/** - PubMed search and API integration
  - `PubMedAPI.py` - PubMed API wrapper
  - `Caching.py` - Caching utilities
  - `citedIn.py` - The Concrete Citation analysis
  
- **OMIM/** - OMIM data integration
  - `omim_drugbank.R` - OMIM and DrugBank integration
  - `create_variant_scores_plot.R` - Visualization scripts

## Citation

Kringel D, Lötsch J, Himmelspach A. Rare variants in the general population as a bioinformatic gateway to pain genetics and analgesic drug discovery. (submitted)


