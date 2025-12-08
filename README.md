# IPR013902 Analysis

This repository contains scripts and outputs for enumerating the distribution of InterPro entry `IPR013902`-annotated proteins across fungal UniProtKB reference proteomes and rendering a hierarchical taxonomy tree for reporting.

## Data Sources
- UniProtKB release: `2025_04`
- Proteomes stream endpoint: `https://rest.uniprot.org/proteomes/stream`
- UniProtKB search endpoint: `https://rest.uniprot.org/uniprotkb/search`
- Scope: `proteome_type:1` (reference proteomes), `taxonomy_id:4751` (Fungi), InterPro `IPR013902`
- Query date: `2025-10-18`

## Quickstart (Reproduce)
```bash
# 1) Download fungal reference proteomes
curl "https://rest.uniprot.org/proteomes/stream?fields=upid%2Corganism%2Corganism_id%2Cprotein_count%2Cbusco%2Ccpd%2Clineage%2Cmnemonic%2Ccomponents%2Cgenome_assembly%2Cgenome_representation&format=tsv&query=%28%28proteome_type%3A1%29+AND+%28taxonomy_id%3A4751%29%29" > fungal_ref_proteomes.tsv

# 2) Count IPR013902 proteins per proteome (creates fungal_ref_proteomes_with_interpro.tsv)
python add_interpro_counts.py

# 3) Generate taxonomy tree output
python taxonomy_tree.py fungal_ref_proteomes_with_interpro.tsv taxonomy_tree_output.txt
```

## Repository Structure
- `docs/METHODS.md`: Detailed methods paragraph and exact commands
- `fungal_ref_proteomes.tsv`: Fungal reference proteomes list (downloaded 2025-10-18)
- `fungal_ref_proteomes_with_interpro.tsv`: IPR013902 counts per proteome
- `add_interpro_counts.py`: UniProtKB counting (uses `x-total-results`)
- `taxonomy_tree.py`: Hierarchical taxonomy rendering for IPR013902 distribution
- `taxonomy_tree_output.txt`: Final tree used in the paper

## Results
- The taxonomy tree output is available in `taxonomy_tree_output.txt` on the feature branch above.