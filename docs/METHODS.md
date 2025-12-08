# Methods: UniProt IPR013902 Counting and Taxonomy Tree Generation

Custom Python scripts queried the UniProt REST API (UniProtKB release 2025_04) to enumerate, for each of the 1,430 fungal UniProtKB reference proteomes (`proteome_type:1`, `taxonomy_id:4751`), the number of proteins annotated with InterPro entry `IPR013902`. Reference proteomes were retrieved via the proteomes stream endpoint (`https://rest.uniprot.org/proteomes/stream`) with fields `upid`, `organism`, `organism_id`, `protein_count`, `busco`, `cpd`, `lineage`, `mnemonic`, `components`, `genome_assembly`, and `genome_representation`. For each proteome, counts of `IPR013902`-annotated proteins were obtained by querying the UniProtKB search endpoint (`https://rest.uniprot.org/uniprotkb/search`) with `proteome:UP000xxxxx AND database:interpro AND IPR013902` and reading the `x-total-results` header (`size=0`, canonical sequences; includes reviewed and unreviewed entries). All queries were executed on 2025-10-18. Results were aggregated into `fungal_ref_proteomes_with_interpro.tsv` and used to generate a hierarchical fungal taxonomy distribution tree saved as `taxonomy_tree_output.txt`. The scripts and outputs are available at `https://github.com/lilindu/cryo-tdk-killers`.

## Commands executed on 2025-10-18

```bash
curl "https://rest.uniprot.org/proteomes/stream?fields=upid%2Corganism%2Corganism_id%2Cprotein_count%2Cbusco%2Ccpd%2Clineage%2Cmnemonic%2Ccomponents%2Cgenome_assembly%2Cgenome_representation&format=tsv&query=%28%28proteome_type%3A1%29+AND+%28taxonomy_id%3A4751%29%29" > fungal_ref_proteomes.tsv

python add_interpro_counts.py

python taxonomy_tree.py fungal_ref_proteomes_with_interpro.tsv taxonomy_tree_output.txt
```

## Implementation references

- Counting query construction and use of `x-total-results`: `add_interpro_counts.py:69-75`
- Tree generation and CLI usage: `taxonomy_tree.py:693-731`