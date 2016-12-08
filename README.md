# Annotations for each scripts

## Raw data processing and filtering

* `filter_probe_vs_probe.R`: ensure every protein is a target for at least one of the two probes

* `filter_20uM_competition.R`: ensure every protein is enriched by corresponding probe with SILAC ratio greater than 5

* `filter_200uM_competition.R`: ensure every protein is enriched by corresponding probe with SILAC ratio greater than 5

* `filter_competed.R`: ensure every protein is competed (SILAC ratio greater than 3) by at least one of tested competitors

* `filter_miscleavaged_peptide.R`: remove miscleavaged tryptic peptide in isoTOP-ABPP 

* `filter_PDB.R`: filter PDB files for structural analysis

## Meta-analysis and cross-referencing

* `mapping_drugbank.R`: map probe targets in [Drugbank](https://www.drugbank.ca/) database

* `mapping_iBAQ.R`: map probe targets in [iBAQ](http://www.mcponline.org/content/11/3/M111.014050.long) database

* `reverse_mapping_iBAQ.R`: map protein entries in iBAQ dataset to probe targets dataset

* `mapping_protein_family_separate.R`: assign protein families to target list

* `get_transmembrane.R`: extract soluble/membrane information from Uniprot databse and assign to targets

* `get_functionalSite.py`: extract functional site information from Uniprot database

* `functionalSite_analysis.R`: analyze spatial distance between functional sites fetched from Uniprot and probe labeled peptides

* `fpocket_analysis.R`: analyze fpocket results to obtain the overlap between fpocket predicted binding sites and probe labeled peptides

# External databases for cross-reference

* `drugbank_all_target_polypeptide_ids.fasta`: drug target identifiers (All group, version 4.2), from [Drugbank](https://www.drugbank.ca/)

* `drugbank_nutraceutical_target_polypeptide_ids.fasta`: drug target identifiers (Nutraceutical group, version 4.2), from [Drugbank](https://www.drugbank.ca/)

* `iBAQ.csv`: iBAQ values in HEK293t cell line, from Supplemental Data of [Geiger, T., Wehner, A., Schaab, C., Cox, J., and Mann, M. (2012). Mol Cell Proteomics 11, M111.014050.](http://www.mcponline.org/content/11/3/M111.014050.long)