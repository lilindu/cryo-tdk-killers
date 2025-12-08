#!/usr/bin/env python3
"""
Taxonomic Tree Generator for IPR013902 Domain Distribution

This script processes fungal proteome data to generate a hierarchical tree showing 
the distribution of IPR013902-containing proteins across different taxonomic levels.
The script analyzes multiple data sources including:
- Taxonomic lineage information
- Organism names for genus/species extraction
- IPR013902 protein counts per proteome
- Proteome identifiers for species tracking

Usage:
    python taxonomy_tree.py [input_file] [output_file]
    python taxonomy_tree.py -h/--help

Arguments:
    input_file: Path to the TSV file containing fungal proteome data with IPR013902 counts
    output_file: Optional output file path (if not provided, prints to stdout)
"""

import argparse
import csv
import sys
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Set


class TaxonomyNode:
    """
    Represents a node in the taxonomic tree.
    
    Each node contains:
    - name: taxonomic name
    - level: taxonomic level (0=root, 1=kingdom, etc.)
    - children: child nodes
    - protein_count: number of IPR013902-containing proteins
    - total_proteomes: total number of proteomes at this level
    """
    
    def __init__(self, name: str, level: int = 0, organism_name: str = None, lineage: str = None):
        self.name = name
        self.level = level
        self.organism_name = organism_name  # Store full organism name for species
        self.lineage = lineage  # Store lineage for positional genus detection
        self.children: Dict[str, 'TaxonomyNode'] = {}
        self.protein_count = 0
        self.total_proteomes = 0
        self.proteome_ids: Set[str] = set()
        self.species_with_proteins: Set[str] = set()  # Track species containing IPR013902 proteins
    
    def add_child(self, name: str, organism_name: str = None, lineage: str = None) -> 'TaxonomyNode':
        """Add a child node and return it."""
        if name not in self.children:
            self.children[name] = TaxonomyNode(name, self.level + 1, organism_name, lineage)
        return self.children[name]
    
    def add_proteome_data(self, proteome_id: str, protein_count: int):
        """Add proteome data to this node."""
        self.proteome_ids.add(proteome_id)
        self.protein_count += protein_count
        self.total_proteomes = len(self.proteome_ids)
        
        # Track species containing IPR013902 proteins
        if protein_count > 0:
            self.species_with_proteins.add(proteome_id)
    
    def get_species_with_proteins_count(self) -> int:
        """Get the number of species containing IPR013902 proteins."""
        return len(self.species_with_proteins)


def parse_lineage(lineage_str: str) -> List[str]:
    """
    Parse a taxonomic lineage string into a list of taxonomic levels.
    Filters out non-ICNafp conforming rankless names between phylum and class.
    Genus and species come from the Organism column, not from lineage.
    
    Args:
        lineage_str: Comma-separated taxonomic lineage
        
    Returns:
        List of taxonomic names from broad to specific, with non-ICNafp conforming ranks removed.
        Returns empty list if no ICNafp-conforming phylum is found.
    """
    if not lineage_str or lineage_str.strip() == '':
        return []
    
    # Split by comma and clean up each level
    levels = [level.strip() for level in lineage_str.split(',')]
    
    # Remove empty levels, common prefixes, and non-ICNafp conforming ranks
    cleaned_levels = []
    phylum_found = False
    class_found = False
    subphylum_found = False
    
    # Higher-than-phylum ranks that are always valid
    higher_than_phylum_ranks = ['Eukaryota', 'Opisthokonta', 'Fungi']
    
    for level in levels:
        if level and level not in ['cellular organisms']:
            # Check if this is a higher-than-phylum rank
            if level in higher_than_phylum_ranks:
                cleaned_levels.append(level)
                continue
            
            # Check if this is a phylum (ends with 'mycota')
            if level.lower().endswith('mycota'):
                phylum_found = True
                cleaned_levels.append(level)
                continue
            
            # If no phylum found yet, skip all other ranks (they would be below phylum level)
            if not phylum_found:
                continue
            
            # Check if this is a subphylum (ends with 'mycotina')
            if level.lower().endswith('mycotina'):
                subphylum_found = True
                cleaned_levels.append(level)
                continue
            
            # Check if this is a class (ends with 'mycetes')
            if level.lower().endswith('mycetes'):
                class_found = True
                cleaned_levels.append(level)
                continue
            
            # If we're between phylum and class, apply ICNafp filtering
            if phylum_found and not class_found:
                # Special case: if we've found subphylum and this looks like genus/species, 
                # it should be filtered from lineage (will come from Organism column)
                if subphylum_found and (level.endswith(' incertae sedis') or 
                                      (len(level.split()) == 1 and level[0].isupper()) or
                                      (len(level.split()) == 2 and level.split()[0][0].isupper())):
                    # This looks like genus or species in lineage - filter it out
                    continue
                
                # Only filter ranks between phylum and class
                if is_icnafp_conforming_rank(level):
                    cleaned_levels.append(level)
                else:
                    # Skip non-conforming rankless names like "saccharomyceta", "leotiomyceta"
                    pass
            else:
                # After class: allow ICNafp conforming ranks
                if is_icnafp_conforming_rank(level):
                    cleaned_levels.append(level)
                else:
                    # This might be genus/species in lineage (though they should come from Organism column)
                    # Allow them through since filtering should only apply between phylum and class
                    cleaned_levels.append(level)
    
    # If no ICNafp-conforming phylum was found, return empty list
    # This handles cases like Microsporidia which is not an ICNafp-conforming phylum name
    if not phylum_found:
        return []
    
    return cleaned_levels


def build_taxonomy_tree(input_file: str) -> TaxonomyNode:
    """
    Build the taxonomic tree from the input TSV file.
    
    Args:
        input_file: Path to the TSV file with proteome data
        
    Returns:
        Root node of the taxonomic tree
    """
    root = TaxonomyNode("All", 0)
    
    with open(input_file, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        
        for row in reader:
            proteome_id = row['Proteome Id']
            lineage = row['Taxonomic lineage']
            organism_name = row['Organism']  # Extract organism name
            interpro_count = int(row['interpro_IPR013902_count'])
            
            # Parse the taxonomic lineage
            lineage_levels = parse_lineage(lineage)
            
            # Handle proteomes without phylum affiliation
            if not lineage_levels:
                # Create a special category for fungi incertae sedis (uncertain taxonomic placement)
                no_phylum_node = root.add_child("Fungi incertae sedis", None, lineage)
                no_phylum_node.add_proteome_data(proteome_id, interpro_count)
                
                # Extract genus and species from the Organism column
                if organism_name:
                    genus_name = extract_genus_from_organism(organism_name)
                    species_name = extract_species_from_organism(organism_name)
                    
                    # Add genus node under the no-phylum category
                    if genus_name:
                        genus_node = no_phylum_node.add_child(genus_name, None, lineage)
                        genus_node.add_proteome_data(proteome_id, interpro_count)
                        
                        # Add species node under genus
                        if species_name:
                            species_node = genus_node.add_child(species_name, organism_name, lineage)
                            species_node.add_proteome_data(proteome_id, interpro_count)
                continue
            
            # Navigate/build the tree structure for normal lineages
            current_node = root
            current_node.add_proteome_data(proteome_id, interpro_count)
            
            # Process lineage levels but stop after family (taxa ending in 'aceae')
            for i, level_name in enumerate(lineage_levels):
                # Add this level to the tree
                current_node = current_node.add_child(level_name, None, lineage)
                current_node.add_proteome_data(proteome_id, interpro_count)
                
                # If this level is a family (ends with 'aceae'), stop processing lineage levels
                if level_name.lower().endswith('aceae'):
                    break
            
            # Extract genus and species from the Organism column
            if organism_name:
                genus_name = extract_genus_from_organism(organism_name)
                species_name = extract_species_from_organism(organism_name)
                
                # Add genus node
                if genus_name:
                    genus_node = current_node.add_child(genus_name, None, lineage)
                    genus_node.add_proteome_data(proteome_id, interpro_count)
                    
                    # Add species node under genus
                    if species_name:
                        species_node = genus_node.add_child(species_name, organism_name, lineage)
                        species_node.add_proteome_data(proteome_id, interpro_count)
    
    return root


def detect_genus_from_lineage(lineage_str: str) -> str:
    """
    Detect the genus name from taxonomic lineage by finding the taxon immediately 
    after the family (which ends in 'aceae').
    
    Args:
        lineage_str: Comma-separated taxonomic lineage
        
    Returns:
        Genus name if found, empty string otherwise
    """
    if not lineage_str or lineage_str.strip() == '':
        return ''
    
    # Split by comma and clean up each level
    levels = [level.strip() for level in lineage_str.split(',')]
    
    # Find the family (ends with 'aceae') and get the next taxon as genus
    for i, level in enumerate(levels):
        if level.lower().endswith('aceae'):
            # The next taxon after family should be genus
            if i + 1 < len(levels):
                potential_genus = levels[i + 1].strip()
                # Make sure it's not empty and doesn't contain special characters that indicate lower ranks
                if (potential_genus and 
                    not potential_genus.startswith('[') and 
                    '/' not in potential_genus and
                    ' ' not in potential_genus):
                    return potential_genus
    
    return ''


def extract_genus_from_organism(organism_name: str) -> str:
    """
    Extract genus name (first word) from the organism name.
    
    Args:
        organism_name: Full organism name from the Organism column
        
    Returns:
        Genus name (first word)
    """
    if not organism_name or organism_name.strip() == '':
        return ''
    
    # Split by space and take the first word
    words = organism_name.strip().split()
    if words:
        return words[0]
    
    return ''


def extract_species_from_organism(organism_name: str) -> str:
    """
    Extract species name (first two words) from the organism name.
    
    Args:
        organism_name: Full organism name from the Organism column
        
    Returns:
        Species name (genus + species epithet)
    """
    if not organism_name or organism_name.strip() == '':
        return ''
    
    # Split by space and take the first two words
    words = organism_name.strip().split()
    if len(words) >= 2:
        return f"{words[0]} {words[1]}"
    elif len(words) == 1:
        return words[0]  # Just genus if only one word
    
    return ''


def is_icnafp_conforming_rank(taxon_name: str, rank: str = None) -> bool:
    """
    Check if a taxonomic rank conforms to International Code of Nomenclature for algae, fungi, and plants (ICNafp, see Article 37.2) ( `https://www.iaptglobal.org/_functions/code/madrid` ).
    
    According to ICNafp rules for fungi:
    • Phyla should end in: mycota
    • Subphyla should end in: mycotina  
    • Classes should end in: mycetes
    • Subclasses should end in: mycetidae
    • Orders should end in: ales
    • Families should end in: aceae
    
    Rankless names (like "saccharomyceta", "leotiomyceta") that don't follow
    these suffix patterns should be filtered out.
    
    Args:
        taxon_name: Name of the taxonomic unit
        rank: Optional detected rank for additional validation
        
    Returns:
        True if the rank conforms to ICNafp rules, False otherwise
    """
    name_lower = taxon_name.lower()
    
    # Valid ICNafp suffixes for fungi
    valid_suffixes = [
        'mycota',      # phylum
        'mycotina',    # subphylum  
        'mycetes',     # class
        'mycetidae',   # subclass
        'ales',        # order
        'aceae'        # family
    ]
    
    # Check if the name ends with any valid ICNafp suffix
    for suffix in valid_suffixes:
        if name_lower.endswith(suffix):
            return True
    
    # Special cases for genus and species (these don't follow suffix rules)
    if rank in ['genus', 'species', 'variety', 'strain']:
        return True
    
    # Names like "saccharomyceta", "leotiomyceta" don't follow ICNafp rules
    # and should be filtered out
    return False


def detect_taxonomic_rank(node: 'TaxonomyNode') -> str:
    """
    Detect taxonomic rank based on International Code of Nomenclature for algae, fungi, and plants (ICNafp) for fungi.
    For genus detection, uses positional logic: the taxon immediately following a family name.
    
    ICNafp Rules for Fungi:
    • Phyla should end in: mycota
    • Subphyla should end in: mycotina  
    • Classes should end in: mycetes
    • Subclasses should end in: mycetidae
    • Orders should end in: ales
    • Families should end in: aceae
    • Genus: the taxon immediately following the family in the taxonomic lineage
    • Species are sometimes broken into varieties, forms and physiological races
    
    Args:
        node: TaxonomyNode containing the taxon name and lineage information
        
    Returns:
        Detected taxonomic rank as a string
    """
    taxon_name = node.name
    name_lower = taxon_name.lower()
    
    # Check for specific endings based on ICNafp rules
    if name_lower.endswith('mycota'):
        return 'phylum'
    elif name_lower.endswith('mycotina'):
        return 'subphylum'
    elif name_lower.endswith('mycetes'):
        return 'class'
    elif name_lower.endswith('mycetidae'):
        return 'subclass'
    elif name_lower.endswith('ales'):
        return 'order'
    elif name_lower.endswith('aceae'):
        return 'family'
    elif name_lower.endswith('ineae'):
        return 'subfamily'  # Common subfamily ending
    elif name_lower.endswith('eae'):
        return 'tribe'  # Common tribe ending
    else:
        # For genera and species, try to infer based on context
        if ' ' in taxon_name:
            # Contains space, likely species (binomial nomenclature)
            if 'var.' in taxon_name or 'strain' in taxon_name or 'f.' in taxon_name:
                return 'variety/strain'
            else:
                return 'species'
        else:
            # Single word that doesn't match ICNafp conventions
            # Use positional logic for genus detection if lineage is available
            if node.lineage:
                genus_name = detect_genus_from_lineage(node.lineage)
                if genus_name == taxon_name:
                    return 'genus'
            
            # For higher taxonomic groups that don't follow ICNafp naming conventions,
            # we should not display rank information
            known_higher_taxa = ['eukaryota', 'opisthokonta', 'fungi', 'dikarya', 'saccharomyceta']
            if (taxon_name.lower() in known_higher_taxa or 
                len(taxon_name) > 15 or  # Very long names are likely higher taxa
                any(x in taxon_name.lower() for x in ['mycota', 'mycotina', 'mycetes', 'mycetidae'])):
                # Higher taxonomic group that doesn't follow ICNafp conventions
                return ''
            else:
                # Uncertain, don't display rank
                return ''


def format_protein_count(count: int) -> str:
    """
    Format protein count with proper singular/plural grammar.
    
    Args:
        count: Number of proteins
        
    Returns:
        Formatted string with proper grammar
    """
    if count == 1:
        return "1 IPR013902 protein"
    else:
        return f"{count} IPR013902 proteins"


def collect_nodes_by_rank(root: TaxonomyNode) -> Dict[str, List[TaxonomyNode]]:
    """
    Collect all nodes grouped by their taxonomic rank for rank-based alignment.
    
    Args:
        root: Root node of the taxonomic tree
        
    Returns:
        Dictionary mapping rank names to lists of nodes with that rank
    """
    nodes_by_rank = defaultdict(list)
    
    def traverse(node):
        # Skip the root node
        if node.level > 0:
            rank = detect_taxonomic_rank(node)
            if rank:  # Only include nodes with detected ranks
                nodes_by_rank[rank].append(node)
        
        # Traverse children
        for child in node.children.values():
            traverse(child)
    
    traverse(root)
    return dict(nodes_by_rank)





def format_tree_output(node: TaxonomyNode, indent: str = "", is_last: bool = True) -> List[str]:
    """
    Format the taxonomic tree for output in a hierarchical structure.
    
    Args:
        node: Current node to format
        indent: Current indentation string
        is_last: Whether this is the last child at this level
        
    Returns:
        List of formatted lines
    """
    lines = []
    
    # Skip the root node in output
    if node.level > 0:
        # Create the tree connector
        connector = "└─ " if is_last else "├─ "
        
        # Format the node name with protein count and taxonomic rank
        protein_text = format_protein_count(node.protein_count)
        rank = detect_taxonomic_rank(node)
        
        # Add rank information to the protein count text
        if rank and rank in ['phylum', 'subphylum', 'class', 'subclass', 'order', 'family', 'subfamily', 'tribe', 'genus', 'species', 'variety/strain']:
            # Special handling for "Fungi incertae sedis" - remove "in this species" text
            if node.name == "Fungi incertae sedis":
                if node.protein_count == 1:
                    rank_text = "1 IPR013902 protein"
                else:
                    rank_text = f"{node.protein_count} IPR013902 proteins"
            else:
                if node.protein_count == 1:
                    rank_text = f"1 IPR013902 protein in this {rank}"
                else:
                    rank_text = f"{node.protein_count} IPR013902 proteins in this {rank}"
        else:
            # Fallback to original format for unrecognized ranks or empty rank
            rank_text = protein_text
            
        line = f"{indent}{connector}{node.name} ({rank_text})"
        
        lines.append(line)
        
        # Update indent for children
        child_indent = indent + ("   " if is_last else "│  ")
    else:
        # Root node - start with no indent
        child_indent = ""
    
    # Process children
    child_items = list(node.children.items())
    for i, (child_name, child_node) in enumerate(child_items):
        is_last_child = (i == len(child_items) - 1)
        lines.extend(format_tree_output(child_node, child_indent, is_last_child))
    
    return lines


def generate_summary_stats(root: TaxonomyNode) -> Dict[str, int]:
    """
    Generate summary statistics for the taxonomic tree.
    
    Args:
        root: Root node of the taxonomic tree
        
    Returns:
        Dictionary with summary statistics
    """
    stats = {
        'total_proteomes': root.total_proteomes,
        'total_proteins': root.protein_count,
        'proteomes_with_domain': 0,
        'major_groups': len(root.children)
    }
    
    # Count proteomes with at least one IPR013902 protein
    def count_proteomes_with_domain(node):
        proteomes_with_domain = set()
        
        def traverse(n):
            for proteome_id in n.proteome_ids:
                # Check if this proteome has any proteins with the domain
                # by looking at the original data
                pass  # This would require access to original data
            
            for child in n.children.values():
                traverse(child)
        
        traverse(node)
        return len(proteomes_with_domain)
    
    return stats


def generate_phylum_subphylum_summary(root: TaxonomyNode) -> None:
    """
    Generate and print hierarchical summary statistics at phylum and subphylum levels.
    
    Args:
        root: Root node of the taxonomic tree
    """
    print(f"\n📈 Summary:")
    print(f"   Total species analyzed: {root.total_proteomes}")
    print(f"   Total IPR013902 proteins found: {root.protein_count}")
    
    # Build hierarchical phylum-subphylum mapping by traversing the tree
    phylum_subphylum_map = {}
    no_phylum_node = None
    
    # Find Fungi node and check for incertae sedis under it
    def find_fungi_node(node):
        """Recursively find the Fungi node in the tree"""
        if node.name == "Fungi":
            return node
        for child in node.children.values():
            result = find_fungi_node(child)
            if result:
                return result
        return None
    
    fungi_node = find_fungi_node(root)
    if fungi_node and "Fungi incertae sedis" in fungi_node.children:
        no_phylum_node = fungi_node.children["Fungi incertae sedis"]
    
    def traverse_for_hierarchy(node, current_phylum=None):
        """Traverse tree to build phylum-subphylum relationships"""
        rank = detect_taxonomic_rank(node)
        
        # Skip fungi incertae sedis as it's handled separately
        if node.name == "Fungi incertae sedis":
            return
        
        if rank == 'phylum':
            current_phylum = node
            if node.name not in phylum_subphylum_map:
                phylum_subphylum_map[node.name] = {'node': node, 'subphyla': []}
        elif rank == 'subphylum' and current_phylum:
            if current_phylum.name in phylum_subphylum_map:
                phylum_subphylum_map[current_phylum.name]['subphyla'].append(node)
        
        # Continue traversing children
        for child in node.children.values():
            traverse_for_hierarchy(child, current_phylum)
    
    traverse_for_hierarchy(root)
    
    # Display hierarchical distribution
    if phylum_subphylum_map or no_phylum_node:
        print(f"\n🔬 Taxonomic distribution (Phylum → Subphylum):")
        
        # Sort phyla by protein count (descending)
        sorted_phyla = sorted(phylum_subphylum_map.items(), 
                            key=lambda x: x[1]['node'].protein_count, reverse=True)
        
        for phylum_name, phylum_data in sorted_phyla:
             phylum_node = phylum_data['node']
             subphyla = phylum_data['subphyla']
             
             species_with_proteins = phylum_node.get_species_with_proteins_count()
             if species_with_proteins > 0:
                 protein_word = "protein" if phylum_node.protein_count == 1 else "proteins"
                 print(f"   {phylum_name} ({phylum_node.total_proteomes} species): {phylum_node.protein_count} IPR013902 {protein_word} in {species_with_proteins} species")
             else:
                 protein_word = "protein" if phylum_node.protein_count == 1 else "proteins"
                 print(f"   {phylum_name} ({phylum_node.total_proteomes} species): {phylum_node.protein_count} IPR013902 {protein_word}")
             
             # Sort subphyla by protein count (descending)
             subphyla.sort(key=lambda x: x.protein_count, reverse=True)
             
             # Display subphyla with indentation
             for subphylum in subphyla:
                 subphylum_species_with_proteins = subphylum.get_species_with_proteins_count()
                 if subphylum_species_with_proteins > 0:
                     protein_word = "protein" if subphylum.protein_count == 1 else "proteins"
                     print(f"      └─ {subphylum.name} ({subphylum.total_proteomes} species): {subphylum.protein_count} IPR013902 {protein_word} in {subphylum_species_with_proteins} species")
                 else:
                     protein_word = "protein" if subphylum.protein_count == 1 else "proteins"
                     print(f"      └─ {subphylum.name} ({subphylum.total_proteomes} species): {subphylum.protein_count} IPR013902 {protein_word}")
             
             # Add spacing between phyla for readability
             print()
        
        # Display fungi incertae sedis at the end
        if no_phylum_node:
            incertae_sedis_species_with_proteins = no_phylum_node.get_species_with_proteins_count()
            if incertae_sedis_species_with_proteins > 0:
                protein_word = "protein" if no_phylum_node.protein_count == 1 else "proteins"
                print(f"   Fungi incertae sedis (fungi with uncertain taxonomic placement) ({no_phylum_node.total_proteomes} species): {no_phylum_node.protein_count} IPR013902 {protein_word} in {incertae_sedis_species_with_proteins} species")
            else:
                protein_word = "protein" if no_phylum_node.protein_count == 1 else "proteins"
                print(f"   Fungi incertae sedis (fungi with uncertain taxonomic placement) ({no_phylum_node.total_proteomes} species): {no_phylum_node.protein_count} IPR013902 {protein_word}")
            
            # Show genera under this category
            genera = [child for child in no_phylum_node.children.values()]
            if genera:
                genera.sort(key=lambda x: x.protein_count, reverse=True)
                for genus in genera:
                    genus_species_with_proteins = genus.get_species_with_proteins_count()
                    if genus_species_with_proteins > 0:
                        protein_word = "protein" if genus.protein_count == 1 else "proteins"
                        print(f"      └─ {genus.name} ({genus.total_proteomes} species): {genus.protein_count} IPR013902 {protein_word} in {genus_species_with_proteins} species")
                    else:
                        protein_word = "protein" if genus.protein_count == 1 else "proteins"
                        print(f"      └─ {genus.name} ({genus.total_proteomes} species): {genus.protein_count} IPR013902 {protein_word}")


def main():
    """
    Main function to process the taxonomic data and generate the tree output.
    """
    # Set up command line argument parsing
    parser = argparse.ArgumentParser(
        description='Generate a hierarchical tree showing the distribution of IPR013902-containing proteins across fungal taxonomic levels.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python taxonomy_tree.py                                    # Use default input file
  python taxonomy_tree.py data.tsv                          # Specify input file
  python taxonomy_tree.py data.tsv output.txt               # Specify input and output files
  python taxonomy_tree.py -h                                # Show this help message

Input file format:
  The input TSV file should contain columns for:
  - Proteome ID
  - Organism name
  - Taxonomic lineage
  - IPR013902 protein count

Data sources analyzed:
  - Taxonomic lineage information for hierarchical classification
  - Organism names for genus/species extraction
  - IPR013902 protein counts per proteome
  - Proteome identifiers for species tracking
        """
    )
    
    parser.add_argument(
        'input_file',
        nargs='?',
        default='fungal_ref_proteomes_with_interpro.tsv',
        help='Path to the TSV file containing fungal proteome data with IPR013902 counts (default: fungal_ref_proteomes_with_interpro.tsv)'
    )
    
    parser.add_argument(
        'output_file',
        nargs='?',
        default=None,
        help='Optional output file path. If not provided, prints to stdout'
    )
    
    # Parse arguments
    args = parser.parse_args()
    input_file = args.input_file
    output_file = args.output_file
    
    try:
        print(f"📊 Processing taxonomic data from: {input_file}")
        
        # Build the taxonomic tree
        root = build_taxonomy_tree(input_file)
        
        # Find the Fungi node in the tree structure
        def find_fungi_node(node):
            """Recursively find the Fungi node in the tree"""
            if node.name == "Fungi":
                return node
            for child in node.children.values():
                result = find_fungi_node(child)
                if result:
                    return result
            return None
        
        fungi_node = find_fungi_node(root)
        
        # Find and move "Fungi incertae sedis" node from root to under Fungi
        incertae_sedis_node = None
        if "Fungi incertae sedis" in root.children:
            incertae_sedis_node = root.children["Fungi incertae sedis"]
            # Store the counts before removing from root
            incertae_sedis_protein_count = incertae_sedis_node.protein_count
            incertae_sedis_proteome_count = incertae_sedis_node.total_proteomes
            incertae_sedis_species_with_proteins = incertae_sedis_node.species_with_proteins.copy()
            
            # Remove from root
            del root.children["Fungi incertae sedis"]
            
            # Add to Fungi node if it exists
            if fungi_node:
                fungi_node.children["Fungi incertae sedis"] = incertae_sedis_node
                # Update protein count in Fungi node
                fungi_node.protein_count += incertae_sedis_protein_count
                fungi_node.total_proteomes += incertae_sedis_proteome_count
                # Update species with proteins in Fungi node
                fungi_node.species_with_proteins.update(incertae_sedis_species_with_proteins)
                
                # Ensure root node maintains correct total counts for summary
                # (The removal above would have reduced root counts, so we need to add them back)
                root.protein_count += incertae_sedis_protein_count
                root.total_proteomes += incertae_sedis_proteome_count
                root.species_with_proteins.update(incertae_sedis_species_with_proteins)
        
        # Generate the tree output starting from Fungi node
        tree_lines = []
        if fungi_node:
            # Start with Fungi as the root with a horizontal line
            tree_lines.append(f"─ {fungi_node.name} ({fungi_node.protein_count} IPR013902 proteins)")
            
            # Process children of Fungi node with proper alignment
            # The "F" in "Fungi" is at position 2 (after "─ "), so children should align there
            child_items = list(fungi_node.children.items())
            
            # Separate regular phyla from incertae sedis
            regular_children = [(name, node) for name, node in child_items if name != "Fungi incertae sedis"]
            incertae_sedis_child = [(name, node) for name, node in child_items if name == "Fungi incertae sedis"]
            
            # Process regular children first
            for i, (child_name, child_node) in enumerate(regular_children):
                is_last_child = (i == len(regular_children) - 1) and not incertae_sedis_child
                tree_lines.extend(format_tree_output(child_node, "  ", is_last_child))
            
            # Process incertae sedis at the bottom
            for i, (child_name, child_node) in enumerate(incertae_sedis_child):
                is_last_child = True  # Always last since it's at the bottom
                tree_lines.extend(format_tree_output(child_node, "  ", is_last_child))
        else:
            # Fallback to original behavior if Fungi node not found
            child_items = list(root.children.items())
            for i, (child_name, child_node) in enumerate(child_items):
                is_last_child = (i == len(child_items) - 1)
                tree_lines.extend(format_tree_output(child_node, "", is_last_child))
        
        # Add header with blank line
        output_lines = [
            "Taxonomy tree illustrating the distribution IPR013902 domain proteins:",
            ""  # Add blank line after title
        ]
        
        # Add the rank-aligned tree structure
        output_lines.extend(tree_lines)
        
        # Output results
        if output_file:
            with open(output_file, 'w', encoding='utf-8') as f:
                for line in output_lines:
                    f.write(line + '\n')
            print(f"✅ Taxonomic tree saved to: {output_file}")
        else:
            for line in output_lines:
                print(line)
        
        # Print phylum and subphylum level summary statistics
        generate_phylum_subphylum_summary(root)
        
    except FileNotFoundError:
        print(f"❌ Error: Input file '{input_file}' not found.")
        sys.exit(1)
    except Exception as e:
        print(f"❌ Error processing data: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()