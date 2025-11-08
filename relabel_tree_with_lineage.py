#!/usr/bin/env python3
"""
Relabel Newick tree leaves by replacing Taxon mnemonic codes with lineage-based labels.

This script reads a Newick tree file and a TSV file (fungal reference proteomes),
builds a mapping from the "Taxon mnemonic" (e.g., DEBHA) to a formatted lineage label
constructed from the "Taxonomic lineage" and the species epithet from "Organism".

Example transformation:
  Leaf label containing "DEBHA" -> "Fungi_Dikarya_Ascomycota_..._Debaryomycetaceae_hansenii"

CLI usage provides help messages and options to customize behavior.
"""

import argparse
import re
import sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional


def read_tsv_build_mapping(
    tsv_path: Path,
    mnemonic_col: str = "Taxon mnemonic",
    lineage_col: str = "Taxonomic lineage",
    organism_col: str = "Organism",
    start_rank: Optional[str] = None,
) -> Dict[str, str]:
    """
    Build a mapping from Taxon mnemonic to formatted lineage label.

    - Reads a tab-separated file containing at least the columns:
      mnemonic_col, lineage_col, organism_col.
    - Parses the lineage into ranks. If `start_rank` is provided, lineage will be
      sliced starting from that rank; otherwise all ranks are used.
      - Appends the lowercased species epithet (from organism) at the end.

    Args:
        tsv_path: Path to the TSV file.
        mnemonic_col: Column name for Taxon mnemonic codes.
        lineage_col: Column name for Taxonomic lineage.
        organism_col: Column name for Organism name.
        start_rank: Rank to start lineage from (case-insensitive match).

    Returns:
        Dict mapping mnemonic (e.g., DEBHA) -> formatted lineage label.
    """
    mapping: Dict[str, str] = {}
    if not tsv_path.exists():
        raise FileNotFoundError(f"TSV file not found: {tsv_path}")

    with tsv_path.open("r", encoding="utf-8") as fh:
        header_line = fh.readline()
        if not header_line:
            raise ValueError("TSV file is empty or missing header")
        headers = header_line.rstrip("\n").split("\t")
        try:
            idx_mn = headers.index(mnemonic_col)
            idx_lin = headers.index(lineage_col)
            idx_org = headers.index(organism_col)
        except ValueError as e:
            raise ValueError(
                f"Required column not found in TSV header: {e}. Headers: {headers}"
            )

        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= max(idx_mn, idx_lin, idx_org):
                # Skip malformed lines
                continue
            mnemonic = parts[idx_mn].strip()
            lineage = parts[idx_lin].strip()
            organism = parts[idx_org].strip()

            if not mnemonic:
                continue
            formatted = format_lineage_label(lineage, organism, start_rank=start_rank)
            if formatted:
                mapping[mnemonic] = formatted

    return mapping


def sanitize_rank_name(name: str) -> str:
    """
    Sanitize a taxonomic rank name for use in Newick labels.

    - Removes parentheses and their contents.
    - Replaces spaces with underscores.
    - Removes disallowed punctuation except underscore and hyphen.

    Args:
        name: Raw rank name.

    Returns:
        Sanitized rank name suitable for Newick labels.
    """
    # Remove parenthetical notes
    name = re.sub(r"\(.*?\)", "", name)
    # Replace spaces with underscores
    name = name.replace(" ", "_")
    # Keep alphanumerics, underscore, hyphen; replace others with underscore
    name = re.sub(r"[^A-Za-z0-9_\-]", "_", name)
    # Collapse multiple underscores
    name = re.sub(r"_+", "_", name).strip("_")
    return name


def extract_species_epithet(organism: str) -> str:
    """
    Extract the species epithet (second token) from an organism string.

    Example: "Debaryomyces hansenii (strain CBS 767)" -> "hansenii"

    Args:
        organism: Organism string from TSV.

    Returns:
        Lowercased species epithet if available, else empty string.
    """
    # Remove parentheses content
    base = re.sub(r"\(.*?\)", "", organism).strip()
    tokens = base.split()
    if len(tokens) >= 2:
        return tokens[1].lower()
    return ""


def format_lineage_label(lineage: str, organism: str, start_rank: Optional[str] = None) -> str:
    """
    Convert a taxonomic lineage string and organism to a formatted label.

    - Splits lineage by semicolons (or commas if present).
    - Starts from the specified `start_rank` (case-insensitive), if found.
    - Sanitizes each rank and joins with underscores.
    - Appends the species epithet token at the end if available.

    Args:
        lineage: Taxonomic lineage field.
        organism: Organism field for species epithet.
        start_rank: Rank to start from (e.g., "Fungi").

    Returns:
        Formatted label string.
    """
    # Split lineage by common delimiters
    ranks = [r.strip() for r in re.split(r"[;]\s*", lineage) if r.strip()]
    if not ranks:
        return ""

    # Optionally slice lineage starting at the requested rank
    if start_rank:
        sanitized_start = sanitize_rank_name(start_rank)
        start_idx = 0
        for i, r in enumerate(ranks):
            r_san = sanitize_rank_name(r)
            if (
                r.lower() == start_rank.lower()
                or r_san.lower() == sanitized_start.lower()
                or start_rank.lower() in r.lower()
            ):
                start_idx = i
                break
        ranks = ranks[start_idx:]

    # Extract genus and species epithet from organism
    base_org = re.sub(r"\(.*?\)", "", organism).strip()
    org_tokens = base_org.split()
    genus = org_tokens[0] if org_tokens else ""
    epithet = org_tokens[1].lower() if len(org_tokens) >= 2 else ""

    # Filter out genus/species ranks from lineage to avoid duplication.
    # Explanation: Many lineage strings already include genus and sometimes species
    # (e.g., "Debaryomyces hansenii"). Since we append the species epithet at the
    # end of the label for consistency, we strip any lineage components that match
    # the genus token or the full "genus species" string, and also avoid appending
    # a standalone species epithet found earlier in the lineage. This prevents
    # duplicated genus/species information in the final label.
    filtered_ranks: List[str] = []
    for r in ranks:
        r_clean = r.strip()
        r_lower = r_clean.lower()
        genus_lower = genus.lower()
        species_full_lower = f"{genus_lower} {epithet}" if genus and epithet else ""

        # Skip if rank includes genus or species info
        if genus and (r_lower == genus_lower or genus_lower in r_lower):
            continue
        if species_full_lower and species_full_lower in r_lower:
            continue
        if epithet and epithet in r_lower:
            continue
        filtered_ranks.append(r_clean)

    # Sanitize ranks
    sanitized = [sanitize_rank_name(r) for r in filtered_ranks if r]

    # Append species epithet at the end
    if epithet:
        sanitized.append(epithet)

    # Join with underscores
    label = "_".join([s for s in sanitized if s])
    return label


def relabel_newick_leaves(newick_text: str, mapping: Dict[str, str]) -> str:
    """
    Relabel leaf names in a Newick string by replacing taxon mnemonic occurrences.

    - Detect leaf tokens (names immediately following '(' or ',' and before ':'/','/')')
    - If a mnemonic key from mapping is found within the leaf token, replace that
      mnemonic substring with the formatted lineage label, preserving any suffix
      (e.g., "/205-277").

    Args:
        newick_text: Original Newick content as a single-line string.
        mapping: Dict of mnemonic -> lineage label.

    Returns:
        Modified Newick string with relabeled leaves.
    """
    out_chars: List[str] = []
    i = 0
    n = len(newick_text)
    expecting_label = False

    # Precompute keys sorted by length (longest first) to avoid partial matches
    keys_sorted = sorted(mapping.keys(), key=len, reverse=True)

    while i < n:
        ch = newick_text[i]
        out_chars.append(ch)

        if ch == '(' or ch == ',':
            # Next token could be a leaf name
            expecting_label = True
            i += 1
            # Peek ahead: if next is '(' or ')' or ',' treat accordingly
            if i < n and newick_text[i] in '(),':
                expecting_label = False
            else:
                # Capture the token until a delimiter ':', ',', or ')'
                start = i
                while i < n and newick_text[i] not in ':,)':
                    i += 1
                token = newick_text[start:i]

                # Attempt replacement if any mnemonic key is present
                replaced = token
                for key in keys_sorted:
                    if key and key in token:
                        replaced = token.replace(key, mapping[key])
                        break

                out_chars.append(replaced)

                expecting_label = False
                # Do not advance i here; the loop will handle delimiter at position i
                continue
        else:
            expecting_label = False
            i += 1

    return "".join(out_chars)


def main() -> None:
    """
    Command-line interface to relabel Newick leaves using lineage mapping from TSV.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Relabel Newick leaves by replacing Taxon mnemonic codes with "
            "formatted lineage labels derived from a TSV file."
        )
    )
    parser.add_argument(
        "--input-tree",
        required=True,
        help="Path to input Newick tree file",
    )
    parser.add_argument(
        "--output-tree",
        required=True,
        help="Path to write the relabeled Newick tree",
    )
    parser.add_argument(
        "--tsv",
        default=str(Path("/Users/lilindu/Documents/Du_Lab_papers_grants/cryo-tdk-killers/fungal_ref_proteomes_with_interpro.tsv")),
        help="Path to TSV file containing Taxon mnemonic, Taxonomic lineage, and Organism columns",
    )
    parser.add_argument(
        "--mnemonic-column",
        default="Taxon mnemonic",
        help="Column name for Taxon mnemonic (default: 'Taxon mnemonic')",
    )
    parser.add_argument(
        "--lineage-column",
        default="Taxonomic lineage",
        help="Column name for Taxonomic lineage (default: 'Taxonomic lineage')",
    )
    parser.add_argument(
        "--organism-column",
        default="Organism",
        help="Column name for Organism (default: 'Organism')",
    )
    parser.add_argument(
        "--start-rank",
        default="",
        help=(
            "Optional: start lineage from this rank (case-insensitive). "
            "Leave empty to use the full lineage."
        ),
    )

    args = parser.parse_args()

    input_tree = Path(args.input_tree)
    output_tree = Path(args.output_tree)
    tsv_path = Path(args.tsv)

    # Build mapping from mnemonic to lineage label
    mapping = read_tsv_build_mapping(
        tsv_path,
        mnemonic_col=args.mnemonic_column,
        lineage_col=args.lineage_column,
        organism_col=args.organism_column,
        start_rank=(args.start_rank if args.start_rank.strip() else None),
    )

    if not mapping:
        print("Warning: No mapping entries constructed from TSV; output may be unchanged.")

    # Read input Newick
    try:
        newick_text = input_tree.read_text(encoding="utf-8").strip()
    except FileNotFoundError:
        print(f"Error: Input tree file not found: {input_tree}")
        sys.exit(1)

    # Relabel leaves
    relabeled = relabel_newick_leaves(newick_text, mapping)

    # Write output
    output_tree.write_text(relabeled + ("\n" if not relabeled.endswith("\n") else ""), encoding="utf-8")
    print(f"Relabeled tree written to: {output_tree}")


if __name__ == "__main__":
    main()