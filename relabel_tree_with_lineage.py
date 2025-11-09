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

    - Splits lineage by semicolons or commas.
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
    # Split lineage by common delimiters (semicolon or comma)
    ranks = [r.strip() for r in re.split(r"[;,]\s*", lineage) if r.strip()]
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


def relabel_newick_leaves_with_records(newick_text: str, mapping: Dict[str, str]) -> Tuple[str, List[Dict[str, str]]]:
    """
    Relabel leaf names in a Newick string and collect per-leaf reporting records.

    This function performs the same tokenization and replacement as
    `relabel_newick_leaves`, but additionally collects a record for each leaf
    containing:
    - original leaf node label
    - whether the original leaf node label has mnemonic (yes/no)
    - current leaf node label (after potential relabeling)
    - lineage information in the current leaf node label (mapping value used)

    Args:
        newick_text: Original Newick content as a single-line string.
        mapping: Dict of mnemonic -> lineage label.

    Returns:
        A tuple of (modified Newick string, list of records for TSV reporting).
    """
    out_chars: List[str] = []
    i = 0
    n = len(newick_text)
    expecting_label = False
    records: List[Dict[str, str]] = []

    # Precompute keys sorted by length (longest first) to avoid partial matches
    keys_sorted = sorted(mapping.keys(), key=len, reverse=True)

    while i < n:
        ch = newick_text[i]
        out_chars.append(ch)

        if ch == '(' or ch == ',':
            expecting_label = True
            i += 1
            if i < n and newick_text[i] in '(),':
                expecting_label = False
            else:
                start = i
                while i < n and newick_text[i] not in ':,)':
                    i += 1
                token = newick_text[start:i]

                # Determine if token includes any mnemonic and perform replacement
                has_mnemonic = False
                used_lineage = ""
                replaced = token
                for key in keys_sorted:
                    if key and key in token:
                        has_mnemonic = True
                        # Replace mnemonic with lineage label
                        replaced = token.replace(key, mapping.get(key, key))
                        used_lineage = mapping.get(key, "")
                        break

                out_chars.append(replaced)

                # Collect reporting record
                records.append({
                    "original leaf node label": token,
                    "whether the original leaf node label has mnemonic": "yes" if has_mnemonic else "no",
                    "current leaf node label": replaced,
                    "lineage information in the current leaf node label": used_lineage,
                })

                expecting_label = False
                continue
        else:
            expecting_label = False
            i += 1

    return "".join(out_chars), records


def write_report_tsv(report_path: Path, records: List[Dict[str, str]]) -> None:
    """
    Write a TSV report of leaf relabeling.

    Columns written:
    - original leaf node label
    - whether the original leaf node label has mnemonic
    - current leaf node label
    - lineage information in the current leaf node label

    Args:
        report_path: Path to write the TSV file.
        records: List of record dicts as produced by relabel_newick_leaves_with_records.
    """
    headers = [
        "original leaf node label",
        "whether the original leaf node label has mnemonic",
        "current leaf node label",
        "lineage information in the current leaf node label",
    ]
    lines = ["\t".join(headers)]
    for r in records:
        row = [r.get(h, "") for h in headers]
        lines.append("\t".join(row))
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """
    Command-line interface to relabel Newick leaves using lineage mapping from TSV.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Relabel Newick leaves by replacing Taxon mnemonic codes with "
            "formatted lineage labels derived from a TSV file."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
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
        "--report-tsv",
        default="",
        help=(
            "Optional: path to write TSV report listing per-leaf original/current labels "
            "and lineage information. If omitted, a default path derived from --output-tree "
            "will be used (<output_tree_stem>.relabel_report.tsv)."
        ),
    )
    parser.add_argument(
        "--tsv",
        default="fungal_ref_proteomes_with_interpro.tsv",
        help=(
            "Path to TSV file containing Taxon mnemonic, Taxonomic lineage, and Organism columns. "
            "If not provided, defaults to 'fungal_ref_proteomes_with_interpro.tsv' in the current working directory."
        ),
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

    # Relabel leaves and collect per-leaf records
    relabeled, records = relabel_newick_leaves_with_records(newick_text, mapping)

    # Compute summary statistics
    total_leaves = len(records)
    with_mnemonic = sum(1 for r in records if r["whether the original leaf node label has mnemonic"] == "yes")
    without_mnemonic = total_leaves - with_mnemonic
    successfully_relabeled = sum(
        1 for r in records
        if r["whether the original leaf node label has mnemonic"] == "yes"
        and r["original leaf node label"] != r["current leaf node label"]
        and r["lineage information in the current leaf node label"] != ""
    )
    not_successfully_relabeled = with_mnemonic - successfully_relabeled

    # Determine report TSV path and write report
    if args.report_tsv and args.report_tsv.strip():
        report_path = Path(args.report_tsv)
    else:
        report_path = output_tree.parent / f"{output_tree.stem}.relabel_report.tsv"

    write_report_tsv(report_path, records)

    # Write output tree
    output_tree.write_text(relabeled + ("\n" if not relabeled.endswith("\n") else ""), encoding="utf-8")

    # Print summary to command line
    print("Relabeling summary:")
    print(f"- Total leaf nodes: {total_leaves}")
    print(f"- Leaf nodes without mnemonic: {without_mnemonic}")
    print(f"- Leaf nodes with mnemonic: {with_mnemonic}")
    print(f"  - Successfully relabeled with lineage: {successfully_relabeled}")
    print(f"  - Not successfully relabeled: {not_successfully_relabeled}")
    print(f"TSV report written to: {report_path}")


if __name__ == "__main__":
    main()