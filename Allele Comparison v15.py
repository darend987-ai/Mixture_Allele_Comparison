import os
import glob
import traceback
from datetime import datetime
from decimal import Decimal, InvalidOperation

import pandas as pd

# -------------------------
# Normalization & helpers
# -------------------------
def canonical_allele(val):
    """Normalize allele values to canonical strings (10.0 -> '10', 9.30 -> '9.3')."""
    if pd.isna(val):
        return None
    s = str(val).strip()
    if s == "" or s.lower() == "nan":
        return None
    try:
        d = Decimal(s)
    except (InvalidOperation, ValueError):
        return s
    if d == d.to_integral():
        return str(int(d))
    s_norm = format(d.normalize(), 'f')
    if s_norm.endswith('.'):
        s_norm = s_norm[:-1]
    return s_norm

def allele_sort_key(a):
    """Sort alleles numerically when possible, otherwise lexicographically."""
    try:
        d = Decimal(a)
        return (0, float(d), "")
    except Exception:
        return (1, float('inf'), str(a))

def sanitize_sheet_name(name, max_len=31):
    """Make a valid Excel sheet name: trim and replace problematic chars."""
    invalid_chars = ['\\','/','?','*',':','[',']']
    for ch in invalid_chars:
        name = name.replace(ch, '_')
    name = name.strip()
    if len(name) > max_len:
        name = name[:max_len-3] + '...'
    return name

# -------------------------
# File detection
# -------------------------
def find_most_recent(files):
    return max(files, key=os.path.getmtime)

def find_file_by_keywords(files, keywords):
    matches = [f for f in files if any(k in f.lower() for k in keywords)]
    if not matches:
        return None
    if len(matches) == 1:
        return matches[0]
    # pick most recently modified
    return find_most_recent(matches)

def find_all_evidence(files, keywords):
    matches = [f for f in files if any(k in f.lower() for k in keywords)]
    # sort by filename for deterministic order, but preserve modified-time for tie-breaking if needed
    matches = sorted(matches, key=lambda x: (os.path.getmtime(x), x))
    return matches

# -------------------------
# CSV reading helpers
# -------------------------
def read_wide_profile(csv_path):
    """Read a wide-format CSV and return dict: locus -> set of canonical alleles."""
    df = pd.read_csv(csv_path, dtype=str)
    if df.shape[1] < 1:
        return {}
    locus_col = df.columns[0]
    profile = {}
    for _, row in df.iterrows():
        locus = str(row[locus_col]).strip()
        alleles = []
        for val in row[1:]:
            a = canonical_allele(val)
            if a is not None:
                alleles.append(a)
        profile[locus] = set(alleles)
    return profile

# -------------------------
# Per-evidence comparison
# -------------------------
def compare_profiles(contrib_profile, suspect_profile, evidence_profile):
    """Return per-locus DataFrame and summary DataFrame for one evidence profile."""
    all_loci = sorted(set(contrib_profile.keys()) | set(suspect_profile.keys()) | set(evidence_profile.keys()))
    rows = []
    summary_counts = {
        "Obligate_of_Suspect": 0,
        "Obligate_of_Assumed_Contributor": 0,
        "Shared_All_Three": 0,
        "Missing_from_Suspect": 0,
        "Missing_from_Assumed_Contributor": 0,
        "Foreign_Alleles": 0
    }

    for locus in all_loci:
        c = contrib_profile.get(locus, set())
        s = suspect_profile.get(locus, set())
        e = evidence_profile.get(locus, set())

        must_from_suspect = sorted(list((e & s) - c), key=allele_sort_key)
        must_from_contrib = sorted(list((e & c) - s), key=allele_sort_key)
        shared_all_three = sorted(list(e & s & c), key=allele_sort_key)
        missing_from_evidence_suspect = sorted(list(s - e), key=allele_sort_key)
        missing_from_evidence_contrib = sorted(list(c - e), key=allele_sort_key)
        foreign_alleles = sorted(list(e - (s | c)), key=allele_sort_key)

        rows.append({
            "Locus": locus,
            "Evidence Alleles": ", ".join(sorted(e, key=allele_sort_key)),
            "Suspect Alleles": ", ".join(sorted(s, key=allele_sort_key)),
            "Contributor Alleles": ", ".join(sorted(c, key=allele_sort_key)),
            "Obligate_of_Suspect": ", ".join(must_from_suspect) if must_from_suspect else "-",
            "Obligate_of_Assumed_Contributor": ", ".join(must_from_contrib) if must_from_contrib else "-",
            "Shared_All_Three": ", ".join(shared_all_three) if shared_all_three else "-",
            "Missing_from_Suspect": ", ".join(missing_from_evidence_suspect) if missing_from_evidence_suspect else "-",
            "Missing_from_Assumed_Contributor": ", ".join(missing_from_evidence_contrib) if missing_from_evidence_contrib else "-",
            "Foreign_Alleles": ", ".join(foreign_alleles) if foreign_alleles else "-",
            "Obligate_of_Suspect_Count": len(must_from_suspect),
            "Obligate_of_Assumed_Contributor_Count": len(must_from_contrib),
            "Shared_All_Three_Count": len(shared_all_three),
            "Missing_from_Suspect_Count": len(missing_from_evidence_suspect),
            "Missing_from_Assumed_Contributor_Count": len(missing_from_evidence_contrib),
            "Foreign_Alleles_Count": len(foreign_alleles)
        })

        # increment summary totals
        summary_counts["Obligate_of_Suspect"] += len(must_from_suspect)
        summary_counts["Obligate_of_Assumed_Contributor"] += len(must_from_contrib)
        summary_counts["Shared_All_Three"] += len(shared_all_three)
        summary_counts["Missing_from_Suspect"] += len(missing_from_evidence_suspect)
        summary_counts["Missing_from_Assumed_Contributor"] += len(missing_from_evidence_contrib)
        summary_counts["Foreign_Alleles"] += len(foreign_alleles)

    df_per_locus = pd.DataFrame(rows)
    df_summary = pd.DataFrame({
        "Category": list(summary_counts.keys()),
        "Count": list(summary_counts.values())
    })
    return df_per_locus, df_summary, summary_counts

# -------------------------
# Top-level runner
# -------------------------
def run_multi_evidence(output_basename=None):
    try:
        cwd = os.getcwd()
        csv_files = [os.path.basename(f) for f in glob.glob(os.path.join(cwd, "*.csv"))]
        if not csv_files:
            raise RuntimeError("No CSV files found in current directory.")

        # keywords (case-insensitive)
        suspect_keywords = ["suspect", "sus"]
        assumed_keywords = ["assumed", "contrib", "contributor", "donor"]
        evidence_keywords = ["evidence", "evid", "sample", "item"]

        suspect_file = find_file_by_keywords(csv_files, suspect_keywords)
        assumed_file = find_file_by_keywords(csv_files, assumed_keywords)
        evidence_files = find_all_evidence(csv_files, evidence_keywords)

        # fallback: if no "assumed" match, try "contrib"/"contributor" specifically
        if not assumed_file:
            assumed_file = find_file_by_keywords(csv_files, ["contrib", "contributor", "donor"])

        if not suspect_file or not assumed_file or not evidence_files:
            raise RuntimeError(
                f"Couldn't auto-detect required files.\nDetected CSVs: {csv_files}\n"
                f"Suspect: {suspect_file}\nAssumed contributor: {assumed_file}\nEvidence files: {evidence_files}\n"
                "Ensure you have one suspect file (name contains 'suspect'), one assumed contributor (name contains 'assumed' or 'contrib'), and at least one evidence file (name contains 'evidence' or 'sample')."
            )

        print("Detected files:")
        print("  Suspect:", suspect_file)
        print("  Assumed contributor:", assumed_file)
        print("  Evidence files:", ", ".join(evidence_files))

        # read suspect & assumed once
        suspect_profile = read_wide_profile(suspect_file)
        assumed_profile = read_wide_profile(assumed_file)

        # prepare output filename
        timestamp = datetime.now().strftime("%Y-%m-%d_%H%M%S")
        if output_basename:
            out_name = f"{output_basename}_{timestamp}.xlsx"
        else:
            out_name = f"Allele_Comparison_Results_{timestamp}.xlsx"
        out_path = os.path.abspath(out_name)

        # track master totals across all evidence files
        master_totals = {
            "Obligate_of_Suspect": 0,
            "Obligate_of_Assumed_Contributor": 0,
            "Shared_All_Three": 0,
            "Missing_from_Suspect": 0,
            "Missing_from_Assumed_Contributor": 0,
            "Foreign_Alleles": 0
        }

        # Write all sheets into one workbook
        with pd.ExcelWriter(out_path, engine="openpyxl") as writer:
            for idx, ev_file in enumerate(evidence_files, start=1):
                try:
                    ev_profile = read_wide_profile(ev_file)
                except Exception as e:
                    # write a small sheet to indicate read error for this evidence file
                    err_df = pd.DataFrame({"Error": [f"Failed to read {ev_file}: {e}"]})
                    sheet = sanitize_sheet_name(f"Error_{os.path.splitext(ev_file)[0]}")
                    err_df.to_excel(writer, sheet_name=sheet, index=False)
                    continue

                df_per_locus, df_summary, counts = compare_profiles(assumed_profile, suspect_profile, ev_profile)

                # update master totals
                for k, v in counts.items():
                    master_totals[k] += v

                base_name = os.path.splitext(ev_file)[0]
                per_locus_sheet = sanitize_sheet_name(f"{base_name}_per_locus")
                summary_sheet = sanitize_sheet_name(f"{base_name}_summary")

                # ensure sheet name uniqueness by appending index if needed
                # pandas/openpyxl will raise if duplicate sheet names; handle proactively
                existing = writer.book.sheetnames if hasattr(writer, "book") else []
                if per_locus_sheet in existing:
                    per_locus_sheet = sanitize_sheet_name(f"{base_name}_per_locus_{idx}")
                if summary_sheet in existing:
                    summary_sheet = sanitize_sheet_name(f"{base_name}_summary_{idx}")

                df_per_locus.to_excel(writer, sheet_name=per_locus_sheet, index=False)
                df_summary.to_excel(writer, sheet_name=summary_sheet, index=False)

            # after processing all evidence files, write master summary
            master_df = pd.DataFrame({
                "Category": list(master_totals.keys()),
                "Total_Count": list(master_totals.values())
            })
            master_sheet = sanitize_sheet_name("Master_Summary")
            master_df.to_excel(writer, sheet_name=master_sheet, index=False)

        print(f"\n✅ All done — results written to: {out_path}")
        return out_path

    except Exception as ex:
        tb = traceback.format_exc()
        print("\n❌ Fatal error:\n", tb)
        # attempt to write an error workbook so user has artifact
        ts = datetime.now().strftime("%Y-%m-%d_%H%M%S")
        err_name = os.path.abspath(f"allele_report_error_{ts}.xlsx")
        try:
            err_df = pd.DataFrame({"Error": [str(ex)], "Traceback": [tb]})
            with pd.ExcelWriter(err_name, engine="openpyxl") as writer:
                err_df.to_excel(writer, sheet_name="Error_Log", index=False)
            print(f"⚠️ Wrote error workbook: {err_name}")
        except Exception as e2:
            print("Also failed to write error workbook:", e2)
        raise

# -------------------------
# If executed directly
# -------------------------
if __name__ == "__main__":
    run_multi_evidence()
