import argparse
from pathlib import Path
from typing import List, Tuple
import sys

import numpy as np
import pandas as pd

# Ensure imports work when launched from workspace root.
THIS_DIR = Path(__file__).resolve().parent
if str(THIS_DIR) not in sys.path:
    sys.path.insert(0, str(THIS_DIR))

from dataset import BioSmokeDataset


def _normalize_reference_df(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize reference CSV so SAMPLE_ID is index when available."""
    out = df.copy()
    out.columns = [str(c).strip() for c in out.columns]

    if "SAMPLE_ID" in out.columns:
        out["SAMPLE_ID"] = out["SAMPLE_ID"].astype(str).str.strip()
        out = out.set_index("SAMPLE_ID")
    elif out.columns[0].lower().startswith("unnamed"):
        # Some CSV exports keep row index in the first unnamed column.
        out.iloc[:, 0] = out.iloc[:, 0].astype(str).str.strip()
        out = out.set_index(out.columns[0])

    out.index = out.index.astype(str).str.strip()
    return out


def _build_python_dataframe(
    sample_file: str,
    taxonomy_file: str,
    metabolomics_file: str,
    kegg_file: str,
    mode: str,
    load_kegg: bool,
) -> pd.DataFrame:
    ds = BioSmokeDataset(
        sample_file=sample_file,
        taxonomy_file=taxonomy_file,
        metabolomics_file=metabolomics_file,
        kegg_file=kegg_file,
        load_kegg=load_kegg,
    )

    if mode == "python-default":
        ds.preprocess_taxonomy_data()
        ds.preprocess_metabolomics_data()
        if ds.kegg_data is not None:
            ds.preprocess_kegg_data()
    elif mode == "r-like-processed":
        # Try to mimic R/dataset.R processed pipeline as closely as possible.
        ds.preprocess_taxonomy_data(
            remove_low_expression=True,
            min_abund_threshold=0.01,
            min_prevalence_threshold=0.1,
            remove_outliers=False,
            relative_abund=True,
            transform=True,
            transform_method="log",
        )
        ds.preprocess_metabolomics_data(
            remove_high_missing=False,
            remove_low_expression=True,
            min_abund_threshold=100.0,
            min_prevalence_threshold=0.1,
            remove_outliers=False,
            relative_abund=False,
            pqn_normalization=True,
            transform=True,
            transform_method="log",
            scale=False,
        )
        if ds.kegg_data is not None:
            ds.preprocess_kegg_data(
                remove_low_expression=True,
                min_abund_threshold=0.01,
                min_prevalence_threshold=0.1,
                remove_outliers=False,
                relative_abund=True,
                transform=True,
                transform_method="log",
            )
    else:
        raise ValueError(f"Unsupported mode: {mode}")

    merged = ds.merge_to_dataframe()
    merged.index = merged.index.astype(str).str.strip()
    return merged


def _coerce_numeric(series: pd.Series) -> Tuple[pd.Series, bool]:
    converted = pd.to_numeric(series, errors="coerce")
    is_numeric_like = converted.notna().sum() > 0
    return converted, is_numeric_like


def compare_dataframes(
    py_df: pd.DataFrame,
    ref_df: pd.DataFrame,
    abs_tol: float = 1e-8,
    rel_tol: float = 1e-6,
    topn: int = 20,
) -> Tuple[bool, List[str]]:
    report: List[str] = []

    def _prefix_of(col: str) -> str:
        if str(col).startswith("tax_"):
            return "tax"
        if str(col).startswith("met_"):
            return "met"
        if str(col).startswith("kegg_"):
            return "kegg"
        return "label_or_other"

    py_rows = set(py_df.index)
    ref_rows = set(ref_df.index)
    only_py_rows = sorted(py_rows - ref_rows)
    only_ref_rows = sorted(ref_rows - py_rows)

    py_cols = set(py_df.columns)
    ref_cols = set(ref_df.columns)
    only_py_cols = sorted(py_cols - ref_cols)
    only_ref_cols = sorted(ref_cols - py_cols)

    report.append(f"Python shape: {py_df.shape}, Reference shape: {ref_df.shape}")
    report.append(f"Rows only in Python: {len(only_py_rows)}")
    report.append(f"Rows only in Reference: {len(only_ref_rows)}")
    report.append(f"Cols only in Python: {len(only_py_cols)}")
    report.append(f"Cols only in Reference: {len(only_ref_cols)}")

    if only_py_rows:
        report.append("Sample rows only in Python (top 20): " + ", ".join(only_py_rows[:20]))
    if only_ref_rows:
        report.append("Sample rows only in Reference (top 20): " + ", ".join(only_ref_rows[:20]))
    if only_py_cols:
        report.append("Columns only in Python (top 20): " + ", ".join(only_py_cols[:20]))
    if only_ref_cols:
        report.append("Columns only in Reference (top 20): " + ", ".join(only_ref_cols[:20]))

    common_rows = sorted(py_rows & ref_rows)
    common_cols = [c for c in ref_df.columns if c in py_cols]

    if not common_rows or not common_cols:
        report.append("No comparable common rows/columns.")
        return False, report

    py_aligned = py_df.loc[common_rows, common_cols].copy()
    ref_aligned = ref_df.loc[common_rows, common_cols].copy()

    numeric_summary = []
    object_mismatch_total = 0

    for col in common_cols:
        py_num, py_num_like = _coerce_numeric(py_aligned[col])
        ref_num, ref_num_like = _coerce_numeric(ref_aligned[col])

        if py_num_like and ref_num_like:
            valid_mask = py_num.notna() & ref_num.notna()
            if valid_mask.sum() == 0:
                continue

            p = py_num[valid_mask].to_numpy(dtype=float)
            r = ref_num[valid_mask].to_numpy(dtype=float)
            diff = np.abs(p - r)
            tol = abs_tol + rel_tol * np.maximum(np.abs(p), np.abs(r))
            mismatch_mask = diff > tol
            mismatch_count = int(mismatch_mask.sum())
            if mismatch_count > 0:
                numeric_summary.append(
                    (
                        col,
                        mismatch_count,
                        float(diff.max(initial=0.0)),
                        float(diff.mean() if diff.size else 0.0),
                    )
                )
        else:
            py_obj = py_aligned[col].fillna("<NA>").astype(str).str.strip()
            ref_obj = ref_aligned[col].fillna("<NA>").astype(str).str.strip()
            mismatch_count = int((py_obj != ref_obj).sum())
            object_mismatch_total += mismatch_count
            if mismatch_count > 0:
                numeric_summary.append((col, mismatch_count, float("nan"), float("nan")))

    numeric_summary.sort(key=lambda x: x[1], reverse=True)

    total_mismatch_cells = sum(item[1] for item in numeric_summary)
    report.append(f"Total mismatched cells in common matrix: {total_mismatch_cells}")

    if numeric_summary:
        report.append(f"Top {min(topn, len(numeric_summary))} mismatched columns:")
        for col, count, max_abs, mean_abs in numeric_summary[:topn]:
            if np.isnan(max_abs):
                report.append(f"  - {col}: mismatches={count} (non-numeric/string compare)")
            else:
                report.append(
                    f"  - {col}: mismatches={count}, max_abs_diff={max_abs:.6g}, mean_abs_diff={mean_abs:.6g}"
                )

    # Summarize mismatch burden by feature prefix.
    prefix_mismatch = {"tax": 0, "met": 0, "kegg": 0, "label_or_other": 0}
    for col, count, _, _ in numeric_summary:
        prefix_mismatch[_prefix_of(col)] += int(count)
    report.append(
        "Mismatch by prefix: "
        + ", ".join([f"{k}={v}" for k, v in prefix_mismatch.items()])
    )

    prefix_only_py = {"tax": 0, "met": 0, "kegg": 0, "label_or_other": 0}
    for c in only_py_cols:
        prefix_only_py[_prefix_of(c)] += 1
    prefix_only_ref = {"tax": 0, "met": 0, "kegg": 0, "label_or_other": 0}
    for c in only_ref_cols:
        prefix_only_ref[_prefix_of(c)] += 1
    report.append(
        "Only-in-Python columns by prefix: "
        + ", ".join([f"{k}={v}" for k, v in prefix_only_py.items()])
    )
    report.append(
        "Only-in-Reference columns by prefix: "
        + ", ".join([f"{k}={v}" for k, v in prefix_only_ref.items()])
    )

    all_consistent = (
        len(only_py_rows) == 0
        and len(only_ref_rows) == 0
        and len(only_py_cols) == 0
        and len(only_ref_cols) == 0
        and total_mismatch_cells == 0
        and object_mismatch_total == 0
    )

    if all_consistent:
        report.append("RESULT: CONSISTENT")
    else:
        report.append("RESULT: NOT CONSISTENT")

    return all_consistent, report


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Check whether Python-built merged dataset matches dataset/merged_dataset_processed.csv"
    )
    parser.add_argument("--sample-file", default="dataset/id_sample.xlsx")
    parser.add_argument("--taxonomy-file", default="dataset/taxonomy_Species_abund.txt")
    parser.add_argument("--metabolomics-file", default="dataset/metabolome_annotated_data.csv")
    parser.add_argument("--kegg-file", default="dataset/KEGG/4_KOEntry/KEGG_KOEntry_abund.txt")
    parser.add_argument("--reference", default="dataset/merged_dataset_processed.csv")
    parser.add_argument(
        "--mode",
        choices=["python-default", "r-like-processed"],
        default="r-like-processed",
        help="Preprocessing mode used to generate Python dataframe for comparison.",
    )
    parser.add_argument("--disable-kegg", action="store_true", help="Disable KEGG import explicitly.")
    parser.add_argument("--abs-tol", type=float, default=1e-8)
    parser.add_argument("--rel-tol", type=float, default=1e-6)
    parser.add_argument("--topn", type=int, default=20)
    parser.add_argument("--report-out", default="results/python_r_consistency_report.txt")

    args = parser.parse_args()

    ref_path = Path(args.reference)
    if not ref_path.exists():
        raise FileNotFoundError(f"Reference file not found: {ref_path}")

    ref_df = pd.read_csv(ref_path, low_memory=False)
    ref_df = _normalize_reference_df(ref_df)

    inferred_need_kegg = any(str(c).startswith("kegg_") for c in ref_df.columns)
    load_kegg = (not args.disable_kegg) and inferred_need_kegg

    py_df = _build_python_dataframe(
        sample_file=args.sample_file,
        taxonomy_file=args.taxonomy_file,
        metabolomics_file=args.metabolomics_file,
        kegg_file=args.kegg_file,
        mode=args.mode,
        load_kegg=load_kegg,
    )

    ok, report = compare_dataframes(
        py_df=py_df,
        ref_df=ref_df,
        abs_tol=args.abs_tol,
        rel_tol=args.rel_tol,
        topn=args.topn,
    )

    report_text = "\n".join(report)
    print(report_text)

    out_path = Path(args.report_out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(report_text, encoding="utf-8")
    print(f"\nReport saved to: {out_path}")

    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
