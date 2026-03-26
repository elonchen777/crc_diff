#!/usr/bin/env python3
"""
Build CRC OS/PFS survival tables from clinical Excel data.

Outputs two datasets:
1) Full table with source fields + derived survival fields
2) Minimal table with patient ID and survival fields: SAMPLE_ID, OS_time_days, OS_event, PFS_time_days, PFS_event
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


DEFAULT_INPUT = Path("dataset/id_sample.xlsx")
DEFAULT_SHEET = "Sheet1"
DEFAULT_OUT_DIR = Path("dataset/survival")

# Column mapping for the uploaded clinical sheet
COL_ID = "SAMPLE_ID"
COL_GROUP = "分组"
COL_ENROLL = "手术日期"  # used as enrollment date
COL_DEATH_FLAG = "是否死亡(live,0;death,1)"
COL_DEATH_DATE = "死亡日期"
COL_RECUR_FLAG = "是否复发(Neg,0;Pos,1)"
COL_RECUR_DATE = "复发时间"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Generate CRC OS/PFS full and min-4 survival datasets")
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT, help="Input Excel file path")
    parser.add_argument("--sheet", type=str, default=DEFAULT_SHEET, help="Input sheet name")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR, help="Output directory")
    parser.add_argument(
        "--group-value",
        type=str,
        default="CRC",
        help="Keep rows where group column contains this value (case-insensitive)",
    )
    parser.add_argument(
        "--no-group-filter",
        action="store_true",
        help="Disable filtering by group column",
    )
    return parser.parse_args()


def parse_mixed_date(series: pd.Series) -> pd.Series:
    """Parse mixed-format date strings into datetime."""
    text = series.astype(str).str.strip()
    text = text.replace(
        {
            "": np.nan,
            "nan": np.nan,
            "NaT": np.nan,
            "NS": np.nan,
            "x": np.nan,
            "?": np.nan,
            "未": np.nan,
            "无": np.nan,
        }
    )
    text = text.str.replace("年", "-", regex=False).str.replace("月", "-", regex=False).str.replace("日", "", regex=False)
    text = text.str.replace("/", "-", regex=False).str.replace(".", "-", regex=False)
    return pd.to_datetime(text, errors="coerce")


def parse_binary_flag(series: pd.Series) -> pd.Series:
    """Convert event flag column to 0/1 integers, non-numeric values become 0."""
    return pd.to_numeric(series, errors="coerce").fillna(0).astype(int)


def build_survival_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.Timestamp]:
    work = df.copy()

    work["enroll_date"] = parse_mixed_date(work[COL_ENROLL])
    work["death_date"] = parse_mixed_date(work[COL_DEATH_DATE])
    work["recur_date"] = parse_mixed_date(work[COL_RECUR_DATE])

    death_flag = parse_binary_flag(work[COL_DEATH_FLAG])

    # Administrative censor date: latest date across all death/recurrence dates
    followup_cutoff = pd.concat([work["death_date"], work["recur_date"]], axis=0).max()

    # OS definitions
    os_event = ((death_flag == 1) | work["death_date"].notna()).astype(int)
    os_end_date = work["death_date"].where(work["death_date"].notna(), followup_cutoff)
    os_time_days = (os_end_date - work["enroll_date"]).dt.days

    # PFS definitions: min(recurrence date, death date, followup cutoff) - enrollment date
    min_event_date = pd.concat([work["recur_date"], work["death_date"]], axis=1).min(axis=1)
    pfs_end_date = min_event_date.where(min_event_date.notna(), followup_cutoff)
    pfs_event = min_event_date.notna().astype(int)
    pfs_time_days = (pfs_end_date - work["enroll_date"]).dt.days

    full = pd.DataFrame(
        {
            "SAMPLE_ID": work[COL_ID],
            "入组日期": work["enroll_date"].dt.date,
            "是否死亡_原始": work[COL_DEATH_FLAG],
            "死亡日期_原始": work[COL_DEATH_DATE],
            "是否复发_原始": work[COL_RECUR_FLAG],
            "复发日期_原始": work[COL_RECUR_DATE],
            "死亡事件_OS(0/1)": os_event,
            "OS终点日期": os_end_date.dt.date,
            "OS_time_days": os_time_days,
            "PFS事件(0/1)": pfs_event,
            "PFS终点日期": pfs_end_date.dt.date,
            "PFS_time_days": pfs_time_days,
            "随访截止日期_全队列统一": followup_cutoff.date() if pd.notna(followup_cutoff) else pd.NaT,
        }
    )

    # Keep valid rows only
    full = full[full["入组日期"].notna()].copy()
    full = full[(full["OS_time_days"] >= 0) & (full["PFS_time_days"] >= 0)].copy()

    min4 = full[["SAMPLE_ID", "OS_time_days", "死亡事件_OS(0/1)", "PFS_time_days", "PFS事件(0/1)"]].rename(
        columns={
            "死亡事件_OS(0/1)": "OS_event",
            "PFS事件(0/1)": "PFS_event",
        }
    )

    return full, min4, followup_cutoff


def save_tables(full: pd.DataFrame, min4: pd.DataFrame, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)

    full_csv = out_dir / "crc_os_pfs_full.csv"
    full_xlsx = out_dir / "crc_os_pfs_full.xlsx"
    min4_csv = out_dir / "crc_os_pfs_min4.csv"
    min4_xlsx = out_dir / "crc_os_pfs_min4.xlsx"

    full.to_csv(full_csv, index=False, encoding="utf-8-sig")
    full.to_excel(full_xlsx, index=False)
    min4.to_csv(min4_csv, index=False, encoding="utf-8-sig")
    min4.to_excel(min4_xlsx, index=False)

    print(f"Saved: {full_csv}")
    print(f"Saved: {full_xlsx}")
    print(f"Saved: {min4_csv}")
    print(f"Saved: {min4_xlsx}")


def main() -> None:
    args = parse_args()

    df = pd.read_excel(args.input, sheet_name=args.sheet)

    if not args.no_group_filter and COL_GROUP in df.columns:
        df = df[df[COL_GROUP].astype(str).str.contains(args.group_value, case=False, na=False)].copy()

    full, min4, followup_cutoff = build_survival_tables(df)
    save_tables(full, min4, args.out_dir)

    print(f"followup_cutoff = {followup_cutoff}")
    print(f"n = {len(full)}")
    print(f"OS events = {int(full['死亡事件_OS(0/1)'].sum())}")
    print(f"PFS events = {int(full['PFS事件(0/1)'].sum())}")


if __name__ == "__main__":
    main()
