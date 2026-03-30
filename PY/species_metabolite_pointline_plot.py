#!/usr/bin/env python3
"""Plot species-metabolite regression panels by group."""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

CURRENT_DIR = Path(__file__).resolve().parent
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

from dataset import BioSmokeDataset
from split_group import prepare_crc_diff_groups
from species_pointline_plot import (
    GROUP_TITLE,
    GROUP_ORDER,
    PLOT_STYLE,
    _canonicalize_label,
    _find_species_column,
    _italicize_species,
    _normalize_name,
    _plot_group_panel,
)


def _find_metabolite_column(df: pd.DataFrame, target_name: str) -> str:
    metabolite_cols = [col for col in df.columns if col.startswith("met_")]
    target_key = _canonicalize_label(target_name)

    exact_matches = [col for col in metabolite_cols if _canonicalize_label(col) == target_key]
    if len(exact_matches) == 1:
        return exact_matches[0]
    if len(exact_matches) > 1:
        raise ValueError(f"Metabolite {target_name} matches multiple columns: {exact_matches}")

    fallback_matches = [
        col for col in metabolite_cols
        if target_key and target_key in _canonicalize_label(col)
    ]
    if len(fallback_matches) == 1:
        return fallback_matches[0]
    if len(fallback_matches) > 1:
        raise ValueError(f"Metabolite {target_name} matches multiple columns: {fallback_matches}")

    available = ", ".join(metabolite_cols[:20])
    raise ValueError(
        f"Cannot find metabolite column for {target_name}. Available examples: {available}"
    )


def _metabolite_label(name: str) -> str:
    text = str(name).strip().strip('"')
    text = re.sub(r"^(met_|tax_|kegg_)", "", text, flags=re.IGNORECASE)
    return text.replace("_", " ")


def _prepare_species_metabolite_group_data(
    df: pd.DataFrame,
    species: str,
    metabolite: str,
    transform: str,
) -> tuple[pd.DataFrame, str, str, str, str]:
    x_col = _find_species_column(df, species)
    y_col = _find_metabolite_column(df, metabolite)

    plot_df = df[[x_col, y_col, "group"]].copy()
    plot_df[x_col] = pd.to_numeric(plot_df[x_col], errors="coerce")
    plot_df[y_col] = pd.to_numeric(plot_df[y_col], errors="coerce")
    plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna()

    plot_df = plot_df[plot_df["group"].isin(GROUP_ORDER)].copy()
    plot_df["group"] = pd.Categorical(plot_df["group"], categories=GROUP_ORDER, ordered=True)

    if transform == "log1p":
        plot_df[x_col] = np.log1p(plot_df[x_col].clip(lower=0))
        plot_df[y_col] = np.log1p(plot_df[y_col].clip(lower=0))
        x_label = f"{_italicize_species(species)} (log1p abundance)"
        y_label = f"{_metabolite_label(metabolite)} (log1p abundance)"
    elif transform == "none":
        x_label = _italicize_species(species)
        y_label = _metabolite_label(metabolite)
    else:
        raise ValueError("transform only supports log1p or none")

    return plot_df, x_col, y_col, x_label, y_label


def _save_single_group_plot(
    plot_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    species: str,
    metabolite: str,
    group_name: str,
    output_dir: Path,
    x_label: str,
    y_label: str,
) -> Path:
    sns.set_theme(style="white", rc=PLOT_STYLE)
    fig, ax = plt.subplots(figsize=(6.9, 6.0), dpi=160)
    _plot_group_panel(ax, plot_df, x_col, y_col, group_name, x_label, y_label)
    display_group = GROUP_TITLE.get(group_name, group_name)
    ax.set_title(
        f"{display_group}: {_italicize_species(species)} vs {_metabolite_label(metabolite)}",
        pad=14,
    )
    fig.tight_layout()

    safe_species = _normalize_name(species)
    safe_metabolite = _normalize_name(metabolite)
    plot_file = output_dir / (
        f"{safe_species}_vs_{safe_metabolite}_{group_name.lower()}_pointline.png"
    )
    fig.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return plot_file


def build_species_metabolite_plot(
    species: str,
    metabolite: str,
    output_dir: str,
    transform: str = "none",
) -> Path:
    ds = BioSmokeDataset()
    ds.preprocess_taxonomy_data(transform=True)
    ds.preprocess_metabolomics_data(transform=True)
    df = ds.merge_to_dataframe()
    df = prepare_crc_diff_groups(df)

    plot_df, x_col, y_col, x_label, y_label = _prepare_species_metabolite_group_data(
        df,
        species,
        metabolite,
        transform,
    )

    if len(plot_df) < 3:
        raise ValueError("Less than 3 points are available for plotting")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    safe_species = _normalize_name(species)
    safe_metabolite = _normalize_name(metabolite)

    data_file = output_path / f"{safe_species}_vs_{safe_metabolite}_pointline_data.csv"
    grouped_file = output_path / f"{safe_species}_vs_{safe_metabolite}_pointline_grouped.png"

    plot_df.to_csv(data_file, index=False)

    sns.set_theme(style="white", rc=PLOT_STYLE)
    fig, axes = plt.subplots(1, 3, figsize=(18.5, 6.0), dpi=160, sharex=True, sharey=True)
    if len(GROUP_ORDER) == 1:
        axes = [axes]

    for idx, group_name in enumerate(GROUP_ORDER):
        _plot_group_panel(
            axes[idx],
            plot_df,
            x_col,
            y_col,
            group_name,
            x_label,
            y_label,
            show_xlabel=True,
            show_ylabel=(idx == 0),
        )

    fig.suptitle(
        f"{_italicize_species(species)} vs {_metabolite_label(metabolite)}",
        fontsize=17,
        fontweight="bold",
        y=1.03,
    )
    fig.tight_layout()
    fig.savefig(grouped_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    single_files = []
    for group_name in GROUP_ORDER:
        single_files.append(
            _save_single_group_plot(
                plot_df,
                x_col,
                y_col,
                species,
                metabolite,
                group_name,
                output_path,
                x_label,
                y_label,
            )
        )

    print(f"Saved grouped plot: {grouped_file}")
    for plot_file in single_files:
        print(f"Saved single-group plot: {plot_file}")
    print(f"Saved data: {data_file}")
    return grouped_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot two species against one metabolite with grouped regression panels"
    )
    parser.add_argument("--species-a", default="Peptostreptococcus_stomatis", help="First species")
    parser.add_argument(
        "--species-b", default="Faecalibacterium_prausnitzii", help="Second species"
    )
    parser.add_argument("--metabolite", default="Cholesterol", help="Target metabolite")
    parser.add_argument(
        "--output-dir",
        default="results/species_metabolite_pointline",
        help="Output directory",
    )
    parser.add_argument(
        "--transform",
        choices=["log1p", "none"],
        default="none",
        help="Value transform before plotting",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_species_metabolite_plot(
        species=args.species_a,
        metabolite=args.metabolite,
        output_dir=args.output_dir,
        transform=args.transform,
    )
    build_species_metabolite_plot(
        species=args.species_b,
        metabolite=args.metabolite,
        output_dir=args.output_dir,
        transform=args.transform,
    )


if __name__ == "__main__":
    main()
