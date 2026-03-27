#!/usr/bin/env python3
"""绘制两个物种丰度的散点回归图，并标注相关性与置信区间。"""

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
from scipy import stats

CURRENT_DIR = Path(__file__).resolve().parent
if str(CURRENT_DIR) not in sys.path:
    sys.path.insert(0, str(CURRENT_DIR))

from dataset import BioSmokeDataset
from split_group import prepare_crc_diff_groups

GROUP_COLORS = {
    "CRC_poordiff": "#FF6B6B",
    "CRC_welldiff": "#FFD166",
    "CTRL": "#06D6A0",
}
GROUP_ORDER = ["CTRL", "CRC_welldiff", "CRC_poordiff"]
GROUP_TITLE = {
    "CTRL": "CTRL",
    "CRC_welldiff": "CRC-Well",
    "CRC_poordiff": "CRC-Poor",
}

PLOT_STYLE = {
    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "axes.edgecolor": "#202020",
    "axes.labelcolor": "#202020",
    "xtick.color": "#202020",
    "ytick.color": "#202020",
    "axes.grid": False,
    "font.size": 11,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
}


def _normalize_name(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", str(name).lower())


def _italicize_species(name: str) -> str:
    text = str(name).strip().strip('"')
    text = re.sub(r"^(tax_|met_|kegg_)", "", text, flags=re.IGNORECASE)
    text = re.sub(r"^[a-z]__", "", text, flags=re.IGNORECASE)
    text = text.replace("_", " ")
    text = text.replace("\\", r"\\")
    text = text.replace(" ", r"\ ")
    return fr"$\it{{{text}}}$"


def _canonicalize_label(label: str) -> str:
    text = str(label).strip().strip('"')
    text = re.sub(r"^(tax_|met_|kegg_)", "", text, flags=re.IGNORECASE)
    text = re.sub(r"^[a-z]__", "", text, flags=re.IGNORECASE)
    return _normalize_name(text)


def _find_species_column(df: pd.DataFrame, target_name: str) -> str:
    species_cols = [col for col in df.columns if col.startswith("tax_")]
    target_key = _canonicalize_label(target_name)

    exact_matches = [col for col in species_cols if _canonicalize_label(col) == target_key]
    if len(exact_matches) == 1:
        return exact_matches[0]
    if len(exact_matches) > 1:
        raise ValueError(f"物种 {target_name} 匹配到多个列: {exact_matches}")

    fallback_matches = [
        col for col in species_cols
        if target_key and target_key in _canonicalize_label(col)
    ]
    if len(fallback_matches) == 1:
        return fallback_matches[0]
    if len(fallback_matches) > 1:
        raise ValueError(f"物种 {target_name} 匹配到多个列: {fallback_matches}")

    available = ", ".join(species_cols[:20])
    raise ValueError(
        f"未找到物种 {target_name} 对应的列。可用物种列示例: {available}"
    )


def _pearson_ci(x: pd.Series, y: pd.Series, confidence: float = 0.95) -> tuple[float, float, float, float]:
    r, p_value = stats.pearsonr(x, y)
    n = len(x)

    if n <= 3 or np.isclose(abs(r), 1.0):
        return r, p_value, np.nan, np.nan

    clipped_r = float(np.clip(r, -0.999999, 0.999999))
    z = np.arctanh(clipped_r)
    se = 1.0 / np.sqrt(n - 3)
    z_crit = stats.norm.ppf((1.0 + confidence) / 2.0)
    ci_low = np.tanh(z - z_crit * se)
    ci_high = np.tanh(z + z_crit * se)
    return r, p_value, ci_low, ci_high


def _prepare_group_data(
    df: pd.DataFrame,
    species_a: str,
    species_b: str,
    transform: str = "log1p",
) -> tuple[pd.DataFrame, str, str]:
    x_col = _find_species_column(df, species_a)
    y_col = _find_species_column(df, species_b)

    plot_df = df[[x_col, y_col, "group"]].copy()
    plot_df[x_col] = pd.to_numeric(plot_df[x_col], errors="coerce")
    plot_df[y_col] = pd.to_numeric(plot_df[y_col], errors="coerce")
    plot_df = plot_df.replace([np.inf, -np.inf], np.nan).dropna()

    plot_df = plot_df[plot_df["group"].isin(GROUP_ORDER)].copy()
    plot_df["group"] = pd.Categorical(plot_df["group"], categories=GROUP_ORDER, ordered=True)

    if transform == "log1p":
        plot_df[x_col] = np.log1p(plot_df[x_col].clip(lower=0))
        plot_df[y_col] = np.log1p(plot_df[y_col].clip(lower=0))
        x_label = f"{_italicize_species(species_a)} (log1p abundance)"
        y_label = f"{_italicize_species(species_b)} (log1p abundance)"
    elif transform == "none":
        x_label = _italicize_species(species_a)
        y_label = _italicize_species(species_b)
    else:
        raise ValueError("transform 仅支持 log1p 或 none")

    return plot_df, x_col, y_col, x_label, y_label


def _plot_group_panel(
    ax: plt.Axes,
    plot_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    group_name: str,
    x_label: str,
    y_label: str,
    show_xlabel: bool = True,
    show_ylabel: bool = True,
) -> None:
    group_df = plot_df[plot_df["group"] == group_name].copy()
    color = GROUP_COLORS[group_name]

    if len(group_df) < 3:
        ax.text(0.5, 0.5, "数据不足", transform=ax.transAxes, ha="center", va="center", fontsize=12)
        ax.set_axis_off()
        return

    r, p_value, ci_low, ci_high = _pearson_ci(group_df[x_col], group_df[y_col])

    sns.regplot(
        data=group_df,
        x=x_col,
        y=y_col,
        ci=95,
        n_boot=1000,
        color=color,
        scatter_kws={
            "s": 44,
            "alpha": 0.82,
            "color": color,
            "edgecolor": "white",
            "linewidths": 0.45,
        },
        line_kws={"color": color, "lw": 2.6},
        ax=ax,
    )

    ax.set_title(f"{GROUP_TITLE[group_name]} (n={len(group_df)})", pad=12)
    ax.set_xlabel(x_label if show_xlabel else "")
    ax.set_ylabel(y_label if show_ylabel else "")

    annotation = (
        f"r = {r:.3f}\n"
        f"95% CI = [{ci_low:.3f}, {ci_high:.3f}]\n"
        f"p = {p_value:.2e}"
    )
    ax.text(
        0.04,
        0.96,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10.5,
        bbox={"facecolor": "white", "alpha": 0.88, "edgecolor": color, "linewidth": 0.9, "pad": 5},
    )

    ax.grid(False)
    ax.set_facecolor("white")
    ax.set_frame_on(True)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.25)
        spine.set_color("#202020")
    ax.tick_params(axis="both", labelsize=10)


def _save_single_group_plot(
    plot_df: pd.DataFrame,
    x_col: str,
    y_col: str,
    species_a: str,
    species_b: str,
    group_name: str,
    output_dir: Path,
    x_label: str,
    y_label: str,
) -> Path:
    sns.set_theme(style="white", rc=PLOT_STYLE)
    fig, ax = plt.subplots(figsize=(6.9, 6.0), dpi=160)
    _plot_group_panel(ax, plot_df, x_col, y_col, group_name, x_label, y_label)
    ax.set_title(
        f"{GROUP_TITLE[group_name]}: {_italicize_species(species_a)} vs {_italicize_species(species_b)}",
        pad=14,
    )
    fig.tight_layout()

    safe_a = _normalize_name(species_a)
    safe_b = _normalize_name(species_b)
    plot_file = output_dir / f"{safe_a}_vs_{safe_b}_{group_name.lower()}_pointline.png"
    fig.savefig(plot_file, dpi=300, bbox_inches="tight")
    plt.close(fig)
    return plot_file


def build_species_pair_plot(
    species_a: str,
    species_b: str,
    output_dir: str,
    transform: str = "log1p",
) -> Path:
    ds = BioSmokeDataset()
    df = ds.merge_to_dataframe()
    df = prepare_crc_diff_groups(df)

    plot_df, x_col, y_col, x_label, y_label = _prepare_group_data(df, species_a, species_b, transform)

    if len(plot_df) < 3:
        raise ValueError("可用于绘图的数据点少于 3 个，无法计算相关性")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    safe_a = _normalize_name(species_a)
    safe_b = _normalize_name(species_b)

    data_file = output_path / f"{safe_a}_vs_{safe_b}_pointline_data.csv"
    combined_file = output_path / f"{safe_a}_vs_{safe_b}_pointline_grouped.png"

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
        f"{_italicize_species(species_a)} vs {_italicize_species(species_b)}",
        fontsize=17,
        fontweight="bold",
        y=1.03,
    )
    fig.tight_layout()
    fig.savefig(combined_file, dpi=300, bbox_inches="tight")
    plt.close(fig)

    single_files = []
    for group_name in GROUP_ORDER:
        single_files.append(
            _save_single_group_plot(
                plot_df,
                x_col,
                y_col,
                species_a,
                species_b,
                group_name,
                output_path,
                x_label,
                y_label,
            )
        )

    print(f"已保存总览图像: {combined_file}")
    for plot_file in single_files:
        print(f"已保存分组图像: {plot_file}")
    print(f"已保存数据: {data_file}")
    return combined_file


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="绘制两个物种的点线相关图")
    parser.add_argument(
        "--species-a",
        default="Peptostreptococcus_stomatis",
        help="第一个物种名称",
    )
    parser.add_argument(
        "--species-b",
        default="Faecalibacterium_prausnitzii",
        help="第二个物种名称",
    )
    parser.add_argument(
        "--output-dir",
        default="results/species_pointline",
        help="输出目录",
    )
    parser.add_argument(
        "--transform",
        choices=["log1p", "none"],
        default="log1p",
        help="绘图前的数值变换方式",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    build_species_pair_plot(
        species_a=args.species_a,
        species_b=args.species_b,
        output_dir=args.output_dir,
        transform=args.transform,
    )


if __name__ == "__main__":
    main()