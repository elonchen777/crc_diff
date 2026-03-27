import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.spatial.distance import pdist, squareform


def pcoa_from_bray(features: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Return 2D PCoA coordinates and explained variance from Bray-Curtis distances."""
    x = features.to_numpy(dtype=float)
    x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)

    dist_vec = pdist(x, metric="braycurtis")
    dist_mat = squareform(dist_vec)

    n = dist_mat.shape[0]
    eye = np.eye(n)
    one = np.ones((n, n)) / n
    j = eye - one
    b = -0.5 * j @ (dist_mat ** 2) @ j

    eigvals, eigvecs = np.linalg.eigh(b)
    order = np.argsort(eigvals)[::-1]
    eigvals = eigvals[order]
    eigvecs = eigvecs[:, order]

    pos = eigvals > 1e-12
    eigvals_pos = eigvals[pos]
    eigvecs_pos = eigvecs[:, pos]
    if eigvals_pos.size < 2:
        coords = np.zeros((n, 2))
        explained = np.array([0.0, 0.0])
        return coords, explained

    coords = eigvecs_pos[:, :2] * np.sqrt(eigvals_pos[:2])
    explained = eigvals_pos[:2] / eigvals_pos.sum()
    return coords, explained


def _build_design_matrix(values: pd.Series, factor_col: str) -> tuple[np.ndarray, int, int]:
    """Build design matrix for PERMANOVA and return X, n_levels, n_unique."""
    if factor_col == "age":
        x = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
        x = (x - np.nanmean(x)) / (np.nanstd(x) + 1e-12)
        x = np.nan_to_num(x, nan=0.0)
        x = x.reshape(-1, 1)
        n_levels = 1
        n_unique = int(np.unique(np.round(x, 6)).size)
    else:
        cat = values.astype("Int64").astype(str)
        x_df = pd.get_dummies(cat, drop_first=True)
        x = x_df.to_numpy(dtype=float)
        n_levels = x_df.shape[1]
        n_unique = int(cat.nunique())
    return x, n_levels, n_unique


def permanova_single_factor(
    features: pd.DataFrame,
    factor_values: pd.Series,
    factor_col: str,
    n_perm: int = 999,
    seed: int = 42,
) -> dict[str, float]:
    """PERMANOVA with one clinical factor, returning pseudo-F, p-value and R2."""
    y = factor_values.reset_index(drop=True)
    x_feat = features.reset_index(drop=True)
    mask = y.notna()
    y = y[mask]
    x_feat = x_feat.loc[mask]

    n = len(y)
    if n < 8:
        return {
            "n_samples": float(n),
            "n_levels": np.nan,
            "df_model": np.nan,
            "df_resid": np.nan,
            "r2": np.nan,
            "f_stat": np.nan,
            "p_value": np.nan,
        }

    x, _, n_unique = _build_design_matrix(y, factor_col)
    if x.shape[1] < 1 or n_unique < 2:
        return {
            "n_samples": float(n),
            "n_levels": float(n_unique),
            "df_model": np.nan,
            "df_resid": np.nan,
            "r2": np.nan,
            "f_stat": np.nan,
            "p_value": np.nan,
        }

    data = np.nan_to_num(x_feat.to_numpy(dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    d = squareform(pdist(data, metric="braycurtis"))

    eye = np.eye(n)
    one = np.ones((n, n)) / n
    j = eye - one
    b = -0.5 * j @ (d ** 2) @ j

    intercept = np.ones((n, 1))
    x_full = np.hstack([intercept, x])
    rank_x = np.linalg.matrix_rank(x_full)
    df_model = rank_x - 1
    df_resid = n - rank_x
    if df_model <= 0 or df_resid <= 0:
        return {
            "n_samples": float(n),
            "n_levels": float(n_unique),
            "df_model": np.nan,
            "df_resid": np.nan,
            "r2": np.nan,
            "f_stat": np.nan,
            "p_value": np.nan,
        }

    def _calc_f_and_r2(design: np.ndarray) -> tuple[float, float]:
        h = design @ np.linalg.pinv(design.T @ design) @ design.T
        ss_model = float(np.trace(h @ b))
        ss_total = float(np.trace(b))
        ss_resid = ss_total - ss_model
        if ss_total <= 0 or ss_resid <= 0:
            return np.nan, np.nan
        ms_model = ss_model / df_model
        ms_resid = ss_resid / df_resid
        f_val = ms_model / ms_resid if ms_resid > 0 else np.nan
        r2_val = ss_model / ss_total
        return f_val, r2_val

    f_obs, r2_obs = _calc_f_and_r2(x_full)
    if np.isnan(f_obs):
        p_val = np.nan
    else:
        rng = np.random.default_rng(seed)
        ge = 0
        valid_perm = 0
        for _ in range(n_perm):
            perm_idx = rng.permutation(n)
            x_perm = np.hstack([intercept, x[perm_idx, :]])
            f_perm, _ = _calc_f_and_r2(x_perm)
            if np.isnan(f_perm):
                continue
            valid_perm += 1
            if f_perm >= f_obs:
                ge += 1
        p_val = (ge + 1) / (valid_perm + 1) if valid_perm > 0 else np.nan

    return {
        "n_samples": float(n),
        "n_levels": float(n_unique),
        "df_model": float(df_model),
        "df_resid": float(df_resid),
        "r2": float(r2_obs),
        "f_stat": float(f_obs),
        "p_value": float(p_val),
    }


def load_and_merge_clinical(relative_csv: Path, clinical_xlsx: Path) -> pd.DataFrame:
    merged = pd.read_csv(relative_csv)
    clin = pd.read_excel(clinical_xlsx)

    if "SAMPLE_ID" not in merged.columns:
        raise ValueError("merged_dataset_relative.csv must contain SAMPLE_ID column")
    if "SAMPLE_ID" not in clin.columns:
        raise ValueError("id_sample.xlsx must contain SAMPLE_ID column")

    nerve_col = next((c for c in clin.columns if "NerveInvasion" in str(c)), None)
    vascular_col = next((c for c in clin.columns if "VascularInvasion" in str(c)), None)
    smoking_col = next((c for c in clin.columns if "0" in str(c) and "1" in str(c) and "����" in str(c)), None)

    needed = [c for c in [nerve_col, vascular_col, smoking_col] if c is not None]
    clin_sub = clin[["SAMPLE_ID", *needed]].copy()

    rename_map = {}
    if nerve_col:
        rename_map[nerve_col] = "nerve_invasion"
    if vascular_col:
        rename_map[vascular_col] = "vascular_invasion"
    if smoking_col:
        rename_map[smoking_col] = "smoking_binary"
    clin_sub = clin_sub.rename(columns=rename_map)

    out = merged.merge(clin_sub, on="SAMPLE_ID", how="left")

    if "smoking_label" in out.columns:
        out["smoking"] = out["smoking_label"]
        if "smoking_binary" in out.columns:
            out["smoking"] = out["smoking"].fillna(out["smoking_binary"])
    elif "smoking_binary" in out.columns:
        out["smoking"] = out["smoking_binary"]
    else:
        out["smoking"] = np.nan

    for col in ["tnm_stage", "differentiation", "smoking", "nerve_invasion", "vascular_invasion"]:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    if "age" in out.columns:
        out["age"] = pd.to_numeric(out["age"], errors="coerce")

    return out


def plot_factor_panels(
    df: pd.DataFrame,
    feature_sets: dict[str, list[str]],
    factor_col: str,
    output_path: Path,
) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(18, 5), constrained_layout=True)
    omics_names = ["Species", "Metabolites", "KO genes"]
    keys = ["tax", "met", "kegg"]

    for ax, key, title in zip(axes, keys, omics_names):
        cols = feature_sets[key]
        work = df[[factor_col, *cols]].copy()
        work = work.dropna(subset=[factor_col])
        if work.empty:
            ax.set_title(f"{title}: no data")
            ax.axis("off")
            continue

        coords, explained = pcoa_from_bray(work[cols])
        pcoa_df = pd.DataFrame(
            {
                "PCoA1": coords[:, 0],
                "PCoA2": coords[:, 1],
                factor_col: work[factor_col].to_numpy(),
            }
        )

        if factor_col == "age":
            points = ax.scatter(
                pcoa_df["PCoA1"],
                pcoa_df["PCoA2"],
                c=pcoa_df[factor_col],
                cmap="viridis",
                s=24,
                alpha=0.85,
                edgecolors="none",
            )
            cb = plt.colorbar(points, ax=ax)
            cb.set_label("Age")
        else:
            pcoa_df[factor_col] = pcoa_df[factor_col].astype("Int64").astype(str)
            sns.scatterplot(
                data=pcoa_df,
                x="PCoA1",
                y="PCoA2",
                hue=factor_col,
                palette="Set2",
                s=24,
                alpha=0.85,
                ax=ax,
            )
            ax.legend(title=factor_col, fontsize=8, title_fontsize=9, loc="best")

        ax.set_title(
            f"{title}\nPC1={explained[0]*100:.1f}%, PC2={explained[1]*100:.1f}%",
            fontsize=10,
        )
        ax.set_xlabel("PCoA1")
        ax.set_ylabel("PCoA2")

    fig.suptitle(f"CRC only: Bray-Curtis PCoA by clinical factor: {factor_col}", fontsize=13)
    fig.savefig(output_path, dpi=300)
    plt.close(fig)


def plot_permanova_summary(stats_df: pd.DataFrame, outdir: Path) -> None:
    work = stats_df.copy()
    work["r2_percent"] = work["r2"] * 100
    work["neglog10_p"] = -np.log10(work["p_value"].clip(lower=1e-300))

    r2_mat = work.pivot(index="factor", columns="omics", values="r2_percent")
    p_mat = work.pivot(index="factor", columns="omics", values="neglog10_p")

    plt.figure(figsize=(8, 5))
    sns.heatmap(r2_mat, annot=True, fmt=".2f", cmap="YlGnBu", cbar_kws={"label": "R2 (%)"})
    plt.title("CRC only: PERMANOVA effect size (R2%)")
    plt.xlabel("Omics")
    plt.ylabel("Clinical factor")
    plt.tight_layout()
    plt.savefig(outdir / "permanova_r2_heatmap_crc.png", dpi=300)
    plt.close()

    plt.figure(figsize=(8, 5))
    sns.heatmap(p_mat, annot=True, fmt=".2f", cmap="OrRd", cbar_kws={"label": "-log10(p)"})
    plt.title("CRC only: PERMANOVA significance (-log10 p)")
    plt.xlabel("Omics")
    plt.ylabel("Clinical factor")
    plt.tight_layout()
    plt.savefig(outdir / "permanova_significance_heatmap_crc.png", dpi=300)
    plt.close()

    plt.figure(figsize=(10, 5))
    sns.scatterplot(
        data=work,
        x="factor",
        y="omics",
        size="r2_percent",
        hue="p_value",
        palette="viridis_r",
        sizes=(60, 500),
        edgecolor="black",
    )
    plt.xticks(rotation=30, ha="right")
    plt.title("CRC only: PERMANOVA summary (size=R2%, color=p-value)")
    plt.xlabel("Clinical factor")
    plt.ylabel("Omics")
    plt.tight_layout()
    plt.savefig(outdir / "permanova_bubble_crc.png", dpi=300)
    plt.close()


def plot_explained_variance_bubble(explained_df: pd.DataFrame, outdir: Path) -> None:
    """Bubble plot for explained variance (R2%) by omics group and clinical factor."""
    work = explained_df.copy()

    plt.figure(figsize=(10, 5))
    ax = sns.scatterplot(
        data=work,
        x="factor",
        y="omics",
        size="explained_variance_percent",
        hue="explained_variance_percent",
        palette="YlGnBu",
        sizes=(80, 1200),
        edgecolor="black",
    )
    plt.xticks(rotation=30, ha="right")
    plt.title("CRC only: explained variance bubble plot (R2%)")
    plt.xlabel("Clinical factor")
    plt.ylabel("Omics")

    for _, row in work.iterrows():
        ax.text(
            x=row["factor"],
            y=row["omics"],
            s=f"{row['explained_variance_percent']:.2f}",
            ha="center",
            va="center",
            fontsize=8,
            color="black",
        )

    plt.tight_layout()
    plt.savefig(outdir / "explained_variance_bubble_crc.png", dpi=300)
    plt.close()


def _sig_label(p: float) -> str:
    if pd.isna(p):
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="CRC-only clinical factors with species/metabolites/KO composition via Bray-Curtis."
    )
    parser.add_argument(
        "--relative",
        default="dataset/merged_dataset_relative.csv",
        help="Path to merged relative abundance table",
    )
    parser.add_argument(
        "--clinical",
        default="dataset/id_sample.xlsx",
        help="Path to clinical xlsx",
    )
    parser.add_argument(
        "--outdir",
        default="results/pcoa_plots/clinical_factors_bray_crc",
        help="Output directory",
    )
    parser.add_argument(
        "--n-perm",
        type=int,
        default=999,
        help="Permutation count for PERMANOVA",
    )
    args = parser.parse_args()

    relative_path = Path(args.relative)
    clinical_path = Path(args.clinical)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_and_merge_clinical(relative_path, clinical_path)

    if "crc_label" in df.columns:
        df = df[df["crc_label"] == 1].copy()
    else:
        df = df[df["SAMPLE_ID"].astype(str).str.contains("CRC")].copy()
    df = df.reset_index(drop=True)

    feature_sets = {
        "tax": [c for c in df.columns if c.startswith("tax_")],
        "met": [c for c in df.columns if c.startswith("met_")],
        "kegg": [c for c in df.columns if c.startswith("kegg_")],
    }

    factors = ["age", "tnm_stage", "differentiation", "nerve_invasion", "vascular_invasion", "smoking"]
    existing_factors = [f for f in factors if f in df.columns]

    summary_rows = []
    permanova_rows = []
    omics_map = {
        "tax": "Species",
        "met": "Metabolites",
        "kegg": "KO genes",
    }

    for factor in existing_factors:
        non_null = int(df[factor].notna().sum())
        unique_n = int(df[factor].nunique(dropna=True))
        summary_rows.append({"factor": factor, "non_null_samples": non_null, "n_unique": unique_n})

        out_png = outdir / f"bray_pcoa_{factor}.png"
        plot_factor_panels(df, feature_sets, factor, out_png)

        for k, cols in feature_sets.items():
            stats = permanova_single_factor(
                features=df[cols],
                factor_values=df[factor],
                factor_col=factor,
                n_perm=args.n_perm,
                seed=42,
            )
            permanova_rows.append(
                {
                    "omics": omics_map[k],
                    "factor": factor,
                    **stats,
                }
            )

    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(outdir / "clinical_factor_summary.csv", index=False)

    stats_df = pd.DataFrame(permanova_rows)
    if not stats_df.empty:
        stats_df["significance"] = stats_df["p_value"].apply(_sig_label)
        stats_df = stats_df.sort_values(["factor", "omics"]).reset_index(drop=True)
        stats_df.to_csv(outdir / "permanova_crc_clinical_factors.csv", index=False)

        explained_df = stats_df[
            ["omics", "factor", "n_samples", "df_model", "r2", "p_value", "significance"]
        ].copy()
        explained_df["explained_variance_percent"] = explained_df["r2"] * 100
        denom = explained_df["n_samples"] - explained_df["df_model"] - 1
        explained_df["adjusted_r2"] = np.where(
            denom > 0,
            1 - (1 - explained_df["r2"]) * (explained_df["n_samples"] - 1) / denom,
            np.nan,
        )
        explained_df["adjusted_explained_variance_percent"] = explained_df["adjusted_r2"] * 100
        explained_df = explained_df.sort_values(["factor", "omics"]).reset_index(drop=True)
        explained_df.to_csv(outdir / "explained_variance_by_group_factor_crc.csv", index=False)

        plot_permanova_summary(stats_df, outdir)
        plot_explained_variance_bubble(explained_df, outdir)

    df[["SAMPLE_ID", *existing_factors]].to_csv(outdir / "clinical_factors_merged.csv", index=False)

    print("Analysis done.")
    print(f"CRC samples kept: {len(df)}")
    print(f"Output directory: {outdir}")
    print("Generated files:")
    for p in sorted(outdir.glob("*")):
        print(f" - {p.name}")


if __name__ == "__main__":
    main()
