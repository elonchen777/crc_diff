from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import chi2, chi2_contingency, mannwhitneyu


sns.set_theme(
    style="whitegrid",
    context="talk",
    rc={
        "figure.dpi": 120,
        "savefig.dpi": 400,
        "axes.titleweight": "bold",
        "axes.labelsize": 13,
        "axes.titlesize": 15,
        "legend.frameon": False,
        "grid.linestyle": "--",
        "grid.alpha": 0.18,
    },
)


def format_p(p):
    if p is None or not np.isfinite(p):
        return "NA"
    if p < 1e-4:
        return "<1e-4"
    return f"{p:.3g}"


def parse_mixed_date(series: pd.Series) -> pd.Series:
    s = series.astype(str).str.strip()
    s = s.replace({
        "": np.nan,
        "NA": np.nan,
        "NaN": np.nan,
        "nan": np.nan,
        "NS": np.nan,
        "x": np.nan,
        "?": np.nan,
        "未": np.nan,
        "无": np.nan,
    })
    s = s.str.replace("年", "-", regex=False).str.replace("月", "-", regex=False).str.replace("日", "", regex=False)
    s = s.str.replace("/", "-", regex=False).str.replace(".", "-", regex=False)
    return pd.to_datetime(s, errors="coerce")


def km_curve(times: np.ndarray, events: np.ndarray):
    order = np.argsort(times)
    t = times[order]
    e = events[order]

    uniq_event_times = np.unique(t[e == 1])
    if len(uniq_event_times) == 0:
        return np.array([0.0]), np.array([1.0])

    surv = 1.0
    xs = [0.0]
    ys = [1.0]

    for ti in uniq_event_times:
        at_risk = np.sum(t >= ti)
        d_i = np.sum((t == ti) & (e == 1))
        if at_risk > 0:
            surv *= (1.0 - d_i / at_risk)
        xs.extend([ti, ti])
        ys.extend([ys[-1], surv])

    return np.array(xs), np.array(ys)


def weighted_logrank_test(times, events, groups, weight="logrank"):
    df = pd.DataFrame({"time": times, "event": events, "group": groups}).dropna()
    if df.empty or df["group"].nunique() != 2:
        return np.nan

    g_levels = list(df["group"].unique())
    g0, g1 = g_levels[0], g_levels[1]

    event_times = np.sort(df.loc[df["event"] == 1, "time"].unique())
    if len(event_times) == 0:
        return np.nan

    O1 = 0.0
    E1 = 0.0
    V1 = 0.0

    # Kaplan-Meier for pooled sample for Peto/Tarone weights
    if weight in ("breslow", "tarone"):
        xs_pool, ys_pool = km_curve(df["time"].values, df["event"].values)

        def s_before(t):
            idx = np.searchsorted(xs_pool, t, side="right") - 1
            if idx < 0:
                return 1.0
            return float(ys_pool[idx])

    for t in event_times:
        n1 = np.sum((df["time"] >= t) & (df["group"] == g1))
        n0 = np.sum((df["time"] >= t) & (df["group"] == g0))
        n = n1 + n0
        d1 = np.sum((df["time"] == t) & (df["event"] == 1) & (df["group"] == g1))
        d0 = np.sum((df["time"] == t) & (df["event"] == 1) & (df["group"] == g0))
        d = d1 + d0

        if n <= 1 or d == 0:
            continue

        if weight == "logrank":
            w = 1.0
        elif weight == "breslow":
            w = s_before(t)
        elif weight == "tarone":
            w = np.sqrt(max(s_before(t), 0.0))
        else:
            w = 1.0

        e1 = d * (n1 / n)
        v1 = (n1 * n0 * d * (n - d)) / (n * n * (n - 1))

        O1 += w * d1
        E1 += w * e1
        V1 += (w * w) * v1

    if V1 <= 0:
        return np.nan

    z2 = (O1 - E1) ** 2 / V1
    return float(chi2.sf(z2, df=1))


def shared_xmax(df, time_col, group_col):
    values = []
    for g in pd.unique(df[group_col].dropna()):
        sub = df[df[group_col] == g]
        if len(sub) == 0:
            continue
        max_time = pd.to_numeric(sub[time_col], errors="coerce").max()
        if pd.notna(max_time):
            values.append(float(max_time))
    if not values:
        return None
    return max(0.0, min(values))


def make_km_plot(
    df,
    time_col,
    event_col,
    group_col,
    title,
    subtitle,
    out_png,
    palette_override=None,
    legend_title="Group",
    y_label="Survival probability",
):
    fig, (ax, ax_risk) = plt.subplots(
        2,
        1,
        figsize=(8.6, 7.2),
        sharex=True,
        gridspec_kw={"height_ratios": [4.6, 1.4], "hspace": 0.06},
    )
    uniq_groups = [g for g in pd.unique(df[group_col]) if pd.notna(g)]
    if len(uniq_groups) == 0:
        plt.close(fig)
        return
    if len(uniq_groups) == 2:
        group_order = ["Low", "High"] if set(uniq_groups) == {"Low", "High"} else sorted(uniq_groups)
    else:
        group_order = sorted(uniq_groups)

    if palette_override is not None:
        palette = {g: palette_override.get(g, "#4C4C4C") for g in group_order}
    else:
        color_seq = ["#2166AC", "#B2182B", "#1B9E77", "#D95F02"]
        palette = {g: color_seq[i % len(color_seq)] for i, g in enumerate(group_order)}

    for g in group_order:
        sub = df[df[group_col] == g]
        if sub.empty:
            continue
        # Extract times and events
        times = sub[time_col].values.astype(float)
        events = sub[event_col].values.astype(int)

        # Plot KM curve
        x, y = km_curve(times, events)
        ax.step(x, y, where="post", label=g, color=palette[g], linewidth=2.4)

        # Plot census marks (vertical ticks for censored patients)
        censored_times = times[events == 0]
        if len(censored_times) > 0:
            # Find the survival probability at each censored time
            # We use the 'y' array which has doubled points [t_i, t_i] for steps
            census_y = []
            for ct in censored_times:
                idx = np.searchsorted(x, ct, side="right") - 1
                if idx < 0:
                    census_y.append(1.0)
                else:
                    census_y.append(y[idx])
            ax.scatter(
                censored_times,
                census_y,
                marker="|",
                s=45,
                color=palette[g],
                alpha=0.85,
                linewidth=0.75,
                label=None,
                zorder=3,
            )

    # x_max = shared_xmax(df, time_col, group_col)
    x_max = 1500
    x_pad = max((x_max if x_max is not None else 0) * 0.04, 1.0)
    ax.set_xlim(-x_pad, x_max if x_max is not None else None)
    ax.set_ylim(0.4, 1)
    ax.set_xlabel("")
    ax.tick_params(axis="x", labelbottom=False)
    ax.set_ylabel(y_label)
    handles, labels = ax.get_legend_handles_labels()
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.legend(
        handles,
        labels,
        title=legend_title,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.03),
        ncol=max(1, len(labels)),
        fontsize=9,
        title_fontsize=9,
        frameon=False,
        handlelength=1.8,
        columnspacing=1.5,
    )
    fig.text(0.5, 0.945, subtitle, ha="center", va="center", fontsize=10.5, color="#303030")

    # Risk set sizes panel: y-axis groups, x-axis days
    if x_max is None or x_max <= 0:
        tick_times = np.array([0.0])
    else:
        tick_times = np.linspace(0, x_max, 6)

    y_positions = np.arange(len(group_order))
    ax_risk.set_yticks(y_positions)
    ax_risk.set_yticklabels(group_order, fontsize=9)
    ax_risk.set_xlabel("days")
    ax_risk.set_xlim(-x_pad, x_max if x_max is not None else None)
    ax_risk.set_ylim(-0.7, len(group_order) - 0.3)
    ax_risk.set_xticks(tick_times)
    ax_risk.grid(axis="x", alpha=0.25)
    ax_risk.grid(axis="y", alpha=0.1)

    for i, g in enumerate(group_order):
        g_time = pd.to_numeric(df.loc[df[group_col] == g, time_col], errors="coerce").values
        g_time = g_time[np.isfinite(g_time)]
        for t in tick_times:
            n_risk = int(np.sum(g_time >= t))
            ax_risk.text(
                t,
                i,
                str(n_risk),
                ha="center",
                va="center",
                fontsize=9,
                color="#222222",
            )

    for spine in ["top", "right"]:
        ax_risk.spines[spine].set_visible(False)

    # Use constrained_layout or manual adjustment if tight_layout warns about GridSpec/RiskSet
    # Given the rect=[0, 0, 1, 0.93] to accommodate the top legend/subtitle
    # We can use subplots_adjust if tight_layout fails to handle the height_ratios correctly
    try:
        fig.tight_layout(rect=[0, 0, 1, 0.93])
    except:
        fig.subplots_adjust(top=0.9, bottom=0.1, left=0.12, right=0.95, hspace=0.1)

    fig.savefig(out_png)
    plt.close(fig)


def make_recur_plots(df, species, out_prefix, group_col="abundance_group", group_title="Abundance group", palette_override=None):
    uniq_groups = [g for g in pd.unique(df[group_col]) if pd.notna(g)]
    if len(uniq_groups) == 0:
        return

    if len(uniq_groups) == 2 and set(uniq_groups) == {"Low", "High"}:
        group_order = ["Low", "High"]
    else:
        group_order = sorted(uniq_groups)

    if palette_override is not None:
        palette = {g: palette_override.get(g, "#4C4C4C") for g in group_order}
    else:
        color_seq = ["#2166AC", "#B2182B", "#1B9E77", "#D95F02"]
        palette = {g: color_seq[i % len(color_seq)] for i, g in enumerate(group_order)}

    # Curve 1: cumulative recurrence rate over time (1 - KM of recurrence-free survival)
    rec_curve_dat = df.dropna(subset=[group_col, "PFS_time_days", "recur_event"]).copy()
    if rec_curve_dat[group_col].nunique() >= 2:
        p_logrank_rec = weighted_logrank_test(
            rec_curve_dat["PFS_time_days"].values,
            rec_curve_dat["recur_event"].values,
            rec_curve_dat[group_col].values,
            weight="logrank",
        )
        p_breslow_rec = weighted_logrank_test(
            rec_curve_dat["PFS_time_days"].values,
            rec_curve_dat["recur_event"].values,
            rec_curve_dat[group_col].values,
            weight="breslow",
        )
        p_tarone_rec = weighted_logrank_test(
            rec_curve_dat["PFS_time_days"].values,
            rec_curve_dat["recur_event"].values,
            rec_curve_dat[group_col].values,
            weight="tarone",
        )

        fig, ax = plt.subplots(figsize=(8.6, 6.2))
        for g in group_order:
            sub = rec_curve_dat[rec_curve_dat[group_col] == g]
            if sub.empty:
                continue
            x, y_surv = km_curve(sub["PFS_time_days"].values.astype(float), sub["recur_event"].values.astype(int))
            y_cum = 1 - y_surv
            ax.step(x, y_cum, where="post", label=g, color=palette[g], linewidth=2.4)

        rec_sub = (
            f"log-rank p={format_p(p_logrank_rec)} | "
            f"breslow p={format_p(p_breslow_rec)} | "
            f"tarone p={format_p(p_tarone_rec)}"
        )
        x_max = shared_xmax(rec_curve_dat, "PFS_time_days", group_col)
        ax.set_xlim(0, x_max if x_max is not None else None)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Time (days)")
        ax.set_ylabel("Cumulative recurrence rate")
        fig.suptitle(f"{species} cumulative recurrence curve", y=0.98, fontsize=16, fontweight="bold")
        ax.set_title(rec_sub, fontsize=11, pad=8, color="#303030")
        ax.legend(
            title=group_title,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0.0,
            fontsize=9,
            title_fontsize=9,
            frameon=False,
            handlelength=1.8,
        )
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
        fig.tight_layout(rect=[0, 0, 0.84, 0.95])
        fig.savefig(f"{out_prefix}_recurrence_cum_curve.png")
        plt.close(fig)

    rate_df = (
        df.dropna(subset=[group_col, "recur_event"]) 
        .groupby(group_col, as_index=False)
        .agg(n=("recur_event", "size"), recur_n=("recur_event", "sum"))
    )
    rate_df["recurrence_rate"] = rate_df["recur_n"] / rate_df["n"]

    p_chi = np.nan
    if len(rate_df[group_col].dropna().unique()) == 2:
        tab = pd.crosstab(df[group_col], df["recur_event"])
        try:
            p_chi = chi2_contingency(tab)[1]
        except Exception:
            p_chi = np.nan

    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    sns.barplot(data=rate_df, x=group_col, y="recurrence_rate", hue=group_col, palette=palette, legend=False, ax=ax, order=group_order)
    for _, r in rate_df.iterrows():
        g_name = r[group_col]
        x_pos = group_order.index(g_name) if g_name in group_order else 0
        ax.text(
            x=x_pos,
            y=r["recurrence_rate"] + 0.02,
            s=f"{int(r['recur_n'])}/{int(r['n'])}",
            ha="center",
            fontsize=11,
        )
    ax.set_ylim(0, 1)
    fig.suptitle(f"{species} recurrence rate by abundance group", y=0.98, fontsize=16, fontweight="bold")
    ax.set_title(f"Chi-square p = {format_p(p_chi)}", fontsize=11, pad=8, color="#303030")
    ax.set_xlabel(group_title)
    ax.set_ylabel("Recurrence rate")
    for spine in ["top", "right"]:
        ax.spines[spine].set_visible(False)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(f"{out_prefix}_recurrence_rate.png")
    plt.close(fig)

    tm = df[(df["recur_event"] == 1) & np.isfinite(df["recur_time_days"]) & (df["recur_time_days"] >= 0)].copy()
    if tm[group_col].nunique() >= 2 and len(tm) >= 4:
        g0, g1 = group_order[0], group_order[1]
        low_vals = tm.loc[tm[group_col] == g0, "recur_time_days"].values
        high_vals = tm.loc[tm[group_col] == g1, "recur_time_days"].values
        p_w = np.nan
        if len(low_vals) > 0 and len(high_vals) > 0:
            try:
                p_w = mannwhitneyu(low_vals, high_vals, alternative="two-sided")[1]
            except Exception:
                p_w = np.nan

        # Curve 2: ECDF of recurrence time among recurrent patients
        fig, ax = plt.subplots(figsize=(8.6, 6.2))
        for g in group_order:
            vals = np.sort(tm.loc[tm[group_col] == g, "recur_time_days"].values)
            if len(vals) == 0:
                continue
            y = np.arange(1, len(vals) + 1) / len(vals)
            ax.step(vals, y, where="post", color=palette[g], linewidth=2.4, label=g)

        x_max = shared_xmax(tm, "recur_time_days", group_col)
        ax.set_xlim(0, x_max if x_max is not None else None)
        ax.set_ylim(0, 1)
        ax.set_xlabel("Recurrence time (days)")
        ax.set_ylabel("Empirical cumulative proportion")
        fig.suptitle(f"{species} recurrence time curve", y=0.98, fontsize=16, fontweight="bold")
        ax.set_title(f"Wilcoxon p = {format_p(p_w)}", fontsize=11, pad=8, color="#303030")
        ax.legend(
            title=group_title,
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            borderaxespad=0.0,
            fontsize=9,
            title_fontsize=9,
            frameon=False,
            handlelength=1.8,
        )
        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)
        fig.tight_layout(rect=[0, 0, 0.84, 0.95])
        fig.savefig(f"{out_prefix}_recurrence_time_curve.png")
        plt.close(fig)


def main():
    species_list = [
        "Peptostreptococcus_stomatis",
        "Faecalibacterium_prausnitzii",
    ]

    abundance_file = Path("dataset/merged_dataset_relative.csv")
    survival_file = Path("dataset/survival/crc_os_pfs_full.csv")
    out_dir = Path("results/R_plots/single_species_survival")
    out_dir.mkdir(parents=True, exist_ok=True)

    ab = pd.read_csv(abundance_file)
    sv = pd.read_csv(survival_file)

    os_event_col = "OS_event" if "OS_event" in sv.columns else "死亡事件_OS(0/1)"
    pfs_event_col = "PFS_event" if "PFS_event" in sv.columns else "PFS事件(0/1)"

    if "SAMPLE_ID" not in ab.columns or "SAMPLE_ID" not in sv.columns:
        raise ValueError("SAMPLE_ID not found in abundance or survival data")

    joined = sv.merge(ab, on="SAMPLE_ID", how="inner")

    joined["recur_event"] = pd.to_numeric(joined["是否复发_原始"], errors="coerce")
    joined.loc[~joined["recur_event"].isin([0, 1]), "recur_event"] = np.nan

    enroll_date = parse_mixed_date(joined["入组日期"]) if "入组日期" in joined.columns else pd.Series(pd.NaT, index=joined.index)
    recur_date = parse_mixed_date(joined["复发日期_原始"]) if "复发日期_原始" in joined.columns else pd.Series(pd.NaT, index=joined.index)
    joined["recur_time_days"] = (recur_date - enroll_date).dt.days

    fallback = (joined["recur_event"] == 1) & (~np.isfinite(joined["recur_time_days"]) | (joined["recur_time_days"] < 0))
    joined.loc[fallback, "recur_time_days"] = joined.loc[fallback, "PFS_time_days"]

    stats_rows = []

    for sp in species_list:
        sp_col = f"tax_s__{sp}"
        if sp_col not in joined.columns:
            print(f"[WARN] Missing species column: {sp_col}")
            continue

        dat = joined[["SAMPLE_ID", "OS_time_days", os_event_col, "PFS_time_days", pfs_event_col, "recur_event", "recur_time_days", sp_col]].copy()
        dat = dat.rename(columns={sp_col: "abundance"})
        dat = dat.rename(columns={os_event_col: "OS_event", pfs_event_col: "PFS_event"})
        dat["abundance"] = pd.to_numeric(dat["abundance"], errors="coerce")

        med = float(np.nanmedian(dat["abundance"].values))
        dat["abundance_group"] = np.where(dat["abundance"] >= med, "High", "Low")

        p_logrank_os = weighted_logrank_test(dat["OS_time_days"].values, dat["OS_event"].values, dat["abundance_group"].values, weight="logrank")
        p_breslow_os = weighted_logrank_test(dat["OS_time_days"].values, dat["OS_event"].values, dat["abundance_group"].values, weight="breslow")
        p_tarone_os = weighted_logrank_test(dat["OS_time_days"].values, dat["OS_event"].values, dat["abundance_group"].values, weight="tarone")

        p_logrank_pfs = weighted_logrank_test(dat["PFS_time_days"].values, dat["PFS_event"].values, dat["abundance_group"].values, weight="logrank")
        p_breslow_pfs = weighted_logrank_test(dat["PFS_time_days"].values, dat["PFS_event"].values, dat["abundance_group"].values, weight="breslow")
        p_tarone_pfs = weighted_logrank_test(dat["PFS_time_days"].values, dat["PFS_event"].values, dat["abundance_group"].values, weight="tarone")

        stats_rows.append({
            "species": sp,
            "median_abundance": med,
            "n": int(len(dat)),
            "n_low": int((dat["abundance_group"] == "Low").sum()),
            "n_high": int((dat["abundance_group"] == "High").sum()),
            "p_logrank_os": p_logrank_os,
            "p_breslow_os": p_breslow_os,
            "p_tarone_os": p_tarone_os,
            "p_logrank_pfs": p_logrank_pfs,
            "p_breslow_pfs": p_breslow_pfs,
            "p_tarone_pfs": p_tarone_pfs,
        })

        os_sub = f"log-rank p={p_logrank_os:.3g} | breslow p={p_breslow_os:.3g} | tarone p={p_tarone_os:.3g}"
        pfs_sub = f"log-rank p={p_logrank_pfs:.3g} | breslow p={p_breslow_pfs:.3g} | tarone p={p_tarone_pfs:.3g}"

        prefix = out_dir / sp
        make_km_plot(
            dat,
            "OS_time_days",
            "OS_event",
            "abundance_group",
            f"{sp} - Overall Survival (KM)",
            os_sub,
            f"{prefix}_KM_OS.png",
            legend_title="Abundance group",
            y_label="OS survival probability",
        )
        make_km_plot(
            dat,
            "PFS_time_days",
            "PFS_event",
            "abundance_group",
            f"{sp} - Progression-Free Survival (KM)",
            pfs_sub,
            f"{prefix}_KM_PFS.png",
            legend_title="Abundance group",
            y_label="PFS survival probability",
        )
        make_recur_plots(dat, sp, str(prefix))

        dat.to_csv(f"{prefix}_patient_level.csv", index=False)

    if stats_rows:
        pd.DataFrame(stats_rows).to_csv(out_dir / "single_species_survival_stats.csv", index=False)

    # Additional analysis: impact of differentiation on survival and recurrence
    if "differentiation" in joined.columns:
        diff_dat = joined[["SAMPLE_ID", "OS_time_days", os_event_col, "PFS_time_days", pfs_event_col, "recur_event", "recur_time_days", "differentiation"]].copy()
        diff_dat = diff_dat.rename(columns={os_event_col: "OS_event", pfs_event_col: "PFS_event"})
        diff_dat["differentiation"] = pd.to_numeric(diff_dat["differentiation"], errors="coerce")
        diff_dat = diff_dat[diff_dat["differentiation"].isin([0, 1])].copy()
        diff_dat["diff_group"] = np.where(diff_dat["differentiation"] == 1, "CRC_poordiff", "CRC_welldiff")

        p_logrank_os_d = weighted_logrank_test(diff_dat["OS_time_days"].values, diff_dat["OS_event"].values, diff_dat["diff_group"].values, weight="logrank")
        p_breslow_os_d = weighted_logrank_test(diff_dat["OS_time_days"].values, diff_dat["OS_event"].values, diff_dat["diff_group"].values, weight="breslow")
        p_tarone_os_d = weighted_logrank_test(diff_dat["OS_time_days"].values, diff_dat["OS_event"].values, diff_dat["diff_group"].values, weight="tarone")

        p_logrank_pfs_d = weighted_logrank_test(diff_dat["PFS_time_days"].values, diff_dat["PFS_event"].values, diff_dat["diff_group"].values, weight="logrank")
        p_breslow_pfs_d = weighted_logrank_test(diff_dat["PFS_time_days"].values, diff_dat["PFS_event"].values, diff_dat["diff_group"].values, weight="breslow")
        p_tarone_pfs_d = weighted_logrank_test(diff_dat["PFS_time_days"].values, diff_dat["PFS_event"].values, diff_dat["diff_group"].values, weight="tarone")

        os_sub_d = (
            f"log-rank p={format_p(p_logrank_os_d)} | "
            f"breslow p={format_p(p_breslow_os_d)} | "
            f"tarone p={format_p(p_tarone_os_d)}"
        )
        pfs_sub_d = (
            f"log-rank p={format_p(p_logrank_pfs_d)} | "
            f"breslow p={format_p(p_breslow_pfs_d)} | "
            f"tarone p={format_p(p_tarone_pfs_d)}"
        )

        make_km_plot(
            diff_dat,
            "OS_time_days",
            "OS_event",
            "diff_group",
            "Differentiation impact on Overall Survival (KM)",
            os_sub_d,
            str(out_dir / "differentiation_KM_OS.png"),
            palette_override={"CRC_poordiff": "#FF6B6B", "CRC_welldiff": "#FFD166"},
            legend_title="Differentiation group",
            y_label="OS survival probability",
        )
        make_km_plot(
            diff_dat,
            "PFS_time_days",
            "PFS_event",
            "diff_group",
            "Differentiation impact on Progression-Free Survival (KM)",
            pfs_sub_d,
            str(out_dir / "differentiation_KM_PFS.png"),
            palette_override={"CRC_poordiff": "#FF6B6B", "CRC_welldiff": "#FFD166"},
            legend_title="Differentiation group",
            y_label="PFS survival probability",
        )

        make_recur_plots(
            diff_dat,
            "Differentiation",
            str(out_dir / "differentiation"),
            group_col="diff_group",
            group_title="Differentiation group",
            palette_override={"CRC_poordiff": "#FF6B6B", "CRC_welldiff": "#FFD166"},
        )

        diff_stats = pd.DataFrame([
            {
                "analysis": "differentiation",
                "n": int(len(diff_dat)),
                "n_well_moderate": int((diff_dat["diff_group"] == "CRC_welldiff").sum()),
                "n_poor": int((diff_dat["diff_group"] == "CRC_poordiff").sum()),
                "p_logrank_os": p_logrank_os_d,
                "p_breslow_os": p_breslow_os_d,
                "p_tarone_os": p_tarone_os_d,
                "p_logrank_pfs": p_logrank_pfs_d,
                "p_breslow_pfs": p_breslow_pfs_d,
                "p_tarone_pfs": p_tarone_pfs_d,
            }
        ])
        diff_stats.to_csv(out_dir / "differentiation_survival_recurrence_stats.csv", index=False)

    print(f"Done. Outputs in: {out_dir}")


if __name__ == "__main__":
    main()
