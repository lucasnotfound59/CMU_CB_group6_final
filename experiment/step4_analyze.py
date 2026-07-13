#!/usr/bin/env python3
"""
Step 4: Analyze correlation between B-factors and docking RMSD.
================================================================
Usage:
    python step4_analyze.py

This script:
    1. Load B-factor data (step2) and cross-docking results (step3)
    2. For each cross-docking pair, identify which binding-site residues
       have the highest B-factors in the RECEPTOR structure
    3. Test: do pairs with higher max/mean B-factor residues show higher RMSD?
    4. Statistical tests: Pearson/Spearman correlation
    5. Generate summary tables and figures

Output:
    output/analysis_summary.csv      — per-pair analysis with B-factor metrics
    output/statistical_tests.txt     — correlation results
    figures/bfactor_vs_rmsd.png      — scatter plot
    figures/bfactor_distribution.png — B-factor distribution by success/failure
    figures/heatmap_bfactor_rmsd.png — heatmap of B-factor vs RMSD per residue
"""

import csv
import os
import sys

import numpy as np
from scipy import stats

sys.path.insert(0, os.path.dirname(__file__))
from config import OUTPUT_DIR, FIGURES_DIR, RMSD_SUCCESS_THRESHOLD, BFACTOR_HIGH_PERCENTILE

try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False
    print("Warning: matplotlib not available, skipping figures")


def load_csv(path: str) -> list[dict]:
    with open(path) as f:
        return list(csv.DictReader(f))


def analyze():
    # === Load data ===
    bfactor_path = os.path.join(OUTPUT_DIR, "bfactor_summary.csv")
    structure_path = os.path.join(OUTPUT_DIR, "bfactor_per_structure.csv")
    docking_path = os.path.join(OUTPUT_DIR, "cross_docking_results.csv")

    if not os.path.exists(bfactor_path):
        print("ERROR: Run step2 first (bfactor_summary.csv not found)")
        sys.exit(1)
    if not os.path.exists(docking_path):
        print("ERROR: Run step3 first (cross_docking_results.csv not found)")
        sys.exit(1)

    bfactors = load_csv(bfactor_path)
    structures = load_csv(structure_path) if os.path.exists(structure_path) else []
    docking = load_csv(docking_path)

    print(f"Loaded {len(bfactors)} residue records, {len(docking)} docking pairs\n")

    # === Build B-factor lookup ===
    # pdb_id -> list of residue dicts
    bf_by_pdb = {}
    for r in bfactors:
        pdb = r["pdb_id"]
        if pdb not in bf_by_pdb:
            bf_by_pdb[pdb] = []
        bf_by_pdb[pdb].append(r)

    native_ligand_by_pdb = {
        r["pdb_id"]: r.get("ligand", "")
        for r in structures
    }

    # === Merge: for each docking pair, get receptor B-factor stats ===
    analysis = []
    for pair in docking:
        rec = pair["receptor"]
        rmsd = pair["rmsd"]
        if rmsd is None or rmsd == "":
            continue
        rmsd = float(rmsd)

        rec_bfactors = bf_by_pdb.get(rec, [])
        if not rec_bfactors:
            continue

        # Receptor B-factor statistics
        avg_bf = [float(r["avg_bfactor"]) for r in rec_bfactors]
        rel_bf = [float(r["relative_bfactor"]) for r in rec_bfactors
                  if r["relative_bfactor"] is not None and r["relative_bfactor"] != ""]

        # Identify high B-factor residues
        if avg_bf:
            threshold = np.percentile(avg_bf, BFACTOR_HIGH_PERCENTILE)
            high_bf_residues = [r for r in rec_bfactors
                                if float(r["avg_bfactor"]) >= threshold]
        else:
            high_bf_residues = []

        analysis.append({
            "receptor": rec,
            "receptor_native_ligand": native_ligand_by_pdb.get(rec, ""),
            "ligand_from": pair["ligand_from"],
            "docked_ligand_native_ligand": native_ligand_by_pdb.get(pair["ligand_from"], ""),
            "rmsd": rmsd,
            "affinity": pair.get("affinity"),
            "receptor_pocket_avg_bfactor": round(np.mean(avg_bf), 2) if avg_bf else None,
            "receptor_pocket_max_bfactor": round(max(avg_bf), 2) if avg_bf else None,
            "receptor_pocket_std_bfactor": round(np.std(avg_bf), 2) if avg_bf else None,
            "receptor_mean_relative_bfactor": round(np.mean(rel_bf), 3) if rel_bf else None,
            "n_high_bfactor_residues": len(high_bf_residues),
            "high_bfactor_residues": ";".join(
                f"{r['resname']}{r['resseq']}" for r in high_bf_residues
            ),
            "success": rmsd < RMSD_SUCCESS_THRESHOLD,
        })

    if not analysis:
        print("ERROR: No valid data pairs for analysis")
        sys.exit(1)

    # === Statistical tests ===
    rmsds = [a["rmsd"] for a in analysis]
    pocket_avgs = [a["receptor_pocket_avg_bfactor"] for a in analysis]
    pocket_maxs = [a["receptor_pocket_max_bfactor"] for a in analysis]
    pocket_stds = [a["receptor_pocket_std_bfactor"] for a in analysis]
    rel_bfs = [a["receptor_mean_relative_bfactor"] for a in analysis
               if a["receptor_mean_relative_bfactor"] is not None]

    print("=" * 60)
    print("STATISTICAL ANALYSIS")
    print("=" * 60)

    # Correlation: pocket avg B-factor vs RMSD
    if len(pocket_avgs) == len(rmsds):
        r_pearson, p_pearson = stats.pearsonr(pocket_avgs, rmsds)
        r_spearman, p_spearman = stats.spearmanr(pocket_avgs, rmsds)
        print(f"\nPocket avg B-factor vs RMSD:")
        print(f"  Pearson:  r = {r_pearson:.4f}, p = {p_pearson:.4f}")
        print(f"  Spearman: ρ = {r_spearman:.4f}, p = {p_spearman:.4f}")

    # Correlation: pocket max B-factor vs RMSD
    if len(pocket_maxs) == len(rmsds):
        r_max, p_max = stats.pearsonr(pocket_maxs, rmsds)
        print(f"\nPocket max B-factor vs RMSD:")
        print(f"  Pearson:  r = {r_max:.4f}, p = {p_max:.4f}")

    # Correlation: pocket std B-factor vs RMSD
    if len(pocket_stds) == len(rmsds):
        r_std, p_std = stats.pearsonr(pocket_stds, rmsds)
        print(f"\nPocket std B-factor vs RMSD:")
        print(f"  Pearson:  r = {r_std:.4f}, p = {p_std:.4f}")

    # Success rate comparison
    success_rmsds = [a["rmsd"] for a in analysis if a["success"]]
    fail_rmsds = [a["rmsd"] for a in analysis if not a["success"]]
    print(f"\nCross-docking success rate: {len(success_rmsds)}/{len(analysis)} "
          f"({100*len(success_rmsds)/len(analysis):.1f}%)")
    native_rows = [a for a in analysis if _is_native_interaction(a)]
    if native_rows:
        print(f"Native/self interactions available for labeling: {len(native_rows)}")
    else:
        print("Native/self interactions available for labeling: 0 "
              "(cross-docking excludes receptor == ligand_from)")

    # T-test: B-factor of successful vs failed pairs
    success_bf = [a["receptor_pocket_avg_bfactor"] for a in analysis if a["success"]]
    fail_bf = [a["receptor_pocket_avg_bfactor"] for a in analysis if not a["success"]]
    t_stat, t_pval = None, None
    if len(success_bf) > 1 and len(fail_bf) > 1:
        t_stat, t_pval = stats.ttest_ind(success_bf, fail_bf)
        print(f"\nT-test: B-factor (success vs failure)")
        print(f"  Success: mean={np.mean(success_bf):.2f}, n={len(success_bf)}")
        print(f"  Failure: mean={np.mean(fail_bf):.2f}, n={len(fail_bf)}")
        print(f"  t = {t_stat:.4f}, p = {t_pval:.4f}")

    # === Save analysis CSV ===
    csv_path = os.path.join(OUTPUT_DIR, "analysis_summary.csv")
    fieldnames = ["receptor", "receptor_native_ligand", "ligand_from",
                  "docked_ligand_native_ligand", "rmsd", "affinity",
                  "receptor_pocket_avg_bfactor", "receptor_pocket_max_bfactor",
                  "receptor_pocket_std_bfactor", "receptor_mean_relative_bfactor",
                  "n_high_bfactor_residues", "high_bfactor_residues", "success"]
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(analysis)
    print(f"\nSaved: {csv_path}")

    # === Save statistical summary ===
    txt_path = os.path.join(OUTPUT_DIR, "statistical_tests.txt")
    with open(txt_path, "w") as f:
        f.write("B-factor vs Docking RMSD: Statistical Analysis\n")
        f.write("=" * 50 + "\n\n")
        f.write(f"N pairs: {len(analysis)}\n")
        f.write(f"Success rate (RMSD < {RMSD_SUCCESS_THRESHOLD}Å): "
                f"{len(success_rmsds)}/{len(analysis)}\n\n")
        f.write(f"Pocket avg B-factor vs RMSD:\n")
        f.write(f"  Pearson r = {r_pearson:.4f}, p = {p_pearson:.4f}\n")
        f.write(f"  Spearman ρ = {r_spearman:.4f}, p = {p_spearman:.4f}\n\n")
        if t_stat is not None:
            f.write(f"T-test (success vs failure B-factors):\n")
            f.write(f"  t = {t_stat:.4f}, p = {t_pval:.4f}\n")
        else:
            f.write(f"T-test: skipped (insufficient data in one or both groups)\n")
    print(f"Saved: {txt_path}")

    # === Generate figures ===
    if HAS_MPL:
        generate_figures(analysis, rmsds, pocket_avgs, success_bf, fail_bf)


def generate_figures(analysis, rmsds, pocket_avgs, success_bf, fail_bf):
    """Generate analysis figures."""
    ligands = [a["ligand_from"] for a in analysis]
    unique_ligands = sorted(set(ligands))
    ligand_colors = _make_category_colors(unique_ligands)

    # --- Figure 1: B-factor vs RMSD scatter ---
    fig, ax = plt.subplots(figsize=(8, 6))
    colors = ["green" if a["success"] else "red" for a in analysis]
    ax.scatter(pocket_avgs, rmsds, c=colors, alpha=0.6, s=30, edgecolors="gray", linewidth=0.5)
    ax.set_xlabel("Receptor Pocket Average B-factor (Å²)", fontsize=12)
    ax.set_ylabel("Cross-docking RMSD (Å)", fontsize=12)
    ax.set_title("B-factor vs Docking RMSD", fontsize=14)
    ax.axhline(y=RMSD_SUCCESS_THRESHOLD, color="gray", linestyle="--", alpha=0.5, label="RMSD = 2Å")
    ax.legend()

    # Trend line
    if len(pocket_avgs) > 2:
        z = np.polyfit(pocket_avgs, rmsds, 1)
        p = np.poly1d(z)
        x_line = np.linspace(min(pocket_avgs), max(pocket_avgs), 100)
        ax.plot(x_line, p(x_line), "b--", alpha=0.5, label="Linear fit")
        ax.legend()

    fig_path = os.path.join(FIGURES_DIR, "bfactor_vs_rmsd.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    # --- Figure 1b: B-factor vs RMSD, colored by original ligand structure ---
    fig, ax = plt.subplots(figsize=(10, 7))
    for ligand in unique_ligands:
        ligand_rows = [a for a in analysis if a["ligand_from"] == ligand]
        x = [a["receptor_pocket_avg_bfactor"] for a in ligand_rows]
        y = [a["rmsd"] for a in ligand_rows]
        ax.scatter(
            x, y,
            color=ligand_colors[ligand],
            alpha=0.75,
            s=32,
            edgecolors="gray",
            linewidth=0.4,
            label=ligand,
        )
        for a in ligand_rows:
            if _is_native_interaction(a):
                ax.annotate(
                    _native_ligand_label(a),
                    (
                        a["receptor_pocket_avg_bfactor"],
                        a["rmsd"],
                    ),
                    xytext=(4, 3),
                    textcoords="offset points",
                    fontsize=8,
                    fontweight="bold",
                    alpha=0.9,
                    color="black",
                )

    ax.set_xlabel("Receptor Pocket Average B-factor (Å²)", fontsize=12)
    ax.set_ylabel("Cross-docking RMSD (Å)", fontsize=12)
    ax.set_title("B-factor vs Docking RMSD by Original Ligand Structure", fontsize=14)
    ax.axhline(
        y=RMSD_SUCCESS_THRESHOLD,
        color="gray",
        linestyle="--",
        alpha=0.6,
        label=f"RMSD = {RMSD_SUCCESS_THRESHOLD:g}Å",
    )
    _add_ligand_legend(ax, unique_ligands)

    fig_path = os.path.join(FIGURES_DIR, "bfactor_vs_rmsd_by_ligand.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    # --- Figure 1c: Docking score / affinity vs B-factor ---
    score_rows = [
        a for a in analysis
        if a.get("affinity") not in (None, "")
        and a.get("receptor_pocket_avg_bfactor") is not None
    ]
    if score_rows:
        fig, ax = plt.subplots(figsize=(10, 7))
        for ligand in unique_ligands:
            ligand_rows = [a for a in score_rows if a["ligand_from"] == ligand]
            if not ligand_rows:
                continue
            x = [a["receptor_pocket_avg_bfactor"] for a in ligand_rows]
            y = [float(a["affinity"]) for a in ligand_rows]
            ax.scatter(
                x, y,
                color=ligand_colors[ligand],
                alpha=0.75,
                s=32,
                edgecolors="gray",
                linewidth=0.4,
                label=ligand,
            )
            for a in ligand_rows:
                if _is_native_interaction(a):
                    ax.annotate(
                        _native_ligand_label(a),
                        (
                            a["receptor_pocket_avg_bfactor"],
                            float(a["affinity"]),
                        ),
                        xytext=(4, 3),
                        textcoords="offset points",
                        fontsize=8,
                        fontweight="bold",
                        alpha=0.9,
                        color="black",
                    )

        score_bfactors = [a["receptor_pocket_avg_bfactor"] for a in score_rows]
        scores = [float(a["affinity"]) for a in score_rows]
        if len(score_bfactors) > 2 and len(set(score_bfactors)) > 1:
            z = np.polyfit(score_bfactors, scores, 1)
            p = np.poly1d(z)
            x_line = np.linspace(min(score_bfactors), max(score_bfactors), 100)
            ax.plot(x_line, p(x_line), "k--", alpha=0.5, label="Linear fit")

        ax.set_xlabel("Receptor Pocket Average B-factor (Å²)", fontsize=12)
        ax.set_ylabel("Vina Affinity / Docking Score (kcal/mol)", fontsize=12)
        ax.set_title("Docking Score vs B-factor by Original Ligand Structure", fontsize=14)
        _add_ligand_legend(ax, unique_ligands)

        fig_path = os.path.join(FIGURES_DIR, "bfactor_vs_docking_score_by_ligand.png")
        fig.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.close(fig)
        print(f"Saved: {fig_path}")
    else:
        print("Skipped: no affinity values available for score-vs-B-factor figure")

    # --- Figure 2: B-factor distribution (success vs failure) ---
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist([success_bf, fail_bf], bins=20, label=["Success (RMSD<2Å)", "Failure (RMSD≥2Å)"],
            color=["green", "red"], alpha=0.7, edgecolor="black", linewidth=0.5)
    ax.set_xlabel("Receptor Pocket Average B-factor (Å²)", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("B-factor Distribution: Successful vs Failed Cross-Docking", fontsize=13)
    ax.legend()
    fig_path = os.path.join(FIGURES_DIR, "bfactor_distribution.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    # --- Figure 3: RMSD distribution ---
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.hist(rmsds, bins=30, color="steelblue", edgecolor="black", linewidth=0.5)
    ax.axvline(x=RMSD_SUCCESS_THRESHOLD, color="red", linestyle="--", label=f"RMSD = {RMSD_SUCCESS_THRESHOLD}Å")
    ax.set_xlabel("RMSD (Å)", fontsize=12)
    ax.set_ylabel("Count", fontsize=12)
    ax.set_title("Cross-Docking RMSD Distribution", fontsize=14)
    ax.legend()
    fig_path = os.path.join(FIGURES_DIR, "rmsd_distribution.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    # --- Figure 4: B-factor vs RMSD spread per receptor ---
    # Hypothesis: higher B-factor → wider RMSD variance across ligands
    from collections import defaultdict
    rec_rmsds = defaultdict(list)
    rec_bfactor = {}
    for a in analysis:
        rec = a["receptor"]
        rec_rmsds[rec].append(a["rmsd"])
        if rec not in rec_bfactor:
            rec_bfactor[rec] = a["receptor_pocket_avg_bfactor"]

    rec_list = sorted(rec_bfactor, key=lambda r: rec_bfactor[r])
    bf_vals = [rec_bfactor[r] for r in rec_list]
    rmsd_stds = [np.std(rec_rmsds[r]) for r in rec_list]
    rmsd_means = [np.mean(rec_rmsds[r]) for r in rec_list]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # 4a: B-factor vs RMSD std
    ax1.scatter(bf_vals, rmsd_stds, c="coral", s=50, edgecolors="gray", linewidth=0.5, alpha=0.8)
    for i, r in enumerate(rec_list):
        ax1.annotate(r, (bf_vals[i], rmsd_stds[i]), fontsize=5, alpha=0.6,
                     xytext=(3, 2), textcoords="offset points")
    if len(bf_vals) > 2:
        z = np.polyfit(bf_vals, rmsd_stds, 1)
        p = np.poly1d(z)
        x_line = np.linspace(min(bf_vals), max(bf_vals), 100)
        ax1.plot(x_line, p(x_line), "k--", alpha=0.5, label="Linear fit")
        r_corr, p_corr = stats.spearmanr(bf_vals, rmsd_stds)
        ax1.set_title(f"B-factor vs RMSD Spread (Spearman ρ={r_corr:.3f}, p={p_corr:.3f})", fontsize=12)
        ax1.legend()
    else:
        ax1.set_title("B-factor vs RMSD Spread", fontsize=12)
    ax1.set_xlabel("Receptor Pocket Avg B-factor (Å²)", fontsize=11)
    ax1.set_ylabel("RMSD Std across all ligands (Å)", fontsize=11)

    # 4b: B-factor vs RMSD mean
    ax2.scatter(bf_vals, rmsd_means, c="steelblue", s=50, edgecolors="gray", linewidth=0.5, alpha=0.8)
    for i, r in enumerate(rec_list):
        ax2.annotate(r, (bf_vals[i], rmsd_means[i]), fontsize=5, alpha=0.6,
                     xytext=(3, 2), textcoords="offset points")
    if len(bf_vals) > 2:
        z2 = np.polyfit(bf_vals, rmsd_means, 1)
        p2 = np.poly1d(z2)
        ax2.plot(x_line, p2(x_line), "k--", alpha=0.5, label="Linear fit")
        r_corr2, p_corr2 = stats.spearmanr(bf_vals, rmsd_means)
        ax2.set_title(f"B-factor vs Mean RMSD (Spearman ρ={r_corr2:.3f}, p={p_corr2:.3f})", fontsize=12)
        ax2.legend()
    else:
        ax2.set_title("B-factor vs Mean RMSD", fontsize=12)
    ax2.set_xlabel("Receptor Pocket Avg B-factor (Å²)", fontsize=11)
    ax2.set_ylabel("Mean RMSD across all ligands (Å)", fontsize=11)

    fig.tight_layout()
    fig_path = os.path.join(FIGURES_DIR, "bfactor_vs_rmsd_spread.png")
    fig.savefig(fig_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {fig_path}")

    # --- Figure 5: Self-docking vs Cross-docking RMSD per receptor ---
    self_rows = [a for a in analysis if a["receptor"] == a["ligand_from"]]
    cross_rows = [a for a in analysis if a["receptor"] != a["ligand_from"]]
    if self_rows and cross_rows:
        # Group cross-docking RMSD by receptor
        cross_by_rec = defaultdict(list)
        for a in cross_rows:
            cross_by_rec[a["receptor"]].append(a["rmsd"])
        self_by_rec = {a["receptor"]: a["rmsd"] for a in self_rows}

        common_recs = sorted(set(self_by_rec) & set(cross_by_rec))
        if common_recs:
            fig, ax = plt.subplots(figsize=(12, 6))
            x_pos = np.arange(len(common_recs))

            self_vals = [self_by_rec[r] for r in common_recs]
            cross_means = [np.mean(cross_by_rec[r]) for r in common_recs]
            cross_stds = [np.std(cross_by_rec[r]) for r in common_recs]

            ax.bar(x_pos - 0.2, self_vals, 0.4, color="#2196F3", alpha=0.8, label="Self-docking RMSD")
            ax.bar(x_pos + 0.2, cross_means, 0.4, yerr=cross_stds, capsize=3,
                   color="#FF9800", alpha=0.8, label="Cross-docking mean ± std")

            ax.set_xticks(x_pos)
            ax.set_xticklabels(common_recs, rotation=90, fontsize=6)
            ax.axhline(y=RMSD_SUCCESS_THRESHOLD, color="red", linestyle="--", alpha=0.5,
                       label=f"RMSD = {RMSD_SUCCESS_THRESHOLD}Å")
            ax.set_xlabel("Receptor (PDB ID)", fontsize=11)
            ax.set_ylabel("RMSD (Å)", fontsize=11)
            ax.set_title("Self-Docking vs Cross-Docking RMSD per Receptor", fontsize=13)
            ax.legend(fontsize=9)

            fig.tight_layout()
            fig_path = os.path.join(FIGURES_DIR, "self_vs_cross_docking.png")
            fig.savefig(fig_path, dpi=150, bbox_inches="tight")
            plt.close(fig)
            print(f"Saved: {fig_path}")

            # Print summary stats
            self_success = sum(1 for r in common_recs if self_by_rec[r] < RMSD_SUCCESS_THRESHOLD)
            cross_success = sum(1 for r in common_recs
                                if np.mean(cross_by_rec[r]) < RMSD_SUCCESS_THRESHOLD)
            print(f"  Self-docking success: {self_success}/{len(common_recs)}")
            print(f"  Cross-docking mean success: {cross_success}/{len(common_recs)}")
    else:
        print("Skipped self-vs-cross figure: no self-docking data found")


def _make_category_colors(labels):
    """Return stable matplotlib colors for categorical ligand labels."""
    if not labels:
        return {}

    cmap_name = "tab20" if len(labels) <= 20 else "nipy_spectral"
    cmap = plt.get_cmap(cmap_name, max(len(labels), 1))
    return {label: cmap(i) for i, label in enumerate(labels)}


def _add_ligand_legend(ax, labels):
    """Add a compact outside legend when ligand count is reasonable."""
    if len(labels) <= 20:
        ncol = 1 if len(labels) <= 12 else 2
        ax.legend(
            title="Original ligand structure",
            bbox_to_anchor=(1.02, 1),
            loc="upper left",
            borderaxespad=0,
            fontsize=8,
            ncol=ncol,
        )
    else:
        # Too many ligands for a readable legend; color still encodes ligand identity.
        ax.text(
            0.99,
            0.01,
            f"{len(labels)} ligand structures colored",
            transform=ax.transAxes,
            ha="right",
            va="bottom",
            fontsize=9,
            color="dimgray",
        )


def _native_ligand_label(row):
    """Label a point with the receptor structure and its own co-crystal ligand."""
    ligand = row.get("receptor_native_ligand") or "native?"
    return f"{row['receptor']}:{ligand}"


def _is_native_interaction(row):
    """True when the docked ligand comes from the same structure as the receptor."""
    return row.get("receptor", "").upper() == row.get("ligand_from", "").upper()


if __name__ == "__main__":
    analyze()
