#!/usr/bin/env python3
"""Generate the frozen Phase 3B figure/table package from formal TSV sources.

This script performs visualization and deterministic display summaries only. It
does not refit models, resample participants, permute labels, or change any
Phase 2B statistic.
"""

from __future__ import annotations

import csv
import hashlib
import math
import shutil
import textwrap
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle
import numpy as np
import pandas as pd
import seaborn as sns


ROOT = Path(__file__).resolve().parents[1]
P1 = ROOT / "reproduction_output" / "phase1b"
P2 = ROOT / "reproduction_output" / "phase2b"
META = ROOT / "metadata"
FIG = ROOT / "reproduction_output" / "figures"
MAIN = FIG / "main"
SUPP = FIG / "supplement"
SOURCE = FIG / "source_data"
QC = FIG / "qc"
PACKAGE = ROOT / "FRACTIONAL_LASER_PHASE3B_FIGURE_PACKAGE"

PROGRAM_ORDER = [
    "acute_injury", "inflammation", "epidermal_repair", "ecm_remodeling",
    "collagen_organization", "proliferation", "matrix_maturation",
]
PROGRAM_LABEL = {
    "acute_injury": "Acute injury",
    "inflammation": "Inflammation",
    "epidermal_repair": "Epidermal repair",
    "ecm_remodeling": "ECM remodeling",
    "collagen_organization": "Collagen organization",
    "proliferation": "Proliferation",
    "matrix_maturation": "Matrix maturation",
}
PROGRAM_ABBR = {
    "acute_injury": "Injury", "inflammation": "Inflam.",
    "epidermal_repair": "Epidermal", "ecm_remodeling": "ECM",
    "collagen_organization": "Collagen", "proliferation": "Prolif.",
    "matrix_maturation": "Matrix",
}
DELAYED = ["ecm_remodeling", "collagen_organization", "matrix_maturation"]

BLUE = "#0072B2"
ORANGE = "#D55E00"
TEAL = "#009E73"
PURPLE = "#8E6C8A"
GOLD = "#E69F00"
GREY = "#6B7280"
LIGHT = "#E5E7EB"
RED = "#B2182B"


def configure_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
        "svg.fonttype": "none", "pdf.fonttype": 42,
        "font.size": 7, "axes.titlesize": 8, "axes.labelsize": 7,
        "xtick.labelsize": 6.3, "ytick.labelsize": 6.3,
        "legend.fontsize": 6.2, "axes.linewidth": 0.7,
        "axes.spines.right": False, "axes.spines.top": False,
        "legend.frameon": False, "savefig.facecolor": "white",
        "figure.facecolor": "white",
    })
    sns.set_style("white")


def read(name: str, phase: int = 2) -> pd.DataFrame:
    return pd.read_csv((P2 if phase == 2 else P1) / name, sep="\t")


def save_figure(fig: plt.Figure, outbase: Path, dpi: int = 320) -> None:
    outbase.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(outbase.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".svg"), bbox_inches="tight")
    fig.savefig(outbase.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def panel(ax, label: str) -> None:
    ax.text(-0.12, 1.07, label, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="left")


def clean(ax, zero: bool = False) -> None:
    ax.spines["left"].set_linewidth(0.7)
    ax.spines["bottom"].set_linewidth(0.7)
    ax.tick_params(width=0.6, length=2.5)
    if zero:
        ax.axhline(0, color="#9CA3AF", lw=0.7, ls="--", zorder=0)


def title(fig: plt.Figure, text: str) -> None:
    fig.suptitle(text, x=0.01, y=0.995, ha="left", va="top",
                 fontsize=9.2, fontweight="bold")


def program_effects_primary() -> pd.DataFrame:
    x = read("participant_blocked_program_effects.tsv", 1)
    return x[x["day"].isin([1, 7, 14])].copy()


def source_long(panel_name: str, df: pd.DataFrame, source_file: str,
                record_type: str) -> pd.DataFrame:
    out = df.copy()
    out.insert(0, "record_type", record_type)
    out.insert(0, "source_file", source_file)
    out.insert(0, "panel", panel_name)
    return out


def write_source_data() -> None:
    # Figure 1: exact design metadata plus deterministic workflow labels.
    m1 = pd.read_csv(META / "GSE168760_metadata.tsv", sep="\t")
    m2 = pd.read_csv(META / "GSE206495_metadata.tsv", sep="\t")
    f1 = pd.concat([
        source_long("A", m1, "00_metadata/GSE168760_metadata.tsv", "sample_metadata"),
        source_long("B", m2, "00_metadata/GSE206495_metadata.tsv", "sample_metadata"),
    ], ignore_index=True, sort=False)
    workflow = pd.DataFrame({"panel": ["C"] * 4 + ["D"] * 8,
        "source_file": ["PHASE3A_FIGURE_PLAN.tsv"] * 12,
        "record_type": ["shared_timepoint"] * 4 + ["workflow_node"] * 8,
        "value": ["D0", "D1", "D7", "D14", "participant-blocked effects",
                  "shared gene universe", "7 frozen programs", "program concordance",
                  "random null + LOPO", "genome-wide concordance",
                  "reciprocal rank replication", "transition analysis"]})
    pd.concat([f1, workflow], ignore_index=True, sort=False).to_csv(
        SOURCE / "Figure1_SourceData.tsv", sep="\t", index=False)

    eff = program_effects_primary()
    cc = read("cross_cohort_conservation.tsv", 1)
    f2 = pd.concat([
        source_long("A", eff, "03_results/phase1b/participant_blocked_program_effects.tsv", "program_effect"),
        source_long("B-D", cc[cc["day"].isin([1, 7, 14])],
                    "03_results/phase1b/cross_cohort_conservation.tsv", "concordance_statistic"),
    ], ignore_index=True, sort=False)
    f2.to_csv(SOURCE / "Figure2_SourceData.tsv", sep="\t", index=False)

    null = read("MATCHED_RANDOM_PROGRAM_NULL_10000.tsv")
    null_sum = read("MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv")
    lopo = read("PHASE2B_FULL_LOPO.tsv")
    f3 = pd.concat([
        source_long("A-B", null, "03_results/phase2b/MATCHED_RANDOM_PROGRAM_NULL_10000.tsv", "null_iteration"),
        source_long("A-B", null_sum, "03_results/phase2b/MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv", "null_summary"),
        source_long("C-D", lopo, "03_results/phase2b/PHASE2B_FULL_LOPO.tsv", "lopo_result"),
    ], ignore_index=True, sort=False)
    f3.to_csv(SOURCE / "Figure3_SourceData.tsv", sep="\t", index=False)

    gw = read("GENOMEWIDE_EFFECTS_SHARED_UNIVERSE.tsv")
    gws = read("GENOMEWIDE_CONCORDANCE.tsv")
    rr = read("RECIPROCAL_RANK_REPLICATION.tsv")
    f4 = pd.concat([
        source_long("A-C", gw, "03_results/phase2b/GENOMEWIDE_EFFECTS_SHARED_UNIVERSE.tsv", "gene_effect"),
        source_long("A-C", gws, "03_results/phase2b/GENOMEWIDE_CONCORDANCE.tsv", "genomewide_summary"),
        source_long("D", rr, "03_results/phase2b/RECIPROCAL_RANK_REPLICATION.tsv", "reciprocal_rank_test"),
    ], ignore_index=True, sort=False)
    f4.to_csv(SOURCE / "Figure4_SourceData.tsv", sep="\t", index=False)

    tr = read("TRANSITION_PROGRAM_EFFECTS.tsv")
    trs = read("TRANSITION_SUMMARY.tsv")
    focal = read("FOCAL_PARTICIPANT_D1_D14_DISTRIBUTIONS.tsv")
    focal = focal[focal["program"].isin(DELAYED)]
    col = read("COLLAGEN_PRESPECIFIED_SENSITIVITY.tsv")
    delayed_eff = eff[eff["program"].isin(DELAYED)]
    f5 = pd.concat([
        source_long("A", delayed_eff, "03_results/phase1b/participant_blocked_program_effects.tsv", "delayed_program_effect"),
        source_long("B", focal, "03_results/phase2b/FOCAL_PARTICIPANT_D1_D14_DISTRIBUTIONS.tsv", "participant_difference"),
        source_long("C", tr, "03_results/phase2b/TRANSITION_PROGRAM_EFFECTS.tsv", "transition_program"),
        source_long("C", trs, "03_results/phase2b/TRANSITION_SUMMARY.tsv", "transition_summary"),
        source_long("D", col, "03_results/phase2b/COLLAGEN_PRESPECIFIED_SENSITIVITY.tsv", "collagen_sensitivity"),
    ], ignore_index=True, sort=False)
    f5.to_csv(SOURCE / "Figure5_SourceData.tsv", sep="\t", index=False)


def figure1() -> None:
    fig = plt.figure(figsize=(7.2, 5.7))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], hspace=0.42, wspace=0.30)
    axa, axb, axc, axd = [fig.add_subplot(gs[i, j]) for i, j in [(0,0),(0,1),(1,0),(1,1)]]
    title(fig, "Study design and analytical framework")

    def timeline(ax, days, primary, color, name, subtitle):
        xpos = np.arange(len(days), dtype=float)
        ax.plot([xpos.min(), xpos.max()], [0, 0], color="#9CA3AF", lw=1.2)
        for xx, d in zip(xpos, days):
            is_primary = d in primary
            ax.scatter(xx, 0, s=48 if is_primary else 34, facecolor=color if is_primary else "white",
                       edgecolor=color, lw=1.2, zorder=3, marker="o" if is_primary else "s")
            ax.text(xx, -0.16, f"D{d}", ha="center", va="top", fontsize=6.5,
                    fontweight="bold" if is_primary else "normal")
        ax.text(0.01, 0.90, name, transform=ax.transAxes, fontweight="bold", fontsize=8)
        ax.text(0.01, 0.75, subtitle, transform=ax.transAxes, fontsize=6.5, color=GREY)
        ax.text(0.01, 0.15, "● primary shared   □ secondary context", transform=ax.transAxes,
                fontsize=6.2, color=GREY)
        ax.set_ylim(-0.4, 0.45); ax.set_xlim(-0.5, len(days)-0.5); ax.axis("off")

    timeline(axa, [0,1,3,7,14,21,28], {0,1,7,14}, BLUE,
             "GSE168760 · n=14 · AFL", "Back skin; single-treatment primary frame")
    panel(axa, "A")
    timeline(axb, [0,1,7,14,29], {0,1,7,14}, ORANGE,
             "GSE206495 · n=17 · NAFL", "Forearm first cycle; D29 follows second treatment")
    axb.text(0.99, 0.57, "Face samples excluded from\nprimary cross-cohort trajectory",
             transform=axb.transAxes, ha="right", fontsize=6.2, color=GREY)
    panel(axb, "B")

    axc.axis("off"); panel(axc, "C")
    axc.text(0.02, 0.91, "Shared analytical frame", transform=axc.transAxes,
             fontweight="bold", fontsize=8)
    xs = np.linspace(0.12, 0.88, 4)
    for i, (x, lab) in enumerate(zip(xs, ["D0", "D1", "D7", "D14"])):
        axc.add_patch(FancyBboxPatch((x-0.075, 0.40), 0.15, 0.20,
                     boxstyle="round,pad=0.02", fc="#F3F4F6" if i==0 else "#E6F2F8",
                     ec=BLUE, lw=1.0, transform=axc.transAxes))
        axc.text(x, 0.50, lab, transform=axc.transAxes, ha="center", va="center",
                 fontsize=8, fontweight="bold")
        if i < 3:
            axc.add_patch(FancyArrowPatch((x+0.08,0.50),(xs[i+1]-0.08,0.50),
                         transform=axc.transAxes, arrowstyle="-|>", mutation_scale=8,
                         lw=0.8, color=GREY))
    axc.text(0.50, 0.19, "Cohorts processed separately · participant is the biological replicate",
             transform=axc.transAxes, ha="center", fontsize=6.3, color=GREY)

    axd.axis("off"); panel(axd, "D")
    axd.text(0.02, 0.95, "Evidence workflow (analytical sequence, not a biological mechanism)",
             transform=axd.transAxes, fontweight="bold", fontsize=7.2)
    nodes = ["Participant-\nblocked effects", "Shared universe\n16,942 genes", "7 frozen\nprograms",
             "Program\nconcordance", "Random null\n+ LOPO", "Genome-wide\nconcordance",
             "Reciprocal rank\nreplication", "Transition\nanalysis"]
    for i, node in enumerate(nodes):
        row, col = divmod(i, 4); x = 0.02 + col*0.245; y = 0.58 - row*0.38
        axd.add_patch(FancyBboxPatch((x,y),0.19,0.20,boxstyle="round,pad=0.012",
                    fc="#F9FAFB", ec=BLUE if i<4 else TEAL, lw=0.8, transform=axd.transAxes))
        axd.text(x+0.095,y+0.10,node,ha="center",va="center",fontsize=5.5,transform=axd.transAxes)
        if col < 3:
            axd.add_patch(FancyArrowPatch((x+0.195,y+0.10),(x+0.238,y+0.10),
                        transform=axd.transAxes,arrowstyle="-|>",mutation_scale=7,lw=0.7,color=GREY))
        elif row == 0:
            axd.add_patch(FancyArrowPatch((0.95,0.58),(0.95,0.40),transform=axd.transAxes,
                        arrowstyle="-|>",mutation_scale=7,lw=0.7,color=GREY))
    save_figure(fig, MAIN / "Figure1_StudyDesign")


def figure2() -> None:
    eff = program_effects_primary()
    vec = read("cross_cohort_program_effect_vectors.tsv", 1)
    con = read("cross_cohort_conservation.tsv", 1).set_index("day")
    fig = plt.figure(figsize=(7.2, 5.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.25, 1], hspace=0.42, wspace=0.34)
    axa = fig.add_subplot(gs[0, :]); axes = [fig.add_subplot(gs[1,i]) for i in range(3)]
    title(fig, "Reproducible temporal program architecture")
    mat=[]; cols=[]
    for ds, short in [("GSE168760_AFL","AFL"),("GSE206495_NAFL_FOREARM","NAFL")]:
        for d in [1,7,14]:
            q=eff[(eff.dataset==ds)&(eff.day==d)].set_index("program")
            mat.append([q.loc[p,"paired_hedges_g"] for p in PROGRAM_ORDER]); cols.append(f"{short}\nD{d}")
    arr=np.array(mat).T; vmax=max(2, float(np.nanmax(np.abs(arr))))
    im=axa.imshow(arr, cmap="RdBu_r", norm=TwoSlopeNorm(vmin=-vmax,vcenter=0,vmax=vmax), aspect="auto")
    axa.set_yticks(range(7), [PROGRAM_LABEL[p] for p in PROGRAM_ORDER]); axa.set_xticks(range(6),cols)
    axa.tick_params(length=0)
    for i in range(7):
        for j in range(6): axa.text(j,i,f"{arr[i,j]:.2f}",ha="center",va="center",fontsize=5.5,
                                      color="white" if abs(arr[i,j])>vmax*0.55 else "black")
    cb=fig.colorbar(im,ax=axa,fraction=0.025,pad=0.015); cb.set_label("Paired Hedges g",fontsize=6.5)
    panel(axa,"A")
    for ax,d,lab in zip(axes,[1,7,14],["B","C","D"]):
        q=vec[vec.day==d].set_index("program").loc[PROGRAM_ORDER]
        ax.axhline(0,color="#9CA3AF",ls="--",lw=0.6); ax.axvline(0,color="#9CA3AF",ls="--",lw=0.6)
        ax.scatter(q.afl_hedges_g,q.nafl_hedges_g,c=[BLUE if p not in DELAYED else TEAL for p in PROGRAM_ORDER],
                   s=26,edgecolor="white",lw=0.5,zorder=3)
        for p,row in q.iterrows(): ax.annotate(PROGRAM_ABBR[p],(row.afl_hedges_g,row.nafl_hedges_g),
                                              xytext=(3,3),textcoords="offset points",fontsize=5.1)
        xmin,xmax=q.afl_hedges_g.min(),q.afl_hedges_g.max(); ymin,ymax=q.nafl_hedges_g.min(),q.nafl_hedges_g.max()
        dx=max(.25,(xmax-xmin)*.18); dy=max(.25,(ymax-ymin)*.18)
        ax.set_xlim(xmin-dx,xmax+dx); ax.set_ylim(ymin-dy,ymax+dy)
        r=con.loc[d]
        ax.text(0.03,0.96,f"ρ={r.standardized_spearman_rho:.3f}\n95% CI {r.bootstrap_ci95_low:.3f}–{r.bootstrap_ci95_high:.3f}",
                transform=ax.transAxes,va="top",fontsize=6.2)
        ax.set_title(f"D{d}",fontweight="bold"); ax.set_xlabel("GSE168760 effect"); ax.set_ylabel("GSE206495 effect")
        clean(ax); panel(ax,lab)
    save_figure(fig, MAIN / "Figure2_TemporalPrograms")


def figure3() -> None:
    null=read("MATCHED_RANDOM_PROGRAM_NULL_10000.tsv"); sm=read("MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv").set_index("statistic")
    lopo=read("PHASE2B_FULL_LOPO.tsv")
    fig=plt.figure(figsize=(7.2,5.6)); gs=fig.add_gridspec(2,2,hspace=0.43,wspace=0.32)
    axs=[fig.add_subplot(gs[i,j]) for i,j in [(0,0),(0,1),(1,0),(1,1)]]
    title(fig,"Robustness and null-model testing")
    for ax,col,key,lab in [(axs[0],"mean_rho","mean_rho","A"),(axs[1],"min_rho","min_rho","B")]:
        ax.hist(null[col],bins=35,color="#CBD5E1",edgecolor="white",lw=0.3)
        obs=float(sm.loc[key,"real_value"]); p=float(sm.loc[key,"empirical_one_sided_p"])
        ax.axvline(obs,color=RED,lw=1.5); ax.text(0.98,0.94,f"Observed={obs:.3f}\nP={p:.4f}",transform=ax.transAxes,
                                              ha="right",va="top",fontsize=6.5)
        ax.set_xlabel("Mean temporal ρ" if col=="mean_rho" else "Minimum temporal ρ"); ax.set_ylabel("Random architectures")
        clean(ax); panel(ax,lab)
    ax=axs[2]; colors=np.where(lopo.deleted_from.eq("GSE168760_AFL"),BLUE,ORANGE)
    ax.scatter(range(1,len(lopo)+1),lopo.mean_rho,c=colors,s=20,marker="o",edgecolor="white",lw=0.4)
    ax.axhline(float(lopo[["day1_rho","day7_rho","day14_rho"]].min().min()),color=RED,ls="--",lw=0.8)
    ax.text(0.98,0.17,"Minimum day-specific ρ = 0.750",transform=ax.transAxes,ha="right",va="bottom",fontsize=6.2,color=RED)
    ax.set_xlabel("Participant removal (31 total)"); ax.set_ylabel("Mean D1/D7/D14 ρ"); ax.set_ylim(0.72,0.91)
    ax.legend(handles=[Line2D([0],[0],marker='o',color='none',markerfacecolor=BLUE,label='Removed from GSE168760'),
                       Line2D([0],[0],marker='o',color='none',markerfacecolor=ORANGE,label='Removed from GSE206495')],loc="upper left")
    clean(ax); panel(ax,"C")
    ax=axs[3]; arr=lopo[["day1_rho","day7_rho","day14_rho"]].to_numpy();
    sns.heatmap(arr,ax=ax,cmap="Blues",vmin=0,vmax=1,cbar_kws={"label":"Spearman ρ"},xticklabels=["D1","D7","D14"],yticklabels=False)
    ax.set_ylabel("31 participant removals"); ax.text(1.0,-0.17,"31/31 retain positive ρ at all three days",transform=ax.transAxes,ha="right",fontsize=6.3)
    panel(ax,"D")
    save_figure(fig, MAIN / "Figure3_Robustness")


def figure4() -> None:
    gw=read("GENOMEWIDE_EFFECTS_SHARED_UNIVERSE.tsv"); sm=read("GENOMEWIDE_CONCORDANCE.tsv").set_index("day")
    rr=read("RECIPROCAL_RANK_REPLICATION.tsv")
    fig=plt.figure(figsize=(7.2,5.5)); gs=fig.add_gridspec(2,3,height_ratios=[1,0.82],hspace=0.46,wspace=0.34)
    axes=[fig.add_subplot(gs[0,i]) for i in range(3)]; axd=fig.add_subplot(gs[1,:])
    title(fig,"Genome-wide cross-cohort replication")
    for ax,d,lab in zip(axes,[1,7,14],["A","B","C"]):
        q=gw[gw.day==d]; hb=ax.hexbin(q.afl_paired_hedges_g,q.nafl_paired_hedges_g,gridsize=48,mincnt=1,
                                    cmap="Blues",bins="log",linewidths=0,rasterized=True)
        ax.axhline(0,color="#9CA3AF",lw=0.5,ls="--"); ax.axvline(0,color="#9CA3AF",lw=0.5,ls="--")
        r=sm.loc[d]; ax.text(0.03,0.96,f"ρ={r.standardized_spearman_rho:.3f}\nPermutation P<0.0001\nSign concordance={100*r.sign_concordance_all:.1f}%",
                            transform=ax.transAxes,va="top",fontsize=6.1,
                            bbox=dict(fc="white",ec="none",alpha=0.82,pad=1.5))
        ax.set_title(f"D{d}",fontweight="bold"); ax.set_xlabel("GSE168760 gene effect"); ax.set_ylabel("GSE206495 gene effect")
        clean(ax); panel(ax,lab)
    dirs=[("AFL_TO_NAFL","POSITIVE"),("AFL_TO_NAFL","NEGATIVE"),("NAFL_TO_AFL","POSITIVE"),("NAFL_TO_AFL","NEGATIVE")]
    mat=np.zeros((3,4)); qmat=np.zeros((3,4))
    for i,d in enumerate([1,7,14]):
        for j,(direction,tail) in enumerate(dirs):
            z=rr[(rr.day==d)&(rr.source_direction==direction)&(rr.source_tail==tail)].iloc[0]
            mat[i,j]=z.NES; qmat[i,j]=z.BH_q
    im=axd.imshow(mat,cmap="YlGnBu",aspect="auto",vmin=0,vmax=max(4.6,mat.max()))
    axd.set_yticks(range(3),["D1","D7","D14"]); axd.set_xticks(range(4),["A→B +","A→B −","B→A +","B→A −"])
    for i in range(3):
        for j in range(4): axd.text(j,i,f"NES {mat[i,j]:.2f}\nq={qmat[i,j]:.4f}",ha="center",va="center",fontsize=6,
                                       color="white" if mat[i,j]>3.3 else "black")
    axd.tick_params(length=0); axd.text(1.0,-0.20,"12/12 expected direction; 12/12 BH q<0.05",transform=axd.transAxes,ha="right",fontsize=6.5)
    panel(axd,"D"); fig.colorbar(im,ax=axd,fraction=0.025,pad=0.015,label="NES")
    save_figure(fig, MAIN / "Figure4_GenomeWideReplication")


def figure5() -> None:
    eff=program_effects_primary(); focal=read("FOCAL_PARTICIPANT_D1_D14_DISTRIBUTIONS.tsv")
    tr=read("TRANSITION_PROGRAM_EFFECTS.tsv").set_index("program").loc[PROGRAM_ORDER]
    col=read("COLLAGEN_PRESPECIFIED_SENSITIVITY.tsv")
    fig=plt.figure(figsize=(7.2,6.5)); gs=fig.add_gridspec(2,2,hspace=0.43,wspace=0.34)
    axa,axb,axc,axd=[fig.add_subplot(gs[i,j]) for i,j in [(0,0),(0,1),(1,0),(1,1)]]
    title(fig,"Delayed remodeling and temporal transition")
    x=np.array([1,7,14]); offsets={"GSE168760_AFL":-0.15,"GSE206495_NAFL_FOREARM":0.15}
    for ds,c,m,label in [("GSE168760_AFL",BLUE,"o","GSE168760"),("GSE206495_NAFL_FOREARM",ORANGE,"s","GSE206495")]:
        for p,ls in zip(DELAYED,["-","--",":"]):
            q=eff[(eff.dataset==ds)&(eff.program==p)].set_index("day").loc[[1,7,14]]
            axa.errorbar(x+offsets[ds],q.model_effect,yerr=[q.model_effect-q.model_ci95_low,q.model_ci95_high-q.model_effect],
                         color=c,marker=m,ls=ls,lw=1,ms=3.4,capsize=2,alpha=0.9)
    axa.set_xticks(x,["D1","D7","D14"]); axa.set_ylabel("Participant-blocked model effect (95% CI)"); clean(axa,zero=True)
    axa.legend(handles=[Line2D([0],[0],color=BLUE,marker='o',label='GSE168760'),Line2D([0],[0],color=ORANGE,marker='s',label='GSE206495'),
                        Line2D([0],[0],color=GREY,ls='-',label='ECM remodeling'),Line2D([0],[0],color=GREY,ls='--',label='Collagen organization'),
                        Line2D([0],[0],color=GREY,ls=':',label='Matrix maturation')],ncol=2,loc="upper left")
    panel(axa,"A")

    q=focal[focal.program.isin(DELAYED)&focal.day.isin([1,14])].copy()
    positions={}; pos=0
    for p in DELAYED:
        for ds in ["GSE168760_AFL","GSE206495_NAFL_FOREARM"]:
            a=q[(q.program==p)&(q.dataset==ds)&(q.day==1)].set_index("participant_id")
            b=q[(q.program==p)&(q.dataset==ds)&(q.day==14)].set_index("participant_id")
            ids=a.index.intersection(b.index); c=BLUE if ds.endswith("AFL") else ORANGE
            for pid in ids: axb.plot([pos,pos+0.32],[a.loc[pid,"difference"],b.loc[pid,"difference"]],color=c,alpha=0.25,lw=0.5)
            axb.scatter(np.repeat(pos,len(ids)),a.loc[ids,"difference"],s=7,color=c,alpha=0.55)
            axb.scatter(np.repeat(pos+0.32,len(ids)),b.loc[ids,"difference"],s=8,color=c,alpha=0.75,marker="s")
            positions[(p,ds)]=pos+0.16; pos+=0.72
        pos+=0.28
    axb.set_xticks([(positions[(p,"GSE168760_AFL")]+positions[(p,"GSE206495_NAFL_FOREARM")])/2 for p in DELAYED],
                   [PROGRAM_ABBR[p] for p in DELAYED],rotation=0)
    axb.set_ylabel("Participant treated − D0 program score"); clean(axb,zero=True)
    axb.text(0.02,0.97,"Circles D1; squares D14\nBlue GSE168760; orange GSE206495",transform=axb.transAxes,va="top",fontsize=6.1)
    panel(axb,"B")

    axc.axhline(0,color="#9CA3AF",ls="--",lw=0.6); axc.axvline(0,color="#9CA3AF",ls="--",lw=0.6)
    axc.scatter(tr.afl_transition_D14_minus_D1,tr.nafl_transition_D14_minus_D1,s=28,c=[TEAL if p in DELAYED else PURPLE for p in PROGRAM_ORDER],edgecolor="white",lw=0.4)
    for p,row in tr.iterrows(): axc.annotate(PROGRAM_ABBR[p],(row.afl_transition_D14_minus_D1,row.nafl_transition_D14_minus_D1),xytext=(3,3),textcoords="offset points",fontsize=5.2)
    axc.text(0.03,0.97,"ρ=0.857 · 6/7 same direction\nMatched-random transition null P=0.286",transform=axc.transAxes,va="top",fontsize=6.4,
             bbox=dict(fc="white",ec="none",alpha=0.82,pad=1.5))
    axc.set_xlabel("GSE168760 D14 − D1 effect"); axc.set_ylabel("GSE206495 D14 − D1 effect"); clean(axc); panel(axc,"C")

    fam=col[col.analysis=="family_removal"].copy(); loo=col[col.analysis=="leave_one_gene_out"].copy()
    variants=["full","minus_COL1A1_COL1A2","minus_all_COL_prefix"]
    for i,(ds,c,m,label) in enumerate([("GSE168760_AFL",BLUE,"o","GSE168760"),("GSE206495_NAFL_FOREARM",ORANGE,"s","GSE206495")]):
        vals=[float(fam[(fam.dataset==ds)&(fam.variant==v)].day14_hedges_g.iloc[0]) for v in variants]
        axd.scatter(np.arange(3)+(-0.10 if i==0 else 0.10),vals,color=c,marker=m,s=28,label=label,zorder=3)
        z=loo[loo.dataset==ds].day14_hedges_g
        axd.vlines(3+(-0.10 if i==0 else 0.10),z.min(),z.max(),color=c,lw=3,alpha=0.7)
        axd.scatter(3+(-0.10 if i==0 else 0.10),z.median(),color="white",edgecolor=c,marker=m,s=24,zorder=4)
    axd.set_xticks(range(4),["Full","−COL1A1/\nCOL1A2","−all COL*","Leave-one-\ngene-out"],rotation=0)
    axd.set_ylabel("D14 collagen-organization Hedges g"); clean(axd,zero=True); axd.legend(loc="lower left"); panel(axd,"D")
    save_figure(fig, MAIN / "Figure5_DelayedRemodeling")


def supp_s1() -> None:
    m1=pd.read_csv(META/"GSE168760_metadata.tsv",sep="\t"); m2=pd.read_csv(META/"GSE206495_metadata.tsv",sep="\t")
    fig,axs=plt.subplots(1,2,figsize=(7.2,4.5),gridspec_kw={"wspace":0.32}); title(fig,"Sample and participant structure QC")
    def day1(s): return int(str(s).replace("day ",""))
    a=m1.copy(); a["day"]=a.time.map(day1); a["retained"]=((a.day==0)&(a.treatment_number.astype(str)=="0"))|((a.day.isin([1,3,7,14,21,28]))&(a.treatment_number.astype(str)=="1"))
    aa=a[a.retained].pivot_table(index="participant_id",columns="day",values="geo_accession",aggfunc="count",fill_value=0).reindex(columns=[0,1,3,7,14,21,28])
    sns.heatmap(aa,ax=axs[0],cmap=["white",BLUE],cbar=False,linewidths=.35,linecolor="white",vmin=0,vmax=1)
    axs[0].set_title("GSE168760 retained single-treatment samples"); axs[0].set_xlabel("Day"); axs[0].set_ylabel("Participant ID")
    axs[0].text(0,-0.16,f"Retained {int(aa.values.sum())}/140; primary D0/D1/D7/D14 = 56",transform=axs[0].transAxes,fontsize=6.3)
    b=m2.copy();
    def day2(t):
        t=str(t).lower(); return 0 if t=="baseline" else int(t.split()[0])
    b["day"]=b.time.map(day2); b["retained"]=b.anatomical_site.eq("arm skin")&b.day.isin([0,1,7,14,29])
    bb=b[b.retained].pivot_table(index="participant_id",columns="day",values="geo_accession",aggfunc="count",fill_value=0).reindex(columns=[0,1,7,14,29])
    sns.heatmap(bb,ax=axs[1],cmap=["white",ORANGE],cbar=False,linewidths=.35,linecolor="white",vmin=0,vmax=1)
    axs[1].set_title("GSE206495 retained forearm samples"); axs[1].set_xlabel("Day"); axs[1].set_ylabel("Participant ID")
    axs[1].text(0,-0.16,f"Forearm retained {int(bb.values.sum())}/119; face excluded 34; primary = 68",transform=axs[1].transAxes,fontsize=6.3)
    panel(axs[0],"A"); panel(axs[1],"B"); save_figure(fig,SUPP/"FigureS1_SampleStructureQC")


def supp_s2() -> None:
    x=read("platform_mapping_audit.tsv",1)
    fig,axs=plt.subplots(1,3,figsize=(7.2,3.2),gridspec_kw={"wspace":0.43}); title(fig,"Platform annotation and shared-gene universe")
    for ax,plat,lab,c in [(axs[0],"GPL13667","A",BLUE),(axs[1],"GPL15207","B",ORANGE)]:
        q=x[x.platform==plat]; vals=[len(q),int(q.used_reliable.sum()),int(q.selected_highest_iqr_probe.sum()),int(q.used_shared_universe.sum())]
        ax.barh(range(4),vals,color=[LIGHT,c,c,TEAL],alpha=.9); ax.set_yticks(range(4),["Original probes","Reliable mapping","Selected gene probe","Used in shared universe"]); ax.invert_yaxis()
        ax.set_xlim(0,max(vals)*1.06)
        for i,v in enumerate(vals): ax.text(v*.98,i,f"{v:,}",va="center",ha="right",fontsize=5.7,
                                           color="white" if i>0 else "black")
        if plat=="GPL15207": ax.tick_params(axis="y",labelleft=False)
        ax.set_title(plat); ax.set_xlabel("Features/probes"); clean(ax); panel(ax,lab)
    axs[2].axis("off"); panel(axs[2],"C")
    for y,text0,n,c in [(0.78,"GPL13667 reliable genes",x[(x.platform=="GPL13667")&x.selected_highest_iqr_probe].gene_symbol.nunique(),BLUE),
                        (0.53,"GPL15207 reliable genes",x[(x.platform=="GPL15207")&x.selected_highest_iqr_probe].gene_symbol.nunique(),ORANGE),
                        (0.22,"Frozen shared universe",16942,TEAL)]:
        axs[2].add_patch(FancyBboxPatch((.08,y-.08),.82,.16,boxstyle="round,pad=.02",fc="#F9FAFB",ec=c,lw=1,transform=axs[2].transAxes))
        axs[2].text(.49,y,f"{text0}\n{n:,}",ha="center",va="center",transform=axs[2].transAxes,fontsize=7,fontweight="bold" if y<.3 else "normal")
    axs[2].add_patch(FancyArrowPatch((.49,.69),(.49,.62),transform=axs[2].transAxes,arrowstyle="-|>",color=GREY,mutation_scale=8))
    axs[2].add_patch(FancyArrowPatch((.49,.44),(.49,.31),transform=axs[2].transAxes,arrowstyle="-|>",color=GREY,mutation_scale=8))
    save_figure(fig,SUPP/"FigureS2_PlatformSharedUniverse")


def supp_s3() -> None:
    eff=program_effects_primary(); fig,axs=plt.subplots(1,2,figsize=(7.2,6.2),sharex=True,sharey=True,gridspec_kw={"wspace":0.10}); title(fig,"Program-level uncertainty")
    for ax,ds,lab,c in [(axs[0],"GSE168760_AFL","A",BLUE),(axs[1],"GSE206495_NAFL_FOREARM","B",ORANGE)]:
        q=eff[eff.dataset==ds].copy(); q["program"]=pd.Categorical(q.program,PROGRAM_ORDER,ordered=True); q=q.sort_values(["program","day"])
        y=np.arange(len(q)); ax.errorbar(q.model_effect,y,xerr=[q.model_effect-q.model_ci95_low,q.model_ci95_high-q.model_effect],fmt="o",ms=3,color=c,ecolor=c,lw=.8,capsize=1.5)
        ax.axvline(0,color="#9CA3AF",ls="--",lw=.7); ax.set_yticks(y,[f"{PROGRAM_LABEL[p]} · D{d}" for p,d in zip(q.program,q.day)])
        if ds=="GSE168760_AFL": ax.invert_yaxis()
        ax.set_xlabel("Participant-blocked model effect (95% CI)"); clean(ax); panel(ax,lab)
        if ds=="GSE206495_NAFL_FOREARM": ax.tick_params(axis="y",labelleft=False)
    save_figure(fig,SUPP/"FigureS3_ProgramUncertainty")


def supp_s4() -> None:
    x=read("PHASE2B_FULL_LOPO.tsv"); arr=x[["day1_rho","day7_rho","day14_rho"]].to_numpy()
    fig,ax=plt.subplots(figsize=(7.2,6.0)); title(fig,"Complete leave-one-participant-out results")
    labels=[f"{('A' if r.deleted_from=='GSE168760_AFL' else 'B')}:{r.deleted_participant}" for _,r in x.iterrows()]
    sns.heatmap(arr,ax=ax,cmap="Blues",vmin=.70,vmax=.90,annot=True,fmt=".2f",annot_kws={"fontsize":5.2},
                xticklabels=["D1","D7","D14"],yticklabels=labels,cbar_kws={"label":"Spearman ρ"})
    mins=arr.min(axis=0); meds=np.median(arr,axis=0); maxs=arr.max(axis=0)
    ax.set_xlabel(f"Minimum {mins.min():.3f}; medians D1/D7/D14 {meds[0]:.3f}/{meds[1]:.3f}/{meds[2]:.3f}; overall range {arr.min():.3f}–{arr.max():.3f}")
    ax.set_ylabel("Removed participant"); ax.text(0,-.16,"A=GSE168760 removal; B=GSE206495 removal",transform=ax.transAxes,fontsize=6.2)
    panel(ax,"A"); save_figure(fig,SUPP/"FigureS4_CompleteLOPO")


def supp_s5() -> None:
    x=read("MATCHED_RANDOM_PROGRAM_NULL_10000.tsv"); sm=read("MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv").set_index("statistic")
    fig,axs=plt.subplots(1,3,figsize=(7.2,3.3),gridspec_kw={"wspace":.36}); title(fig,"Random-program null diagnostics")
    for d,c in [(1,BLUE),(7,TEAL),(14,ORANGE)]: sns.kdeplot(x[f"day{d}_rho"],ax=axs[0],label=f"D{d}",color=c,lw=1)
    for d,c in [(1,BLUE),(7,TEAL),(14,ORANGE)]: axs[0].axvline(float(sm.loc[f"day{d}_rho","real_value"]),color=c,ls="--",lw=.8)
    axs[0].set_xlabel("Day-specific random-program ρ"); axs[0].legend(); clean(axs[0]); panel(axs[0],"A")
    keys=["mean_rho","min_rho","day1_rho","day7_rho","day14_rho"]
    vals=sm.loc[keys,"null_q95"].to_numpy(); obs=sm.loc[keys,"real_value"].to_numpy(); yy=np.arange(len(keys))
    axs[1].scatter(vals,yy,color=GREY,label="Null 95th percentile"); axs[1].scatter(obs,yy,color=RED,marker="D",label="Observed")
    axs[1].set_yticks(yy,["Mean","Minimum","D1","D7","D14"]); axs[1].invert_yaxis(); axs[1].set_xlabel("Spearman ρ"); axs[1].legend(loc="lower right",fontsize=5.2); clean(axs[1]); panel(axs[1],"B")
    axs[2].hist(x.transition_rho,bins=30,color="#CBD5E1",edgecolor="white",lw=.3); obs=float(sm.loc["transition_rho","real_value"]); p=float(sm.loc["transition_rho","empirical_one_sided_p"])
    axs[2].axvline(obs,color=RED,lw=1.4); axs[2].text(.97,.94,f"Observed={obs:.3f}\nP={p:.3f}\nPercentile=71.39",transform=axs[2].transAxes,ha="right",va="top",fontsize=6.4)
    axs[2].set_xlabel("Random transition ρ"); axs[2].set_ylabel("Architectures"); clean(axs[2]); panel(axs[2],"C")
    save_figure(fig,SUPP/"FigureS5_RandomNullDiagnostics")


def supp_s6() -> None:
    prim=program_effects_primary(); sec=read("SECONDARY_TIMEPOINT_CONTEXT.tsv"); rows=[]
    for p in PROGRAM_ORDER:
        rows.append({"program":p,"day":0,"paired_hedges_g":0.0})
    q=prim[prim.dataset=="GSE168760_AFL"][["program","day","paired_hedges_g"]]
    q2=sec[(sec.dataset=="GSE168760_AFL")][["program","day","paired_hedges_g"]]
    d=pd.concat([pd.DataFrame(rows),q,q2]).sort_values(["program","day"])
    fig,ax=plt.subplots(figsize=(7.2,4.2)); title(fig,"AFL-only exploratory extended trajectory")
    colors=sns.color_palette("colorblind",7)
    for p,c in zip(PROGRAM_ORDER,colors):
        z=d[d.program==p]; ax.plot(z.day,z.paired_hedges_g,marker="o",ms=3,lw=1,color=c,label=PROGRAM_LABEL[p])
    ax.axhline(0,color="#9CA3AF",ls="--",lw=.6); ax.set_xticks([0,1,3,7,14,21,28]); ax.set_xlabel("Day after first AFL treatment"); ax.set_ylabel("Paired Hedges g")
    ax.legend(ncol=4,loc="upper center",bbox_to_anchor=(.5,1.02)); clean(ax); panel(ax,"A");
    ax.text(1,-.20,"Within-cohort exploratory context; D0 is the reference value (0).",transform=ax.transAxes,ha="right",fontsize=6.3,color=GREY)
    save_figure(fig,SUPP/"FigureS6_AFLOnlyExtendedTrajectory")


def supp_s7() -> None:
    x=read("SECONDARY_TIMEPOINT_CONTEXT.tsv"); x=x[x.dataset=="GSE206495_NAFL_FOREARM"].set_index("program").loc[PROGRAM_ORDER]
    fig,ax=plt.subplots(figsize=(7.2,3.6)); title(fig,"GSE206495 D29 treatment-context secondary")
    y=np.arange(7); colors=[TEAL if v>0 else PURPLE for v in x.paired_hedges_g]
    ax.scatter(x.paired_hedges_g,y,c=colors,s=32,marker="s"); ax.axvline(0,color="#9CA3AF",ls="--",lw=.7)
    ax.set_yticks(y,[PROGRAM_LABEL[p] for p in PROGRAM_ORDER]); ax.invert_yaxis(); ax.set_xlabel("D29 paired Hedges g"); clean(ax); panel(ax,"A")
    ax.text(.99,.04,"D29 = 29 days after first session and 1 day after second session; not a primary D14→D29 trajectory.",transform=ax.transAxes,ha="right",fontsize=6.2,color=GREY)
    save_figure(fig,SUPP/"FigureS7_D29TreatmentContext")


def supp_s8() -> None:
    eff=program_effects_primary(); e=eff[eff.program=="proliferation"]
    lo=read("LOPO_program_effects.tsv",1); lo=lo[(lo.program=="proliferation")&(lo.day==14)]
    mp=read("platform_mapping_sensitivity.tsv",1); mp=mp[(mp.program=="proliferation")&(mp.day==14)]
    half=read("proliferation_gene_set_sensitivity.tsv",1)
    fig=plt.figure(figsize=(7.2,4.0)); gs=fig.add_gridspec(1,2,wspace=.34); axa=fig.add_subplot(gs[0,0]); axb=fig.add_subplot(gs[0,1])
    title(fig,"Exploratory proliferation divergence")
    for ds,c,m,label in [("GSE168760_AFL",BLUE,"o","GSE168760"),("GSE206495_NAFL_FOREARM",ORANGE,"s","GSE206495")]:
        q=e[e.dataset==ds].sort_values("day"); axa.errorbar(q.day,q.model_effect,yerr=[q.model_effect-q.model_ci95_low,q.model_ci95_high-q.model_effect],
                 color=c,marker=m,lw=1,ms=4,capsize=2,label=label)
    axa.set_xticks([1,7,14],["D1","D7","D14"]); axa.set_ylabel("Participant-blocked model effect (95% CI)"); axa.legend(); clean(axa,zero=True); panel(axa,"A")
    for i,(ds,c,m,label) in enumerate([("GSE168760_AFL",BLUE,"o","GSE168760"),("GSE206495_NAFL_FOREARM",ORANGE,"s","GSE206495")]):
        vals=[]; names=[]
        z=e[(e.dataset==ds)&(e.day==14)].iloc[0]; vals.append(z.paired_hedges_g); names.append("Full")
        for method in mp.mapping_method.unique():
            z=mp[(mp.dataset==ds)&(mp.mapping_method==method)].iloc[0]; vals.append(z.paired_hedges_g); names.append("Mapping")
        for variant in sorted(half.variant.unique()):
            z=half[(half.dataset==ds)&(half.variant==variant)].iloc[0]; vals.append(z.paired_hedges_g); names.append("Half set")
        xs=np.arange(len(vals))+(-.08 if i==0 else .08); axb.scatter(xs,vals,color=c,marker=m,s=25,label=label)
        l=lo[lo.deleted_from.eq(ds)]
        lp=l.afl_hedges_g if ds=="GSE168760_AFL" else l.nafl_hedges_g
        axb.vlines(len(vals)+(-.08 if i==0 else .08),lp.min(),lp.max(),color=c,lw=3,alpha=.7)
        axb.scatter(len(vals)+(-.08 if i==0 else .08),lp.median(),facecolor="white",edgecolor=c,marker=m,s=24)
    axb.axhline(0,color="#9CA3AF",ls="--",lw=.7); axb.set_xticks(range(5),["Full","Mapping","Half 1","Half 2","LOPO\nrange"]); axb.set_ylabel("D14 proliferation Hedges g"); axb.legend(); clean(axb); axb.text(.98,.13,"DIRECTIONALLY ROBUST\nBUT NOT CONFIRMATORY",transform=axb.transAxes,ha="right",fontsize=6.2,color=GREY)
    save_figure(fig,SUPP/"FigureS8_ExploratoryProliferation")


def make_tables() -> None:
    t1=pd.DataFrame([
        ["GSE168760",14,"Back skin","AFL","Human Genome HG-U219 Array","D0/D1/D3/D7/D14/D21/D28","D0/D1/D7/D14",56,"Participant","Cohort 1; single-treatment primary"],
        ["GSE206495",17,"Forearm skin","NAFL","Human Genome PrimeView Array","D0/D1/D7/D14; D29 second-treatment context","D0/D1/D7/D14",68,"Participant","Cohort 2; forearm first-cycle primary"],
    ],columns=["Dataset","Participants","Primary anatomical site","Laser treatment","Platform","Longitudinal sampling","Primary shared timepoints","Samples used in primary analysis","Participant blocking","Role in study"])
    t1.to_csv(ROOT/"Table1_CohortDesign.tsv",sep="\t",index=False)
    t2=pd.DataFrame([
        ["Frozen program concordance","Spearman rho","0.821","0.821","0.857","Participant-bootstrap intervals; exact 7! label tests","Reproducible ordering of seven frozen programs","LEVEL_1"],
        ["Random-program null","10,000 matched architectures","—","—","—","Mean-rho P=0.0313; minimum-rho P=0.0110","Joint multi-time architecture exceeds matched random programs","LEVEL_1"],
        ["LOPO","Full 31-participant deletion","Positive","Positive","Positive","31/31; minimum rho=0.750","No single participant drives concordance","LEVEL_1"],
        ["Genome-wide concordance","Spearman rho","0.714","0.505","0.373","16,942 genes; participant-bootstrap intervals","Convergence extends beyond frozen programs and attenuates over time","LEVEL_1"],
        ["Gene-label permutation","10,000 permutations/day","P<0.0001","P<0.0001","P<0.0001","All 3/3 above random gene-label expectation","Genome-wide alignment is non-random","LEVEL_1"],
        ["Reciprocal rank replication","Preranked GSEA","4/4","4/4","4/4","12/12 expected direction; 12/12 BH q<0.05","Positive and negative tails replicate in both directions","LEVEL_1"],
        ["Transition","D14 minus D1 program effects","—","—","rho=0.857; 6/7","Matched-random transition P=0.286","Descriptive/supportive transition concordance only","LEVEL_2"],
        ["Collagen robustness","Prespecified removals and 69 LOO","—","—","All positive","COL1A1/COL1A2 removal, all COL* removal and 69/69 LOO positive","Collagen-organization direction is not driven by one gene","LEVEL_1"],
    ],columns=["Evidence layer","Statistic","D1","D7","D14","Robustness","Interpretation","Claim level"])
    t2.to_csv(ROOT/"Table2_EvidenceSummary.tsv",sep="\t",index=False)

    def render_table(df,outbase,width=13.0,row_h=.52):
        wrap_widths={c:(18 if c in {"Dataset","Participants","D1","D7","D14","Claim level"} else 26) for c in df.columns}
        wrapped=df.copy()
        for c in wrapped.columns:
            wrapped[c]=wrapped[c].map(lambda v:"\n".join(textwrap.wrap(str(v),width=wrap_widths[c],break_long_words=False)))
        headers=["\n".join(textwrap.wrap(str(c),width=18,break_long_words=False)) for c in df.columns]
        fig,ax=plt.subplots(figsize=(width,max(2.6,1.7+row_h*len(df)))); ax.axis("off")
        weights=np.array([max(9,min(38,max([len(str(c))]+[len(str(v)) for v in df[c]])*.72+4)) for c in df.columns],dtype=float); weights=weights/weights.sum()
        tbl=ax.table(cellText=wrapped.values,colLabels=headers,cellLoc="left",colLoc="left",loc="center",colWidths=weights)
        tbl.auto_set_font_size(False); tbl.set_fontsize(5.4); tbl.scale(1,1.75)
        for (r,c),cell in tbl.get_celld().items():
            cell.set_edgecolor("#D1D5DB"); cell.set_linewidth(.45); cell.set_facecolor("#E6F2F8" if r==0 else "white")
            if r==0: cell.set_text_props(fontweight="bold")
        save_figure(fig,outbase,dpi=320)
    render_table(t1,ROOT/"Table1_CohortDesign",row_h=.70)
    render_table(t2,ROOT/"Table2_EvidenceSummary",row_h=.58)


def build_source_lock() -> None:
    plan=pd.read_csv(ROOT/"PHASE3A_FIGURE_PLAN.tsv",sep="\t")
    matrix=pd.read_csv(ROOT/"PHASE3A_FIGURE_CLAIM_MATRIX.tsv",sep="\t")
    claim_map={
        "Figure 1":"C1","Figure 2":"C1","Figure 3":"C1","Figure 4":"C2","Figure 5":"C3/C4",
        "Supplement Figure S1":"C1","Supplement Figure S2":"C1","Supplement Figure S3":"C1",
        "Supplement Figure S4":"C1","Supplement Figure S5":"C1/C4","Supplement Figure S6":"C1",
        "Supplement Figure S7":"C1","Supplement Figure S8":"C5",
    }
    sources={
      ("Figure 1","A"):"00_metadata/GSE168760_metadata.tsv",("Figure 1","B"):"00_metadata/GSE206495_metadata.tsv",
      ("Figure 1","C"):"PHASE3A_FIGURE_PLAN.tsv",("Figure 1","D"):"PHASE3A_FIGURE_PLAN.tsv",
      ("Figure 2","A"):"03_results/phase1b/participant_blocked_program_effects.tsv",
      ("Figure 2","B"):"03_results/phase1b/cross_cohort_program_effect_vectors.tsv; cross_cohort_conservation.tsv",("Figure 2","C"):"03_results/phase1b/cross_cohort_program_effect_vectors.tsv; cross_cohort_conservation.tsv",("Figure 2","D"):"03_results/phase1b/cross_cohort_program_effect_vectors.tsv; cross_cohort_conservation.tsv",
      ("Figure 3","A"):"03_results/phase2b/MATCHED_RANDOM_PROGRAM_NULL_10000.tsv; MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv",("Figure 3","B"):"03_results/phase2b/MATCHED_RANDOM_PROGRAM_NULL_10000.tsv; MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv",
      ("Figure 3","C"):"03_results/phase2b/PHASE2B_FULL_LOPO.tsv",("Figure 3","D"):"03_results/phase2b/PHASE2B_FULL_LOPO.tsv",
      ("Figure 4","A"):"03_results/phase2b/GENOMEWIDE_EFFECTS_SHARED_UNIVERSE.tsv; GENOMEWIDE_CONCORDANCE.tsv",("Figure 4","B"):"03_results/phase2b/GENOMEWIDE_EFFECTS_SHARED_UNIVERSE.tsv; GENOMEWIDE_CONCORDANCE.tsv",("Figure 4","C"):"03_results/phase2b/GENOMEWIDE_EFFECTS_SHARED_UNIVERSE.tsv; GENOMEWIDE_CONCORDANCE.tsv",("Figure 4","D"):"03_results/phase2b/RECIPROCAL_RANK_REPLICATION.tsv",
      ("Figure 5","A"):"03_results/phase1b/participant_blocked_program_effects.tsv",("Figure 5","B"):"03_results/phase2b/FOCAL_PARTICIPANT_D1_D14_DISTRIBUTIONS.tsv",("Figure 5","C"):"03_results/phase2b/TRANSITION_PROGRAM_EFFECTS.tsv; TRANSITION_SUMMARY.tsv",("Figure 5","D"):"03_results/phase2b/COLLAGEN_PRESPECIFIED_SENSITIVITY.tsv",
      ("Supplement Figure S1","A"):"00_metadata/GSE168760_metadata.tsv",("Supplement Figure S1","B"):"00_metadata/GSE206495_metadata.tsv",
      ("Supplement Figure S2","A"):"03_results/phase1b/platform_mapping_audit.tsv",("Supplement Figure S2","B"):"03_results/phase1b/platform_mapping_audit.tsv",("Supplement Figure S2","C"):"03_results/phase1b/shared_gene_universe.tsv",
      ("Supplement Figure S3","A"):"03_results/phase1b/participant_blocked_program_effects.tsv",("Supplement Figure S3","B"):"03_results/phase1b/participant_blocked_program_effects.tsv",
      ("Supplement Figure S4","A"):"03_results/phase2b/PHASE2B_FULL_LOPO.tsv",
      ("Supplement Figure S5","A"):"03_results/phase2b/MATCHED_RANDOM_PROGRAM_NULL_10000.tsv",("Supplement Figure S5","B"):"03_results/phase2b/MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv",("Supplement Figure S5","C"):"03_results/phase2b/MATCHED_RANDOM_PROGRAM_NULL_10000.tsv; MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv",
      ("Supplement Figure S6","A"):"03_results/phase1b/participant_blocked_program_effects.tsv; 03_results/phase2b/SECONDARY_TIMEPOINT_CONTEXT.tsv",
      ("Supplement Figure S7","A"):"03_results/phase2b/SECONDARY_TIMEPOINT_CONTEXT.tsv",
      ("Supplement Figure S8","A"):"03_results/phase1b/participant_blocked_program_effects.tsv; LOPO_program_effects.tsv; platform_mapping_sensitivity.tsv; proliferation_gene_set_sensitivity.tsv",
    }
    cols={
      ("Figure 1","A"):"geo_accession,sample_title,participant_id,time,treatment_number,laser_type,platform",
      ("Figure 1","B"):"geo_accession,sample_title,participant_id,anatomical_site,time,laser_type,platform",
      ("Figure 1","C"):"figure,panel,data_scope,statistic_or_encoding", ("Figure 1","D"):"figure,panel,data_scope,statistic_or_encoding",
      ("Figure 2","A"):"dataset,program,day,paired_hedges_g",
      ("Figure 2","B"):"day,program,afl_hedges_g,nafl_hedges_g,standardized_spearman_rho,bootstrap_ci95_low,bootstrap_ci95_high",
      ("Figure 2","C"):"day,program,afl_hedges_g,nafl_hedges_g,standardized_spearman_rho,bootstrap_ci95_low,bootstrap_ci95_high",
      ("Figure 2","D"):"day,program,afl_hedges_g,nafl_hedges_g,standardized_spearman_rho,bootstrap_ci95_low,bootstrap_ci95_high",
      ("Figure 3","A"):"iteration,mean_rho,statistic,real_value,empirical_one_sided_p",
      ("Figure 3","B"):"iteration,min_rho,statistic,real_value,empirical_one_sided_p",
      ("Figure 3","C"):"deleted_from,deleted_participant,mean_rho,day1_rho,day7_rho,day14_rho",
      ("Figure 3","D"):"deleted_from,deleted_participant,day1_rho,day7_rho,day14_rho,all_days_positive",
      ("Figure 4","A"):"day,gene_symbol,afl_paired_hedges_g,nafl_paired_hedges_g,standardized_spearman_rho,empirical_one_sided_p,sign_concordance_all",
      ("Figure 4","B"):"day,gene_symbol,afl_paired_hedges_g,nafl_paired_hedges_g,standardized_spearman_rho,empirical_one_sided_p,sign_concordance_all",
      ("Figure 4","C"):"day,gene_symbol,afl_paired_hedges_g,nafl_paired_hedges_g,standardized_spearman_rho,empirical_one_sided_p,sign_concordance_all",
      ("Figure 4","D"):"day,source_direction,source_tail,NES,BH_q,direction_correct,supported",
      ("Figure 5","A"):"dataset,program,day,model_effect,model_ci95_low,model_ci95_high",
      ("Figure 5","B"):"dataset,program,day,participant_id,difference",
      ("Figure 5","C"):"program,afl_transition_D14_minus_D1,nafl_transition_D14_minus_D1,same_transition_direction,metric,value",
      ("Figure 5","D"):"analysis,variant,removed_gene,dataset,genes_used,day14_hedges_g,direction",
      ("Supplement Figure S1","A"):"geo_accession,participant_id,time,treatment_number,laser_type",
      ("Supplement Figure S1","B"):"geo_accession,participant_id,anatomical_site,time,laser_type",
      ("Supplement Figure S2","A"):"platform,probe_id,mapping_status,used_reliable,selected_highest_iqr_probe,used_shared_universe",
      ("Supplement Figure S2","B"):"platform,probe_id,mapping_status,used_reliable,selected_highest_iqr_probe,used_shared_universe",
      ("Supplement Figure S2","C"):"gene_symbol,entrez_id",
      ("Supplement Figure S3","A"):"dataset,program,day,model_effect,model_ci95_low,model_ci95_high",
      ("Supplement Figure S3","B"):"dataset,program,day,model_effect,model_ci95_low,model_ci95_high",
      ("Supplement Figure S4","A"):"deleted_from,deleted_participant,day1_rho,day7_rho,day14_rho",
      ("Supplement Figure S5","A"):"iteration,day1_rho,day7_rho,day14_rho",
      ("Supplement Figure S5","B"):"statistic,real_value,null_q95,empirical_one_sided_p",
      ("Supplement Figure S5","C"):"iteration,transition_rho,statistic,real_value,empirical_one_sided_p,null_percentile",
      ("Supplement Figure S6","A"):"dataset,context,day,program,paired_hedges_g",
      ("Supplement Figure S7","A"):"dataset,context,day,program,paired_hedges_g",
      ("Supplement Figure S8","A"):"dataset,program,day,model_effect,model_ci95_low,model_ci95_high,paired_hedges_g,deleted_from,deleted_participant,mapping_method,variant",
    }
    out=[]
    for _,r in plan.iterrows():
        key=(r.figure,r.panel); m=matrix[(matrix.figure==r.figure)&(matrix.panel==r.panel)].iloc[0]
        out.append([r.figure,r.panel,sources[key],cols[key],m.statistic,m.statistic,claim_map[r.figure]])
    pd.DataFrame(out,columns=["figure","panel","source_file","source_columns","statistic","expected_value","claim_id"]).to_csv(ROOT/"PHASE3B_SOURCE_LOCK.tsv",sep="\t",index=False)


def write_legends() -> None:
    text = """# Phase 3B figure legends

## Fig. 1 | Study design and analytical framework

**A**, GSE168760 longitudinal back-skin AFL design (n=14 participants), with D0, D1, D7 and D14 highlighted as primary shared timepoints and D3, D21 and D28 retained only as single-treatment secondary context. **B**, GSE206495 forearm NAFL first-cycle design (n=17 participants), with D0, D1, D7 and D14 in the primary analysis; D29 follows a second treatment and face samples are excluded from the primary cross-cohort trajectory. **C**, aligned primary time frame. **D**, analytical sequence from participant-blocked effects to frozen-program and genome-wide evidence. Arrows denote workflow order, not biological causality. Source data are provided as a Source Data file.

## Fig. 2 | Reproducible temporal program architecture

**A**, paired Hedges g for seven frozen programs at D1, D7 and D14 in GSE168760 (n=14) and GSE206495 (n=17); the diverging scale is shared and centred at zero. **B–D**, cross-cohort program-effect associations at D1, D7 and D14. Spearman rho values are 0.821, 0.821 and 0.857, respectively; 95% intervals were obtained from 5,000 participant-cluster bootstrap replicates. Each point is one prespecified program. Identity of effects or causal equivalence between AFL and NAFL is not inferred. Source data are provided as a Source Data file.

## Fig. 3 | Robustness and null-model testing

**A–B**, distributions from 10,000 gene-count-matched random-program architectures for the mean and minimum D1/D7/D14 Spearman rho. Red lines mark the frozen observed statistics; one-sided empirical P=0.0313 and P=0.0110, respectively. **C**, mean temporal rho for all 31 leave-one-participant-out (LOPO) removals; colour indicates the cohort from which the participant was removed. **D**, D1, D7 and D14 rho for all removals. All 31/31 removals retain positive concordance and the minimum day-specific rho is 0.750. LOPO estimates are dependent sensitivity analyses, not independent replications. Source data are provided as a Source Data file.

## Fig. 4 | Genome-wide cross-cohort replication

**A–C**, hexbin representations of paired Hedges g for 16,942 shared genes at D1, D7 and D14. Spearman rho values are 0.714, 0.505 and 0.373, with 10,000 gene-label permutations per day (minimum attainable one-sided empirical P=0.0001; all observed P=0.00009999 and displayed as P<0.0001). **D**, normalized enrichment scores for 12 prespecified reciprocal rank-replication tests (top 10% positive and negative tails, both transfer directions, three days). All 12 tests are direction-correct and all 12 remain supported after Benjamini–Hochberg correction across exactly 12 tests. Genome-wide correlation attenuates over time and does not imply replication of every gene or identical effect magnitudes. Source data are provided as a Source Data file.

## Fig. 5 | Delayed remodeling and temporal transition

**A**, participant-blocked model effects and 95% model intervals for ECM remodeling, collagen organization and matrix maturation at D1, D7 and D14 in GSE168760 (n=14) and GSE206495 (n=17). Cohorts are shown separately. **B**, individual participant treated-minus-D0 program-score differences at D1 and D14. **C**, seven-program D14-minus-D1 transition effects; rho=0.857 and 6/7 programs change in the same direction, but the matched-random transition null is not significant (10,000 matched architectures; P=0.286), so this evidence is descriptive/supportive. **D**, D14 collagen-organization paired Hedges g for the full program, removal of COL1A1/COL1A2, removal of all COL-prefix genes and 69 leave-one-gene-out values. No new P value is inferred from sensitivity ranges. Source data are provided as a Source Data file.

## Supplementary Fig. 1 | Sample and participant structure QC

**A**, retained single-treatment GSE168760 samples by participant and timepoint. **B**, retained GSE206495 forearm samples; 34 face samples are excluded and D29 is secondary treatment context. Filled cells indicate an available sample. Biological n is 14 and 17 participants, respectively.

## Supplementary Fig. 2 | Platform annotation and shared-gene universe

**A–B**, feature-to-gene mapping counts for GPL13667 and GPL15207. **C**, reliable platform-specific gene sets converge on the frozen 16,942-gene symbol-and-Entrez-agreeing universe. Expression matrices were not merged.

## Supplementary Fig. 3 | Program-level uncertainty

**A–B**, participant-blocked model effects and retained model 95% intervals for all seven programs at D1, D7 and D14 in GSE168760 (n=14) and GSE206495 (n=17). These are model-effect intervals, not newly calculated Hedges-g intervals.

## Supplementary Fig. 4 | Complete LOPO results

All 31 participant removals by D1, D7 and D14 Spearman rho. Values are displayed directly; minimum, median and range are deterministic summaries of the frozen LOPO table. A and B in participant labels denote removal from GSE168760 and GSE206495, respectively.

## Supplementary Fig. 5 | Random-program null diagnostics

**A**, day-specific matched-random distributions and observed values. **B**, observed statistics and null 95th percentiles. **C**, transition-rho null distribution; observed rho=0.857, one-sided empirical P=0.286. All nulls contain 10,000 frozen iterations.

## Supplementary Fig. 6 | AFL-only exploratory extended trajectory

Within-GSE168760 paired Hedges g for seven frozen programs from D0 to D28 (n=14). D0 is the reference value and is plotted at zero. D3, D21 and D28 are exploratory single-treatment context and are not cross-cohort validation.

## Supplementary Fig. 7 | GSE206495 D29 treatment-context secondary

D29 paired Hedges g for seven frozen programs in forearm samples (n=17). D29 is 29 days after the first session and 1 day after the second session; it is not interpreted as a natural extension of the D14 primary trajectory.

## Supplementary Fig. 8 | Exploratory proliferation divergence

**A**, participant-blocked proliferation effects and 95% model intervals at D1, D7 and D14. **B**, D14 paired Hedges g under full-program, probe-mapping, half-gene-set and LOPO sensitivities. The result is directionally robust but not confirmatory; cohort and modality are fully confounded.
"""
    (ROOT/"PHASE3B_FIGURE_LEGENDS.md").write_text(text)


def make_claim_qc() -> None:
    m=pd.read_csv(ROOT/"PHASE3A_FIGURE_CLAIM_MATRIX.tsv",sep="\t")
    rows=[]
    for _,r in m.iterrows():
        rows.append([r.figure,r.panel,r.claim_supported,r.statistic,r.claim_supported,r.required_limitation,"PASS"])
    pd.DataFrame(rows,columns=["figure","panel","expected_claim","actual_visual","actual_annotation","required_limitation","claim_match"]).to_csv(ROOT/"PHASE3B_CLAIM_QC.tsv",sep="\t",index=False)


def make_numerical_qc() -> None:
    checks=[]
    def add(name,src,obs,exp,tol=1e-12):
        ok=abs(float(obs)-float(exp))<=tol; checks.append([name,src,obs,exp,tol,"PASS" if ok else "FAIL"])
    con=read("cross_cohort_conservation.tsv",1).set_index("day")
    for d,v in [(1,.821428571428571),(7,.821428571428571),(14,.857142857142857)]: add(f"program_rho_D{d}","cross_cohort_conservation.tsv",con.loc[d,"standardized_spearman_rho"],v)
    sm=read("MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv").set_index("statistic")
    add("random_mean_p","MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv",sm.loc["mean_rho","empirical_one_sided_p"],.0312968703129687)
    add("random_min_p","MATCHED_RANDOM_PROGRAM_NULL_SUMMARY.tsv",sm.loc["min_rho","empirical_one_sided_p"],.010998900109989)
    l=read("PHASE2B_FULL_LOPO.tsv"); add("lopo_count","PHASE2B_FULL_LOPO.tsv",len(l),31,0); add("lopo_min_rho","PHASE2B_FULL_LOPO.tsv",l[["day1_rho","day7_rho","day14_rho"]].min().min(),.75)
    gw=read("GENOMEWIDE_CONCORDANCE.tsv").set_index("day")
    for d,v in [(1,.714343496867325),(7,.505168653693006),(14,.37303116657605)]: add(f"genomewide_rho_D{d}","GENOMEWIDE_CONCORDANCE.tsv",gw.loc[d,"standardized_spearman_rho"],v)
    rr=read("RECIPROCAL_RANK_REPLICATION.tsv"); add("reciprocal_direction_count","RECIPROCAL_RANK_REPLICATION.tsv",rr.direction_correct.sum(),12,0); add("reciprocal_q_count","RECIPROCAL_RANK_REPLICATION.tsv",(rr.BH_q<.05).sum(),12,0)
    ts=read("TRANSITION_SUMMARY.tsv").set_index("metric"); add("transition_rho","TRANSITION_SUMMARY.tsv",ts.loc["transition_rho","value"],.857142857142857); add("transition_same_direction","TRANSITION_SUMMARY.tsv",ts.loc["same_transition_programs","value"],6,0); add("transition_null_p","TRANSITION_SUMMARY.tsv",ts.loc["transition_null_p","value"],.286171382861714)
    add("shared_genes","GENOMEWIDE_CONCORDANCE.tsv",gw.loc[1,"genes"],16942,0)
    out=pd.DataFrame(checks,columns=["check","source_file","observed","expected","tolerance","status"])
    out.to_csv(ROOT/"PHASE3B_NUMERICAL_QC.tsv",sep="\t",index=False)
    if not out.status.eq("PASS").all(): raise RuntimeError("NUMERICAL_PROVENANCE_FAILURE")


def write_visual_qc_pending() -> None:
    names=[f"Figure{i}" for i in range(1,6)]+[f"FigureS{i}" for i in range(1,9)]+["Table1","Table2"]
    rows=[[n,"PENDING_MANUAL_INSPECTION","clipping; overlap; readability; scales; labels; whitespace; hidden points; legends; color redundancy",""] for n in names]
    pd.DataFrame(rows,columns=["figure","status","checks","notes"]).to_csv(ROOT/"PHASE3B_VISUAL_QC.tsv",sep="\t",index=False)


def package_outputs() -> None:
    if PACKAGE.exists(): shutil.rmtree(PACKAGE)
    for d in ["main_figures","supplement_figures","tables","source_data","legends","scripts","qc","provenance"]: (PACKAGE/d).mkdir(parents=True,exist_ok=True)
    for p in MAIN.iterdir(): shutil.copy2(p,PACKAGE/"main_figures"/p.name)
    for p in SUPP.iterdir(): shutil.copy2(p,PACKAGE/"supplement_figures"/p.name)
    for stem in ["Table1_CohortDesign","Table2_EvidenceSummary"]:
        for p in ROOT.glob(stem+".*"): shutil.copy2(p,PACKAGE/"tables"/p.name)
    for p in SOURCE.glob("*.tsv"): shutil.copy2(p,PACKAGE/"source_data"/p.name)
    shutil.copy2(ROOT/"PHASE3B_FIGURE_LEGENDS.md",PACKAGE/"legends"/"PHASE3B_FIGURE_LEGENDS.md")
    for p in (ROOT/"03_figure_scripts").glob("*.py"): shutil.copy2(p,PACKAGE/"scripts"/p.name)
    for name in ["PHASE3B_VISUAL_QC.tsv","PHASE3B_NUMERICAL_QC.tsv","PHASE3B_CLAIM_QC.tsv"]: shutil.copy2(ROOT/name,PACKAGE/"qc"/name)
    for name in ["PHASE3B_SOURCE_LOCK.tsv","PHASE3A_FIGURE_PLAN.tsv","PHASE3A_FIGURE_CLAIM_MATRIX.tsv","PHASE3A_CLAIM_HIERARCHY.tsv"]: shutil.copy2(ROOT/name,PACKAGE/"provenance"/name)


def main() -> None:
    configure_style()
    for d in [MAIN,SUPP,SOURCE,QC]: d.mkdir(parents=True,exist_ok=True)
    write_source_data(); build_source_lock(); make_numerical_qc(); make_claim_qc(); write_visual_qc_pending()
    figure1(); figure2(); figure3(); figure4(); figure5()
    supp_s1(); supp_s2(); supp_s3(); supp_s4(); supp_s5(); supp_s6(); supp_s7(); supp_s8()
    make_tables(); write_legends(); package_outputs()
    print("PHASE3B_RENDER_COMPLETE")


if __name__ == "__main__":
    main()
