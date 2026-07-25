#!/usr/bin/env python3
"""Build concise supplementary-figure Source Data Excel workbooks."""

from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "resource_data"
OUT.mkdir(exist_ok=True)

# Workbook name -> (sheet, relative TSV, core columns, reader-facing description)
FIGURES = {
    "FigureS1_resource_data.xlsx": [
        (
            "S1a",
            "output/figs1/figs1a-clinical_summary.tsv",
            ["patient", "response", "response_v2", "treatment_sum", "mandard_score"],
            "Patient-level clinical information used for the summary bars in Supplementary Figure 1a.",
        ),
    ],
    "FigureS2_resource_data.xlsx": [
        (
            "S2a",
            "output/figs2/figs2a-signature_contribution.tsv",
            ["sample", "sig_cat", "prop"],
            "Mutational signature contributions shown for each sample in Supplementary Figure 2a.",
        ),
        (
            "S2b",
            "output/figs2/figs2b-immune_editing_score-tn_pair-fresh-t.tsv",
            ["patient_id", "response", "sample_id", "patient", "sample_type", "editing_score"],
            "Immune-editing scores and sample groups shown in Supplementary Figure 2b.",
        ),
    ],
    "FigureS3_resource_data.xlsx": [
        (
            "S3a",
            "output/figs3/figs3a-epi-malig_ratio.tsv",
            ["cell_type", "malig", "non_malig", "epi_major"],
            "Malignant and non-malignant epithelial-cell counts shown in Supplementary Figure 3a.",
        ),
        (
            "S3b",
            "output/figs3/figs3b-epi-cluster_cc-flt.tsv",
            ["patient", "patient_id", "response", "sample", "sample_type", "cell_type", "freq", "pct"],
            "Epithelial-cell cluster counts and percentages shown in Supplementary Figure 3b.",
        ),
        (
            "S3c_interactions",
            "output/figs3/figs3c-cellchat_individual-flt_comm.tsv",
            ["response", "sample_type", "source", "target", "interaction_name", "patient", "ligand", "receptor", "prob", "pval", "n_sample"],
            "Individual ligand-receptor interactions retained for the Supplementary Figure 3c analysis.",
        ),
        (
            "S3c_summary",
            "output/figs3/figs3c-cellchat_individual-flt_comm-stat.tsv",
            ["response", "sample_type", "source", "target", "n_interaction", "prob_mean", "prob_med"],
            "Interaction counts and mean or median probabilities used for Supplementary Figure 3c.",
        ),
        (
            "S3c_nodes",
            "output/figs3/figs3c-cellchat_individual-node_info.tsv",
            ["name", "node_color"],
            "Node labels and colors used in the Supplementary Figure 3c network.",
        ),
        (
            "S3c_edges",
            "output/figs3/figs3c-cellchat_individual-prob_mean_diff-R_vs_NR.tsv",
            ["source", "target", "weight.x", "weight.y", "weight.diff", "weight", "comp_type", "weight_type"],
            "R-versus-NR edge weights plotted in the Supplementary Figure 3c network.",
        ),
    ],
    "FigureS4_resource_data.xlsx": [
        (
            "S4a",
            "output/figs4/figs4a-virus_detection_in_macrotype-bi-cell_cnt.tsv",
            ["celltype", "bi_type", "n_cell"],
            "Virus-detected and non-detected cell counts shown in Supplementary Figure 4a.",
        ),
        (
            "S4b",
            "output/figs4/figs4b-virus_detection_in_macrotype-bi-sample_lvl-pct.tsv",
            ["patient", "response", "sample", "sample_type", "celltype", "Non-detected", "Virus-detected"],
            "Sample-level percentages of virus-detected cells shown in Supplementary Figure 4b.",
        ),
    ],
    "FigureS5_resource_data.xlsx": [
        (
            "S5a",
            "output/figs5/figs5a-exhau_cyto_score.tsv",
            ["celltype", "cell_state", "ex_score", "cyto_score"],
            "Cell-level exhaustion and cytotoxicity scores shown in Supplementary Figure 5a.",
        ),
        (
            "S5d",
            "output/figs5/figs5d-pbulk_exprs_in_t_state-CXCL13.tsv",
            ["patient", "response", "sample_type", "cellgp", "gene", "logcpm"],
            "CXCL13 pseudobulk expression values shown in Supplementary Figure 5d.",
        ),
        (
            "S5e",
            "output/figs5/figs5e-selected_t_cluster_composition.tsv",
            ["patient", "response", "sample", "sample_type", "cell_type", "pct"],
            "Selected T-cell cluster percentages shown in Supplementary Figure 5e.",
        ),
        (
            "S5f",
            "output/figs5/figs5f-t_expand_in_subtype-min10.tsv",
            ["patient", "response", "sample_type", "subtype", "n_expand_cell_per_sample_stype", "n_cell_per_sample_subtype", "pct_by_subtype"],
            "Expanded-cell counts and percentages shown in Supplementary Figure 5f.",
        ),
    ],
    "FigureS6_resource_data.xlsx": [
        (
            "S6a",
            "output/figs6/figs6a-pbulk_exprs_in_t_state-chemokine_receptors.tsv",
            ["patient", "response", "sample_type", "cellgp", "gene", "logcpm"],
            "Chemokine-receptor pseudobulk expression values shown in Supplementary Figure 6a.",
        ),
        (
            "S6b",
            "output/figs6/figs6b-pbulk_exprs_in_whole-chemokine_ligands.tsv",
            ["patient", "response", "sample", "sample_type", "symbol", "logcpm"],
            "Chemokine-ligand pseudobulk expression values shown in Supplementary Figure 6b.",
        ),
        (
            "S6c",
            "output/figs6/figs6c-pbulk_exprs_in_major-CXCL12.tsv",
            ["patient", "response", "sample", "sample_type", "celltype", "symbol", "logcpm"],
            "CXCL12 pseudobulk expression values shown in Supplementary Figure 6c.",
        ),
    ],
    "FigureS7_resource_data.xlsx": [
        (
            "S7a_edges",
            "output/figs7/figs7a-edge_info-r_vs_nr.tsv",
            ["edge_type", "from", "to", "start_cyto", "start_exhau", "end_cyto", "end_exhau", "weight.x", "weight.y", "weight.diff", "weight", "wtype", "comp_type"],
            "R-versus-NR edge weights plotted in the Supplementary Figure 7a network.",
        ),
        (
            "S7a_nodes",
            "output/figs7/figs7a-node_info-source_data.tsv",
            ["node", "exhau", "cyto", "node_type"],
            "Node positions and cell-state labels used in the Supplementary Figure 7a network.",
        ),
    ],
    "FigureS8_resource_data.xlsx": [
        (
            "S8a_expansion",
            "output/figs8/figs8a-expanded_pct_delta-CD8_CXCL13.tsv",
            ["patient", "response", "subtype", "Baseline", "Treat", "pct_expand_diff"],
            "Change in the expanded CD8_CXCL13-cell percentage shown in Supplementary Figure 8a.",
        ),
        (
            "S8a_diversity",
            "output/figs8/figs8a-diversity-CD8_CXCL13-min10.tsv",
            ["patient", "response", "sample_type", "shannon"],
            "CD8_CXCL13 TCR diversity values shown in Supplementary Figure 8a.",
        ),
        (
            "S8b",
            "output/figs8/figs8b-tcr_share_inter_subtype-clone_frac_new.tsv",
            ["sample", "patient", "sample_type", "subtype", "clonotype", "share_type", "clone_frac_new", "cell_state"],
            "Clonotype fractions and sharing categories shown in Supplementary Figure 8b.",
        ),
        (
            "S8e1",
            "output/figs8/figs8e1-pbulk_exprs-weighted_CX3CL1_in_endo.tsv",
            ["patient", "response", "sample", "sample_type", "wt_cx3cl1_exprs_unlog", "logcpm"],
            "Weighted endothelial CX3CL1 expression values shown in Supplementary Figure 8e-1.",
        ),
        (
            "S8e2",
            "output/figs8/figs8e2-pbulk_exprs-cx3cl1_in_TME.tsv",
            ["patient", "response", "sample", "sample_type", "cx3cl1_exprs_unlog", "cx3cl1_exprs_log"],
            "Tumor-microenvironment CX3CL1 expression values shown in Supplementary Figure 8e-2.",
        ),
        (
            "S8f",
            "output/figs8/figs8f-weighted_CX3CL1_exprs-CD8_CX3CR1_diversity.tsv",
            ["patient", "response", "sample", "sample_type", "wt_cx3cl1_exprs_log1p", "cd8_cx3cr1_shannon"],
            "CX3CL1 expression and CD8_CX3CR1 TCR diversity values shown in Supplementary Figure 8f.",
        ),
        (
            "S8g",
            "output/figs8/figs8g-weighted_CX3CL1_exprs-CD8_CX3CR1_expand_pct_in_TNK.tsv",
            ["patient", "response", "sample", "wt_cx3cl1_exprs_log1p", "cd8_cx3cr1_expand_pct_in_TNK"],
            "CX3CL1 expression and expanded CD8_CX3CR1-cell percentages shown in Supplementary Figure 8g.",
        ),
        (
            "S8h",
            "output/figs8/figs8h-cc_endo.tsv",
            ["patient", "response", "sample", "sample_type", "cell_type", "pct"],
            "Endothelial-cell subtype percentages shown in Supplementary Figure 8h.",
        ),
    ],
}


def normalise_labels(df):
    """Use current binary-response terminology without changing detailed response fields."""
    for col in ("response", "comp_type", "weight_type", "wtype"):
        if col in df:
            present = df[col].notna()
            df.loc[present, col] = (
                df.loc[present, col]
                .replace({"PR": "NR"})
                .astype(str)
                .str.replace("PR", "NR", regex=False)
            )
    return df


def preserve_patient_ids(df):
    """Keep three-digit patient IDs as text in Excel (for example, 006)."""
    if "patient_id" in df:
        values = pd.to_numeric(df["patient_id"], errors="coerce")
        if values.notna().all():
            df["patient_id"] = values.astype("Int64").astype("string").str.zfill(3)
        else:
            df["patient_id"] = df["patient_id"].astype("string")
    return df


for workbook, tables in FIGURES.items():
    contents = []
    with pd.ExcelWriter(OUT / workbook, engine="openpyxl") as writer:
        for sheet, rel, columns, description in tables:
            path = ROOT / rel
            if not path.exists():
                raise FileNotFoundError(path)
            df = pd.read_csv(path, sep="\t")
            missing = [column for column in columns if column not in df.columns]
            if missing:
                raise ValueError(f"{path}: missing columns {missing}")
            df = normalise_labels(df.loc[:, columns].copy())
            df = preserve_patient_ids(df)
            df.to_excel(writer, sheet_name=sheet, index=False)
            if "patient_id" in df:
                worksheet = writer.book[sheet]
                patient_id_column = list(df.columns).index("patient_id") + 1
                for cells in worksheet.iter_cols(
                    min_col=patient_id_column,
                    max_col=patient_id_column,
                    min_row=2,
                ):
                    for cell in cells:
                        cell.number_format = "@"
            contents.append((sheet, description))

        pd.DataFrame(contents, columns=["sheet", "description"]).to_excel(
            writer, sheet_name="Contents", index=False
        )

        workbook_object = writer.book
        workbook_object._sheets.insert(
            0, workbook_object._sheets.pop(workbook_object.sheetnames.index("Contents"))
        )
        for worksheet in workbook_object.worksheets:
            worksheet.freeze_panes = "A2"
            worksheet.auto_filter.ref = worksheet.dimensions
            for column in worksheet.columns:
                letter = column[0].column_letter
                width = min(72, max(12, max(len(str(cell.value or "")) for cell in column) + 2))
                worksheet.column_dimensions[letter].width = width

print("Created", len(FIGURES), "workbooks in", OUT)
