#!/usr/bin/env python3

import argparse
import os
import sys
import tempfile

import muon as mu
import kaichi


MODELS = ["umi", "max", "ratio", "poisson_gauss"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="kaichi guide assignment")
    parser.add_argument(
        "--rawcounts.h5mu",
        dest="rawcounts",
        required=True,
        help="path to input h5mu (MuData with crispr modality)",
    )
    parser.add_argument("--output_dir", "-o", required=True, help="output directory")
    parser.add_argument("--name", "-n", required=True, help="output file prefix")
    parser.add_argument(
        "--model",
        required=True,
        choices=MODELS,
        help="assignment model",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    os.makedirs(args.output_dir, exist_ok=True)

    print("Loading input h5mu ...", file=sys.stderr)
    mdata = mu.read_h5mu(args.rawcounts)
    adata_guides = mdata["crispr"]
    print(
        f"  crispr modality: {adata_guides.n_obs} x {adata_guides.n_vars}",
        file=sys.stderr,
    )

    if adata_guides.n_vars == 0:
        raise ValueError("crispr modality is empty")

    # kaichi.assign() takes an h5ad path; write the crispr modality to a temp h5ad.
    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        adata_guides.write_h5ad(tmp_path)
        print(f"Running kaichi.assign(model={args.model!r}) ...", file=sys.stderr)
        out = kaichi.assign(tmp_path, model=args.model)
    finally:
        os.unlink(tmp_path)

    # kaichi returns a minimal .var (guide IDs only). Re-attach input's var
    # metadata (target_gene, sequence, etc.) so the output matches crispat's shape.
    print("Merging input .var metadata ...", file=sys.stderr)
    out.var = adata_guides.var.loc[out.var_names].copy()
    out.uns["guide_assignment_params"] = {"model": args.model}

    # Preserve input obs columns (e.g. batch) not written by kaichi.
    input_extra = adata_guides.obs.drop(
        columns=[c for c in adata_guides.obs.columns if c in out.obs.columns],
        errors="ignore",
    )
    out.obs = out.obs.join(input_extra)

    # kaichi calls this n_guides_detected; rename to n_guides_assigned for consistency.
    if "n_guides_detected" in out.obs.columns:
        out.obs = out.obs.rename(columns={"n_guides_detected": "n_guides_assigned"})

    # target_gene for singly-assigned cells (multi-infected and unassigned get "").
    if "target_gene" in out.var.columns:
        guide_to_gene = out.var["target_gene"].to_dict()
        single = ~out.obs["is_unassigned"].astype(bool) & ~out.obs[
            "is_multi_infected"
        ].astype(bool)
        out.obs["target_gene"] = ""
        out.obs.loc[single, "target_gene"] = (
            out.obs.loc[single, "guide_id"].map(guide_to_gene).fillna("")
        )

    n_assigned = int((~out.obs["is_unassigned"].astype(bool)).sum())
    n_total = out.n_obs
    print(
        f"  {n_assigned}/{n_total} cells assigned ({100 * n_assigned / n_total:.1f}%)",
        file=sys.stderr,
    )

    out_path = os.path.join(args.output_dir, f"{args.name}.h5ad")
    print(f"Writing {out_path} ...", file=sys.stderr)
    out.write_h5ad(out_path)


if __name__ == "__main__":
    main()
