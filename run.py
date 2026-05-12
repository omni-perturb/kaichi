#!/usr/bin/env python3

import argparse
import os
import subprocess
import sys
import tempfile

import anndata as ad
import muon as mu


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

    out_path = os.path.join(args.output_dir, f"{args.name}.h5ad")

    with tempfile.NamedTemporaryFile(suffix=".h5ad", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        adata_guides.write_h5ad(tmp_path)
        print(f"Running kaichi {args.model} ...", file=sys.stderr)
        subprocess.run(
            [
                "kaichi",
                "assign",
                "--counts",
                tmp_path,
                "--model",
                args.model,
                "--output",
                out_path,
            ],
            check=True,
        )
    finally:
        os.unlink(tmp_path)

    # kaichi writes a minimal .var (guide IDs only). Re-attach input's var
    # metadata (target_gene, sequence, etc.) so the output matches crispat's shape.
    print("Merging input .var metadata ...", file=sys.stderr)
    out = ad.read_h5ad(out_path)
    out.var = adata_guides.var.loc[out.var_names].copy()
    out.uns["guide_assignment_params"] = {"model": args.model}

    n_assigned = int((out.obs["guide_identity"].astype(str) != "").sum())
    n_total = out.n_obs
    print(
        f"  {n_assigned}/{n_total} cells assigned ({100 * n_assigned / n_total:.1f}%)",
        file=sys.stderr,
    )

    print(f"Writing {out_path} ...", file=sys.stderr)
    out.write_h5ad(out_path)


if __name__ == "__main__":
    main()
