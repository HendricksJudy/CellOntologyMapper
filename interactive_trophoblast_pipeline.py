#!/usr/bin/env python3
"""
Robust end-to-end pipeline for single-cell trophoblast data (interactive)
"""
from __future__ import annotations

import argparse
import logging
import os
import shutil
from pathlib import Path
from typing import Sequence

import matplotlib.colors as mcolors
import scanpy as sc
import omicverse as ov
from colormaps.colormap import Colormap
import colormaps as cmaps

logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.INFO)
LOGGER = logging.getLogger(__name__)


def rgba_to_hex(colors: Sequence[Sequence[float]], include_alpha: bool = False) -> list[str]:
    return [mcolors.to_hex(rgba, keep_alpha=include_alpha) for rgba in colors]


def safe_split(label: str, sep: str = "_", idx: int = 1) -> str:
    tokens = label.split(sep, maxsplit=idx)
    if len(tokens) <= idx:
        LOGGER.warning("Label '%s' lacks separator '%s'; returning unchanged", label, sep)
        return label
    return tokens[idx]


def ensure_listed_colormap(cmap: Colormap, n: int) -> Colormap:
    if not hasattr(cmap, "colors"):
        LOGGER.info("Resampling continuous colormap '%s' to %d discrete colours", cmap.name, n)
        return cmap.resampled(n)
    return cmap


def palette_from_categories(categories: Sequence[str], base_cmap: Colormap) -> list[str]:
    cats = sorted(categories)
    base_cmap = ensure_listed_colormap(base_cmap, len(cats))
    hex_codes = rgba_to_hex(base_cmap.colors)
    if len(hex_codes) < len(cats):
        raise ValueError("Colour map shorter than number of categories")
    return [hex_codes[cats.index(c)] for c in categories]


def validate_model_dir(p: Path) -> Path:
    p = p.expanduser().resolve()
    if p.exists() and p.is_file():
        msg = f"model_dir path '{p}' exists but is a file."
        if not os.isatty(0):
            raise NotADirectoryError(msg)
        print(msg)
        choice = input("Rename to *.bak (r), delete (d), or abort (a)? [r/d/a]: ").strip().lower()
        if choice == "r":
            bak = p.with_suffix(p.suffix + ".bak")
            shutil.move(p, bak)
            LOGGER.info("Renamed offending file to %s", bak)
        elif choice == "d":
            p.unlink()
            LOGGER.info("Deleted offending file")
        else:
            raise SystemExit("Aborted by user.")
    p.mkdir(parents=True, exist_ok=True)
    return p


def load_adata(path: Path) -> sc.AnnData:
    adata = sc.read(path)
    if "X_umap" not in adata.obsm:
        raise KeyError("UMAP embedding 'X_umap' not found in adata.obsm")
    return adata


def clean_celltype_predictions(adata: sc.AnnData) -> None:
    col = "celltype_predictions"
    if col not in adata.obs:
        raise KeyError(f"Expected column '{col}' in adata.obs")
    adata.obs[col] = adata.obs[col].astype(str).map(safe_split).fillna("Unknown").astype("category")


def plot_umap(
    adata: sc.AnnData,
    colour_key: str,
    title: str,
    palette: list[str] | dict[str, str],
    out_path: Path,
    size: int = 10,
) -> None:
    fig, ax = ov.plt.subplots(figsize=(3, 3))
    ov.pl.embedding(
        adata,
        basis="X_umap",
        color=colour_key,
        title=title,
        palette=palette,
        show=False,
        ax=ax,
        size=size,
    )
    ax.xaxis.set_label_coords(0.05, 0.04)
    ax.yaxis.set_label_coords(0.04, 0.05)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=300, bbox_inches="tight")
    LOGGER.info("Saved %s", out_path)


def run_mapper(
    adata: sc.AnnData,
    cl_json: Path,
    taxonomy_txt: Path,
    model_dir: Path,
    tissue: str = "trophoblast",
    species: str = "Homo sapiens",
    threshold: float = 0.3,
) -> sc.AnnData:
    model_dir = validate_model_dir(model_dir)
    mapper = ov.single.CellOntologyMapper(
        cl_obo_file=cl_json,
        model_name="sentence-transformers/all-MiniLM-L6-v2",
        local_model_dir=model_dir,
    )
    mapper.setup_llm_expansion(
        api_type="ollama",
        model="gemma3n",
        base_url="http://localhost:11434",
        tissue_context=tissue,
        species=species,
        study_context=f"{tissue} development in early pregnancy",
    )
    mapper.load_cell_taxonomy_resource(taxonomy_txt, species_filter=[species, "Mus musculus"])
    mapper.map_adata_with_taxonomy(
        adata,
        cell_name_col="celltype_predictions",
        new_col_name="enhanced_cell_ontology",
        expand_abbreviations=True,
        use_taxonomy=True,
        species=species,
        tissue_context=tissue,
        threshold=threshold,
    )
    expected = {
        "enhanced_cell_ontology",
        "enhanced_cell_ontology_taxonomy_match",
        "enhanced_cell_ontology_ct_id",
        "enhanced_cell_ontology_cl_id",
    }
    if expected - set(adata.obs.columns):
        raise RuntimeError("Mapper did not create expected columns")
    return adata


def prompt_path(question: str, must_exist: bool = True) -> Path:
    p = Path(input(f"{question}: ").strip()).expanduser()
    if must_exist and not p.exists():
        raise FileNotFoundError(p)
    return p


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(description="Reliable single-cell pipeline (interactive)")
    parser.add_argument("h5ad", nargs="?", type=Path)
    parser.add_argument("cl_json", nargs="?", type=Path)
    parser.add_argument("taxonomy", nargs="?", type=Path)
    parser.add_argument("model_dir", nargs="?", type=Path)
    parser.add_argument("--figdir", type=Path, default=Path("figures"))
    args = parser.parse_args(argv)
    if not all((args.h5ad, args.cl_json, args.taxonomy, args.model_dir)):
        LOGGER.warning("Missing CLI arguments – entering interactive mode.")
        if args.h5ad is None:
            args.h5ad = prompt_path("Path to .h5ad file")
        if args.cl_json is None:
            args.cl_json = prompt_path("Path to Cell Ontology JSON")
        if args.taxonomy is None:
            args.taxonomy = prompt_path("Path to Cell-Taxonomy TXT")
        if args.model_dir is None:
            args.model_dir = prompt_path("Directory for sentence-transformer cache", must_exist=False)
    adata = load_adata(args.h5ad)
    LOGGER.info("Loaded AnnData with %d cells and %d variables", adata.n_obs, adata.n_vars)
    clean_celltype_predictions(adata)
    categories = list(adata.obs["celltype_predictions"].cat.categories)
    colour_list = palette_from_categories(categories, cmaps.cet_g_bw_minc_minl)
    adata.uns["celltype_predictions_colors"] = colour_list
    plot_umap(
        adata,
        "celltype_predictions",
        "Trophoblast (author labels)",
        colour_list,
        args.figdir / "umap_trophoblast_author.png",
        size=6,
    )
    adata = run_mapper(adata, args.cl_json, args.taxonomy, args.model_dir)
    map_targets = {
        "enhanced_cell_ontology": "Cell Ontology term",
        "enhanced_cell_ontology_taxonomy_match": "Cell Taxonomy term",
        "enhanced_cell_ontology_ct_id": "Cell Taxonomy ID",
        "enhanced_cell_ontology_cl_id": "Cell Ontology ID",
    }
    for key, label in map_targets.items():
        mapping = (
            adata.obs[["celltype_predictions", key]]
            .drop_duplicates()
            .set_index("celltype_predictions")[key]
            .to_dict()
        )
        palette = {
            mapping[k]: adata.uns["celltype_predictions_colors"][i]
            for i, k in enumerate(mapping.keys())
            if k in mapping
        }
        plot_umap(
            adata,
            key,
            f"Trophoblast – {label}",
            palette,
            args.figdir / f"umap_trophoblast_{key}.png",
            size=6,
        )
    out_csv = args.h5ad.with_name(args.h5ad.stem + "_annotated.csv")
    adata.obs.to_csv(out_csv)
    LOGGER.info("Wrote annotated metadata to %s", out_csv)


if __name__ == "__main__":
    main()
