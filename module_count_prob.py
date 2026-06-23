#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Permutation test for multilayer network modularity using Infomap.

Accepts standard Infomap *Multilayer format directly:

  # A multilayer network
  *Vertices N
  id "name"
  ...
  *Multilayer
  # layer_id node_id layer_id node_id weight
  18 1338 18 1346 0.009477332
  ...

Permutation logic:
  Treat (layer_v, node_v) -- columns 3 & 4 -- as an atomic pair.
  Shuffle these pairs across ALL rows (both intra- and inter-layer edges).
  layer_u, node_u, and weight stay in their original rows.
  Each node always drags its own layer_v with it when reassigned.

Infomap command per run:
  ./Infomap network . --ftree --clu -f undirected -s 3 -N 100

Usage:
  python infomap_perm_nodes.py \\
    --input network.txt \\
    --n_perm 1000 \\
    --outdir results/ \\
    --infomap ./Infomap
"""

import argparse
import os
import random
import subprocess
import tempfile
import numpy as np
from collections import Counter


# -- CLI -----------------------------------------------------------------------
def get_args():
    p = argparse.ArgumentParser()
    p.add_argument("--input",   required=True,
                   help="Infomap *Multilayer .txt file")
    p.add_argument("--n_perm",  type=int, default=1000)
    p.add_argument("--outdir",  default="perm_results")
    p.add_argument("--infomap", default="./Infomap",
                   help="Path to Infomap binary (default: ./Infomap)")
    p.add_argument("--seed",    type=int, default=42)
    return p.parse_args()


# -- PARSER --------------------------------------------------------------------
def read_multilayer_txt(path):
    """
    Parse an Infomap *Multilayer file.

    Sections handled:
      *Vertices  -> node id + quoted name
      *Multilayer (or *Intra / *Inter for legacy format) -> edges

    Returns:
        nodes : dict {node_id (int) -> name (str)}
        edges : list of (layer_u, node_u, layer_v, node_v, weight)
    """
    nodes   = {}
    edges   = []
    section = None

    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            low = line.lower()

            # -- section headers ----------------------------------------------
            if low.startswith("*vertices"):
                section = "vertices"
                continue
            if low.startswith("*multilayer"):
                section = "multilayer"
                continue
            if low.startswith("*intra"):
                section = "intra"
                continue
            if low.startswith("*inter"):
                section = "inter"
                continue

            # -- vertices -----------------------------------------------------
            if section == "vertices":
                parts = line.split(None, 1)
                nid  = int(parts[0])
                name = parts[1].strip('"') if len(parts) > 1 else str(nid)
                nodes[nid] = name

            # -- multilayer edges: layer_u node_u layer_v node_v weight -------
            elif section in ("multilayer", "intra", "inter"):
                parts = line.split()
                if len(parts) < 5:
                    continue
                lu = int(parts[0]);  u = int(parts[1])
                lv = int(parts[2]);  v = int(parts[3])
                w  = float(parts[4])
                edges.append((lu, u, lv, v, w))

    return nodes, edges


def write_multilayer(path, nodes, edges):
    """Write the network back in Infomap *Multilayer format."""
    with open(path, "w") as fh:
        fh.write("# multilayer network (permuted)\n")
        fh.write(f"*Vertices {len(nodes)}\n")
        for nid in sorted(nodes):
            fh.write(f'{nid} "{nodes[nid]}"\n')
        fh.write("*Multilayer\n")
        fh.write("# layer_u node_u layer_v node_v weight\n")
        for (lu, u, lv, v, w) in edges:
            fh.write(f"{lu} {u} {lv} {v} {w}\n")


# -- PERMUTATION ---------------------------------------------------------------
def permute_target_pairs(edges, rng):
    """
    Shuffle (layer_v, node_v) pairs across all rows.
    layer_u, node_u, and weight remain in their original rows.

    This rewires who each source node interacts with while every node
    keeps its own layer label.
    """
    target_pairs = [(lv, v) for (_, _, lv, v, _) in edges]
    rng.shuffle(target_pairs)

    new_edges = []
    for (lu, u, _, _, w), (lv_new, v_new) in zip(edges, target_pairs):
        new_edges.append((lu, u, lv_new, v_new, w))
    return new_edges


# -- INFOMAP -------------------------------------------------------------------
def run_infomap(net_path, outdir, infomap_bin):
    """Run Infomap; return number of top-level modules or None on failure."""
    cmd = [
        infomap_bin, net_path, outdir,
        "--ftree", "--clu", "-f", "undirected", "-s", "3", "-N", "100",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("  [Infomap error]", result.stderr[:400])
        return None

    base = os.path.splitext(os.path.basename(net_path))[0]
    for ext in (".ftree", ".tree"):
        tree = os.path.join(outdir, base + ext)
        if os.path.exists(tree):
            break
    else:
        print(f"  [Warning] .ftree/.tree not found for {net_path}")
        return None

    modules = set()
    with open(tree) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("*"):
                continue
            top = line.split()[0].split(":")[0]
            modules.add(top)
    return len(modules)


# -- MAIN ----------------------------------------------------------------------
def main():
    args = get_args()
    rng  = random.Random(args.seed)
    os.makedirs(args.outdir, exist_ok=True)

    # -- Parse network ---------------------------------------------------------
    print("=== Parsing network ===")
    nodes, edges = read_multilayer_txt(args.input)
    n_intra = sum(1 for (lu, _, lv, _, _) in edges if lu == lv)
    n_inter = len(edges) - n_intra
    print(f"  {len(nodes)} nodes | "
          f"{n_intra} intra-layer edges | "
          f"{n_inter} inter-layer edges | "
          f"{len(edges)} total edges")

    lv_counts = Counter(lv for (_, _, lv, _, _) in edges)
    print(f"  {len(lv_counts)} distinct target layers: "
          + ", ".join(f"L{l}={c}" for l, c in sorted(lv_counts.items())))

    # -- Observed run ----------------------------------------------------------
    print("\n=== Running Infomap on observed network ===")
    obs_dir = os.path.join(args.outdir, "observed")
    os.makedirs(obs_dir, exist_ok=True)
    obs_net = os.path.join(obs_dir, "observed.txt")
    write_multilayer(obs_net, nodes, edges)

    n_obs = run_infomap(obs_net, obs_dir, args.infomap)
    if n_obs is None:
        raise RuntimeError("Infomap failed on the observed network. "
                           "Check the binary path / network format.")
    print(f"  Observed modules: {n_obs}")

    # -- Permutation runs ------------------------------------------------------
    print(f"\n=== Running {args.n_perm} permutations ===")
    perm_counts = []
    perm_dir = os.path.join(args.outdir, "permutations")
    os.makedirs(perm_dir, exist_ok=True)

    for i in range(1, args.n_perm + 1):
        # 1) Shuffle (layer_v, node_v) pairs across rows
        p_edges = permute_target_pairs(edges, rng)

        # 2) Write permuted network to temp file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".txt", dir=perm_dir, delete=False
        ) as tf:
            tf_path = tf.name
        write_multilayer(tf_path, nodes, p_edges)

        # 3) Run Infomap
        n_mod = run_infomap(tf_path, perm_dir, args.infomap)

        # 4) Clean up Infomap output files to save disk space
        base = os.path.splitext(tf_path)[0]
        for ext in (".tree", ".ftree", ".clu", ".map"):
            f = base + ext
            if os.path.exists(f):
                os.remove(f)
        if os.path.exists(tf_path):
            os.remove(tf_path)

        if n_mod is not None:
            perm_counts.append(n_mod)

        if i % 50 == 0 or i == 1:
            mean_str = f"{np.mean(perm_counts):.2f}" if perm_counts else "n/a"
            print(f"  {i}/{args.n_perm} done  "
                  f"(n_valid={len(perm_counts)}, mean={mean_str})")

    perm_counts = np.array(perm_counts, dtype=float)

    # -- Statistics ------------------------------------------------------------
    print("\n=== Statistics ===")
    mu  = np.mean(perm_counts)
    sd  = np.std(perm_counts, ddof=1) if len(perm_counts) > 1 else 0.0
    z   = (n_obs - mu) / sd if sd > 0 else np.nan
    n   = len(perm_counts)

    # +1 correction so p is never exactly 0
    p_right = (np.sum(perm_counts >= n_obs) + 1) / (n + 1)
    p_left  = (np.sum(perm_counts <= n_obs) + 1) / (n + 1)
    p_two   = min(1.0, 2 * min(p_right, p_left))

    print(f"  Observed modules     : {n_obs}")
    print(f"  Null mean +/- SD     : {mu:.2f} +/- {sd:.2f}")
    print(f"  Z-score              : {z:.4f}")
    print(f"  p-value (two-tailed) : {p_two:.4f}")
    print(f"  p-value (right)      : {p_right:.4f}  (obs >= null)")
    print(f"  p-value (left)       : {p_left:.4f}  (obs <= null)")
    print(f"  N successful perms   : {n}")

    # -- Save outputs ----------------------------------------------------------
    np_out = os.path.join(args.outdir, "perm_module_counts.npy")
    np.save(np_out, perm_counts)

    summary = os.path.join(args.outdir, "summary.txt")
    with open(summary, "w") as fh:
        fh.write(f"observed_modules\t{n_obs}\n")
        fh.write(f"null_mean\t{mu:.4f}\n")
        fh.write(f"null_sd\t{sd:.4f}\n")
        fh.write(f"z_score\t{z:.4f}\n")
        fh.write(f"p_two_tailed\t{p_two:.4f}\n")
        fh.write(f"p_right\t{p_right:.4f}\n")
        fh.write(f"p_left\t{p_left:.4f}\n")
        fh.write(f"n_permutations\t{n}\n")

    print(f"\n  Summary    -> {summary}")
    print(f"  Raw counts -> {np_out}")


if __name__ == "__main__":
    main()
