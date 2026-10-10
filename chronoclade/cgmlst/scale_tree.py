"""Compiled NJ, iterative PCoA, indexed tree distances and bounded HTML views."""

from __future__ import annotations
import io
import json
from pathlib import Path
import subprocess
import numpy as np
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree


def leading_pcoa(matrix, *, seed=42, batch_size=128, tolerance=1e-7):
    """Leading axes of the exact centered squared-distance operator.

    The eigensolver is iterative. No sample or distance is approximated/subsetted;
    the uncomputed negative spectrum and positive inertia fraction are unknown.
    """
    from scipy.sparse.linalg import LinearOperator, eigsh

    n = len(matrix)
    if n <= 3:
        centered = np.eye(n) - np.ones((n, n)) / n
        gram = -0.5 * centered @ (matrix**2) @ centered
        values, vectors = np.linalg.eigh(gram)
        order = np.argsort(values)[::-1][:2]
        values, vectors = values[order], vectors[:, order]
    else:

        def apply(vector):
            vector = vector - vector.mean()
            result = np.empty(n)
            for start in range(0, n, batch_size):
                result[start : start + batch_size] = (
                    np.asarray(matrix[start : start + batch_size]) ** 2
                ) @ vector
            result -= result.mean()
            return -0.5 * result

        operator = LinearOperator((n, n), matvec=apply, dtype=np.float64)
        values, vectors = eigsh(
            operator,
            k=2,
            which="LA",
            tol=tolerance,
            v0=np.random.default_rng(seed).normal(size=n),
            maxiter=500,
        )
        order = np.argsort(values)[::-1]
        values, vectors = values[order], vectors[:, order]
    # Axis signs are arbitrary; fixed largest loading convention makes exports repeatable.
    for column in range(vectors.shape[1]):
        anchor = int(np.argmax(np.abs(vectors[:, column])))
        if vectors[anchor, column] < 0:
            vectors[:, column] *= -1
    coords = np.zeros((n, 2))
    coords[:, : len(values)] = vectors * np.sqrt(np.maximum(values, 0))
    return coords, dict(
        pcoa_eigenvalues=values.tolist(),
        pcoa_positive_axis_fraction=None,
        pcoa_solver="iterative ARPACK leading eigenpairs of the exact all-profile centered distance operator",
        pcoa_tolerance=tolerance,
        pcoa_full_spectrum_computed=False,
        pcoa_inertia_note="Total positive inertia and negative spectrum are not calculated by the leading-axis solver.",
    )


class IndexedTree(Tree):
    """Exact LCA branch distances for network roots without repeated tree searches."""

    @classmethod
    def from_tree(cls, tree):
        indexed = cls(root=tree.root, rooted=tree.rooted, id=tree.id, name=tree.name)
        indexed.reindex()
        return indexed

    def reindex(self):
        self._nodes, parents, depths, levels = [], [], [], []
        stack = [(self.root, 0, 0.0, 0)]
        while stack:
            node, parent, depth, level = stack.pop()
            index = len(self._nodes)
            self._nodes.append(node)
            parents.append(parent)
            depths.append(depth)
            levels.append(level)
            for child in reversed(node.clades):
                stack.append((child, index, depth + (child.branch_length or 0.0), level + 1))
        self._node_index = {node: i for i, node in enumerate(self._nodes)}
        self._name_index = {node.name: i for i, node in enumerate(self._nodes) if node.name}
        self._depths = depths
        self._levels = levels
        self._up = [np.asarray(parents, dtype=np.int64)]
        for _ in range(max(levels, default=0).bit_length()):
            self._up.append(self._up[-1][self._up[-1]])

    def distance(self, target1, target2=None):
        def index(target):
            return self._name_index[target] if isinstance(target, str) else self._node_index[target]

        a, b = index(target1), 0 if target2 is None else index(target2)
        start_a, start_b = a, b
        if self._levels[a] < self._levels[b]:
            a, b = b, a
        difference = self._levels[a] - self._levels[b]
        for level in range(difference.bit_length()):
            if difference & (1 << level):
                a = int(self._up[level][a])
        if a != b:
            for ancestors in reversed(self._up):
                if ancestors[a] != ancestors[b]:
                    a, b = int(ancestors[a]), int(ancestors[b])
            a = int(self._up[0][a])
        return self._depths[start_a] + self._depths[start_b] - 2 * self._depths[a]

    def root_with_outgroup(self, *args, **kwargs):
        super().root_with_outgroup(*args, **kwargs)
        self.reindex()

    def prune(self, *args, **kwargs):
        result = super().prune(*args, **kwargs)
        self.reindex()
        return result


def rapidnj_tree(matrix, sample_ids, output, *, executable=None, memory_mb=2048):
    """Stream exact PHYLIP to real RapidNJ; verify restored full tip identity."""
    from .scale_setup import rapidnj_executable

    binary = rapidnj_executable(executable)
    output = Path(output)
    names = [f"CC{i:09d}" for i in range(len(sample_ids))]
    source = output.with_suffix(".phy")
    with source.open("w") as handle:
        handle.write(str(len(names)) + "\n")
        for i, name in enumerate(names):
            handle.write(
                name + " " + " ".join(format(float(value), ".17g") for value in matrix[i]) + "\n"
            )
    command = [binary, str(source), "-i", "pd", "-o", "t", "-m", str(memory_mb)]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    log = output.with_suffix(".rapidnj.log")
    log.write_text(completed.stderr)
    if completed.returncode:
        raise RuntimeError(f"RapidNJ failed ({completed.returncode}); see {log}")
    tree = Phylo.read(io.StringIO(completed.stdout), "newick")
    mapping = dict(zip(names, sample_ids))
    tips = tree.get_terminals()
    if len(tips) != len(names) or {tip.name for tip in tips} != set(names):
        raise RuntimeError("RapidNJ output does not contain every input sample exactly once")
    for tip in tips:
        tip.name = mapping[tip.name]
    tree = IndexedTree.from_tree(tree)
    tree.root_with_outgroup(sample_ids[0])
    Phylo.write(tree, output, "newick")
    source.unlink()  # Interoperability scratch, not quadratic CSV evidence.
    return tree, dict(
        backend="RapidNJ",
        executable=binary,
        command=command,
        executable_sha256=__import__("chronoclade.artifacts", fromlist=["file_sha256"]).file_sha256(
            binary
        ),
        memory_budget_mb=memory_mb,
        stderr_path=str(log),
        sample_count=len(sample_ids),
        distance_precision="PHYLIP 17 significant digits; RapidNJ float32 internals",
    )


def write_collapsible_tree(network, rows, labels, path):
    """Offline lazy tree expansion; full state audit remains the source of truth."""
    path = Path(path)
    audit = network["reconstruction"]
    by_id = {row["sample_id"]: row for row in rows}
    data = []
    for node in audit["nodes"]:
        record = by_id.get(node["name"], {})
        ident = record.get("sample_id")
        data.append(
            dict(
                id=node["id"],
                children=node["children_ids"],
                name=node["name"],
                label=labels.get(ident, ident) if ident else None,
                state=node.get("state"),
                country=record.get("country") or "Country unknown",
                date=str(
                    record.get("collection_date") or record.get("collection_year") or "Date unknown"
                ),
                role=record.get("role") or record.get("origin"),
                branch=node.get("branch_length"),
            )
        )
    payload = json.dumps(
        dict(nodes=data, root=audit["root_id"], colors=network["country_colors"]),
        separators=(",", ":"),
    ).replace("<", "\\u003c")
    path.write_text(
        """<!doctype html><html lang="en"><meta charset="utf-8"><title>Full cgMLST tree</title>
<style>body{font:14px system-ui;margin:24px;color:#27323b}button{background:#fff;border:1px solid #ccd4da;border-radius:5px;padding:6px 10px;margin:4px;cursor:pointer}.branch{border-left:2px solid #ced5dc;margin-left:18px;padding:4px 0 4px 10px}.tip{padding:6px;border-left:4px solid #aaa;margin:4px 0}.input{font-weight:700}.muted{color:#667}</style>
<h1>Complete country-coloured cgMLST NJ tree</h1><p>Expand clades to inspect every genome. All tip IDs and identical-profile multiplicities are retained. Branch values are mismatch fractions; negative branches are retained. Colours show the same representative ancestral country states as the network.</p>
<p id="count"></p><button id="reset">Collapse all</button><div id="tree"></div><script>
const data="""
        + payload
        + """;const nodes=new Map(data.nodes.map(n=>[n.id,n]));const totals=new Map();
for(let i=data.nodes.length-1;i>=0;i--){let n=data.nodes[i];totals.set(n.id,n.children.length?n.children.reduce((s,c)=>s+totals.get(c),0):1)}
document.getElementById('count').textContent=totals.get(data.root)+' genomes in the complete tree';
function render(id){const n=nodes.get(id),box=document.createElement('div');box.className='branch';box.style.borderColor=data.colors[n.state]||'#aaa';
if(!n.children.length){box.className='tip'+(['input','local','query','focal'].includes(n.role)?' input':'');box.style.borderColor=data.colors[n.country]||'#aaa';box.textContent=n.label+' | '+n.country+' | '+n.date+' | ID '+n.name+' | branch '+n.branch;return box}
let button=document.createElement('button'),children=document.createElement('div');button.textContent='▸ '+totals.get(id)+' genomes · '+(n.state||'Unknown country')+' · branch '+n.branch;let expanded=false;button.onclick=()=>{expanded=!expanded;children.replaceChildren();if(expanded)n.children.forEach(c=>children.append(render(c)));button.textContent=(expanded?'▾ ':'▸ ')+totals.get(id)+' genomes · '+(n.state||'Unknown country')+' · branch '+n.branch};box.append(button,children);return box}
function reset(){document.getElementById('tree').replaceChildren(render(data.root))}document.getElementById('reset').onclick=reset;reset();
</script></html>"""
    )
    return str(path)
