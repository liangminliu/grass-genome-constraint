#!/usr/bin/env python3
"""Prune a ROADIES Newick to the MSA panel and write a runtime config.

The 97-tip tree may contain three outgroups. Branch lengths are removed because
msa_pipeline uses this as an alignment guide topology. Original files are kept.
Simple unquoted Newick labels and the bundled config layout are supported.
"""
import argparse
import json
import re
from pathlib import Path


def parse_newick(value):
    text = "".join(value.split())
    if any(mark in text for mark in ("'", '"', "[", "]")):
        raise ValueError("quoted labels and Newick comments need manual handling")
    pos = 0

    def label():
        nonlocal pos
        start = pos
        while pos < len(text) and text[pos] not in ",():;":
            pos += 1
        return text[start:pos]

    def suffix():
        nonlocal pos
        label()  # optional internal-node label
        if pos < len(text) and text[pos] == ":":
            pos += 1
            label()  # branch length

    def subtree():
        nonlocal pos
        if pos >= len(text):
            raise ValueError("truncated Newick")
        if text[pos] == "(":
            pos += 1
            children = [subtree()]
            while pos < len(text) and text[pos] == ",":
                pos += 1
                children.append(subtree())
            if pos >= len(text) or text[pos] != ")":
                raise ValueError("unbalanced Newick parentheses")
            pos += 1
            suffix()
            return children
        name = label()
        if not name:
            raise ValueError("empty Newick tip")
        if pos < len(text) and text[pos] == ":":
            pos += 1
            label()
        return name

    result = subtree()
    if pos >= len(text) or text[pos] != ";" or pos != len(text) - 1:
        raise ValueError("expected one complete Newick tree ending in ';'")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--tree", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    lines = args.config.read_text(encoding="utf-8").splitlines()
    ref = next((m.group(1) for line in lines if (m := re.fullmatch(r"refName:\s*(\S+)\s*", line))), None)
    species = []
    in_species = False
    for line in lines:
        if line == "species:":
            in_species = True
        elif in_species and (m := re.fullmatch(r"\s+-\s+(\S+)\s*", line)):
            species.append(m.group(1))
        elif in_species and line and not line[0].isspace() and not line.startswith("#"):
            in_species = False
    if not ref or not species or len(set([ref, *species])) != len(species) + 1:
        p.error("could not read distinct refName/species IDs from config")
    expected = set([ref, *species])
    parsed = parse_newick(args.tree.read_text(encoding="utf-8"))
    found = []

    def prune(node):
        if isinstance(node, str):
            found.append(node)
            return node if node in expected else None
        children = [kept for child in node if (kept := prune(child)) is not None]
        if not children:
            return None
        if len(children) == 1:
            return children[0]
        return children

    kept = prune(parsed)
    if len(found) != len(set(found)) or not expected.issubset(found) or not isinstance(kept, list):
        p.error(f"tree does not contain distinct configured tips; missing={sorted(expected-set(found))}")

    def render(node):
        return node if isinstance(node, str) else "(" + ",".join(map(render, node)) + ")"

    guide = render(kept) + ";"
    count = 0
    for n, line in enumerate(lines):
        if re.fullmatch(r"speciesTree:\s*.*", line):
            lines[n] = "speciesTree: " + json.dumps(guide)
            count += 1
    if count != 1:
        p.error("expected one speciesTree field in config")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote {args.output}: {len(found)} input tips, {len(expected)} guide tips, "
          f"{len(found)-len(expected)} pruned")


if __name__ == "__main__":
    main()
