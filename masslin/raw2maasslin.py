"""
Convert `taxonomy_Species_abund.txt` and `id_sample.xlsx` into two TSVs
usable by MaAsLin2 (features and metadata).

Produces:
- dataset/maaslin_features.tsv  (samples x features)
- dataset/maaslin_metadata.tsv  (samples x covariates)

Usage:
	python dataset2maasslin.py --tax dataset/taxonomy_Species_abund.txt \
		--meta dataset/id_sample.xlsx --outdir dataset

The script will try to auto-detect the sample ID column in the metadata.
"""

from pathlib import Path
import argparse
import sys

import pandas as pd


def _clean_names(cols):
	"""Sanitize column names: remove tabs/newlines and make unique."""
	out = []
	seen = {}
	for c in cols:
		s = str(c)
		s = s.replace('\t', ' ').replace('\n', ' ').replace('\r', ' ').strip()
		# replace characters that R may treat specially (e.g. '#' is comment.char)
		s = s.replace('#', '_')
		if s == '':
			s = 'NA'
		if s in seen:
			seen[s] += 1
			s = f"{s}__{seen[s]}"
		else:
			seen[s] = 0
		out.append(s)
	return out


def detect_sample_column(meta: pd.DataFrame, sample_ids):
	# Common candidate names
	candidates = [c for c in meta.columns if c.lower() in ("sample","sample_id","sampleid","id","ids","ID","Name")]
	# add the first column as candidate as well
	if len(meta.columns) > 0:
		candidates += [meta.columns[0]]

	# remove duplicates while preserving order
	seen = set()
	candidates = [x for x in candidates if not (x in seen or seen.add(x))]

	best_col = None
	best_match = 0
	for col in candidates:
		vals = meta[col].astype(str).values
		match = len(set(vals) & set(sample_ids))
		if match > best_match:
			best_match = match
			best_col = col

	# fallback: try all columns and pick the one with most overlap
	if best_col is None:
		for col in meta.columns:
			vals = meta[col].astype(str).values
			match = len(set(vals) & set(sample_ids))
			if match > best_match:
				best_match = match
				best_col = col

	return best_col, best_match


def main():
	
	tax_path = Path("dataset/taxonomy_Species_abund.txt")
	meta_path = Path("dataset/id_maaslin.xlsx")
	outdir = Path("dataset/maaslin")
	outdir.mkdir(parents=True, exist_ok=True)

	if not tax_path.exists():
		print(f"taxonomy file not found: {tax_path}", file=sys.stderr)
		sys.exit(1)
	if not meta_path.exists():
		print(f"metadata file not found: {meta_path}", file=sys.stderr)
		sys.exit(1)

	# read taxonomy (features x samples)
	df_tax = pd.read_csv(tax_path, sep="\t", header=0, index_col=0)

	# transpose to samples x features for MaAsLin2
	df_feat = df_tax.T

	# sanitize feature (column) names to avoid embedded tabs/newlines or duplicate names
	df_feat.columns = _clean_names(df_feat.columns)

	# read metadata (try all sheets if needed)
	try:
		df_meta = pd.read_excel(meta_path, engine="openpyxl")
	except Exception:
		df_meta = pd.read_excel(meta_path)

	# drop automatically generated unnamed columns and sanitize metadata column names
	# (these can cause read.table() issues in R)
	unnamed = [c for c in df_meta.columns if str(c).lower().startswith('unnamed')]
	if unnamed:
		print('Dropping metadata columns:', unnamed)
		df_meta = df_meta.drop(columns=unnamed)
	# sanitize remaining metadata column names
	df_meta.columns = _clean_names(df_meta.columns)

	# detect sample id column
	sample_ids = list(df_feat.index.astype(str))
	col, matches = detect_sample_column(df_meta, sample_ids)
	if col is None or matches == 0:
		print("Warning: could not detect sample ID column with overlap; using DataFrame index or first column.")
		if df_meta.index.is_unique and len(set(df_meta.index.astype(str)) & set(sample_ids)) > 0:
			df_meta.index = df_meta.index.astype(str)
		else:
			# use first column as IDs
			first_col = df_meta.columns[0]
			df_meta[first_col] = df_meta[first_col].astype(str)
			df_meta = df_meta.set_index(first_col)
	else:
		df_meta[col] = df_meta[col].astype(str)
		df_meta = df_meta.set_index(col)

	# Ensure indices are strings
	df_meta.index = df_meta.index.astype(str)
	df_feat.index = df_feat.index.astype(str)

	# intersect samples
	common = [s for s in df_feat.index if s in df_meta.index]
	if len(common) == 0:
		# try to match by removing possible prefixes/suffixes (simple heuristics)
		stripped_meta = {s: s.split(".")[0].split("_")[0] for s in df_meta.index}
		stripped_feat = {s: s.split(".")[0].split("_")[0] for s in df_feat.index}
		rev_meta = {}
		for k, v in stripped_meta.items():
			rev_meta.setdefault(v, []).append(k)
		common = []
		for f, fv in stripped_feat.items():
			if fv in rev_meta:
				common.append(f)

	if len(common) == 0:
		print("No overlapping sample IDs found between taxonomy and metadata.", file=sys.stderr)
		print("Taxonomy samples (example):", list(df_feat.index[:10]), file=sys.stderr)
		print("Metadata samples (example):", list(df_meta.index[:10]), file=sys.stderr)
		sys.exit(1)

	# filter and reorder
	df_feat_out = df_feat.loc[common].copy()
	df_meta_out = df_meta.loc[common].copy()

	# write outputs (explicit index label and safe formatting)
	feat_out = outdir / "maaslin_features.tsv"
	meta_out = outdir / "maaslin_metadata.tsv"

	# ensure index label for compatibility with R's read.table(row.names=1)
	df_feat_out.index.name = 'SAMPLE_ID'
	df_meta_out.index.name = 'SAMPLE_ID'

	df_feat_out.to_csv(feat_out, sep="\t", index=True, index_label='SAMPLE_ID', na_rep='')
	df_meta_out.to_csv(meta_out, sep="\t", index=True, index_label='SAMPLE_ID', na_rep='')

	print(f"Wrote features: {feat_out} (samples x features: {df_feat_out.shape})")
	print(f"Wrote metadata: {meta_out} (samples x covariates: {df_meta_out.shape})")


if __name__ == "__main__":
	main()

