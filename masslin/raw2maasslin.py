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


def process_feature_file(tax_path, df_meta, out_path, sep="\t"):
	if not tax_path.exists():
		print(f"File not found: {tax_path}", file=sys.stderr)
		return

	# read features (features x samples or samples x features)
	# Taxonomy usually features x samples
	df_tax = pd.read_csv(tax_path, sep=sep, header=0, index_col=0)

	# If the first metadata sample IDs match columns, it's features x samples
	# otherwise if they match index, it's samples x features.
	# But common case here is features x samples needs transpose.
	# Check if index contains samples
	sample_ids = list(df_meta.index.astype(str))
	overlap_index = len(set(df_tax.index.astype(str)) & set(sample_ids))
	overlap_cols = len(set(df_tax.columns.astype(str)) & set(sample_ids))

	if overlap_cols > overlap_index:
		df_feat = df_tax.T
	else:
		df_feat = df_tax

	# sanitize feature (column) names
	df_feat.columns = _clean_names(df_feat.columns)
	df_feat.index = df_feat.index.astype(str)

	# intersect samples
	common = [s for s in df_feat.index if s in df_meta.index]
	if len(common) == 0:
		# try to match by removing possible prefixes/suffixes
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
		print(f"No overlapping sample IDs found for {tax_path}", file=sys.stderr)
		return

	# filter and reorder
	df_feat_out = df_feat.loc[common].copy()
	df_feat_out.index.name = 'SAMPLE_ID'
	df_feat_out.to_csv(out_path, sep="\t", index=True, index_label='SAMPLE_ID', na_rep='')
	print(f"Wrote features to {out_path} (shape: {df_feat_out.shape})")
	return common


def main():
	
	tax_path = Path("dataset/taxonomy_Species_abund.txt")
	metabolite_path = Path("dataset/metabolome_all_data.csv")
	ko_path = Path("dataset/KEGG/4_KOEntry/KEGG_KOEntry_abund.txt")
	meta_path = Path("dataset/id_sample.xlsx")
	outdir = Path("dataset/maaslin")
	outdir.mkdir(parents=True, exist_ok=True)

	if not meta_path.exists():
		print(f"metadata file not found: {meta_path}", file=sys.stderr)
		sys.exit(1)

	# read metadata
	try:
		df_meta = pd.read_excel(meta_path, engine="openpyxl")
	except Exception:
		df_meta = pd.read_excel(meta_path)

	unnamed = [c for c in df_meta.columns if str(c).lower().startswith('unnamed')]
	if unnamed:
		df_meta = df_meta.drop(columns=unnamed)
	df_meta.columns = _clean_names(df_meta.columns)

	# detect sample id column - we need some sample IDs to match against
	# Just read enough of one file to get sample IDs if possible, or use heuristic
	# For now, initialize df_meta with index if we can't detect.
	# Actually, let's just pick any feature file to get sample IDs for detection
	df_tax_tmp = pd.read_csv(tax_path, sep="\t", header=0, index_col=0, nrows=1)
	sample_ids = list(df_tax_tmp.columns.astype(str))
	
	col, matches = detect_sample_column(df_meta, sample_ids)
	if col is None or matches == 0:
		if df_meta.index.is_unique:
			df_meta.index = df_meta.index.astype(str)
		else:
			first_col = df_meta.columns[0]
			df_meta = df_meta.set_index(first_col)
	else:
		df_meta = df_meta.set_index(col)
	df_meta.index = df_meta.index.astype(str)
	df_meta.index.name = 'SAMPLE_ID'

	# Process Taxonomy
	common_samples = process_feature_file(tax_path, df_meta, outdir / "maaslin_raw_features.tsv", sep="\t")

	# Process Metabolome
	process_feature_file(metabolite_path, df_meta, outdir / "maaslin_raw_metabolite_features.tsv", sep=",")

	# Process KO
	process_feature_file(ko_path, df_meta, outdir / "maaslin_raw_ko_features.tsv", sep="\t")

	# Write metadata (only for common samples of taxonomy by default, or all if preferred)
	# Usually MaAsLin filters metadata to match features anyway.
	if common_samples:
		df_meta_out = df_meta.loc[common_samples].copy()
		meta_out = outdir / "maaslin_raw_metadata.tsv"
		df_meta_out.to_csv(meta_out, sep="\t", index=True, index_label='SAMPLE_ID', na_rep='')
		print(f"Wrote metadata: {meta_out} (shape: {df_meta_out.shape})")


if __name__ == "__main__":
	main()

