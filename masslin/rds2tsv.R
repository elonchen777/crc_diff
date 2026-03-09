fit_outdir <- file.path(outdir, "fit_exports")
if (!dir.exists(fit_outdir)) dir.create(fit_outdir, recursive = TRUE)

# save full R object
# saveRDS(fit, file = file.path(fit_outdir, "maaslin_fit.rds"))
fit <- readRDS(file.path(fit_outdir, "maaslin_fit.rds"))

# iterate over top-level components and write any data.frames or matrices
for (nm in names(fit)) {
   message("Processing fit component: ", nm)
	obj <- fit[[nm]]
	if (is.data.frame(obj) || is.matrix(obj)) {
		fn <- file.path(fit_outdir, paste0(nm, ".tsv"))
		write.table(as.data.frame(obj), file = fn, sep = "\t", row.names = FALSE, quote = FALSE)
		message("Wrote: ", fn)
	} else if (is.list(obj)) {
		# if list contains data.frames, write them too
		for (subnm in names(obj)) {
			subobj <- obj[[subnm]]
			if (is.data.frame(subobj) || is.matrix(subobj)) {
				fn <- file.path(fit_outdir, paste0(nm, "_", subnm, ".tsv"))
				write.table(as.data.frame(subobj), file = fn, sep = "\t", row.names = FALSE, quote = FALSE)
				message("Wrote: ", fn)
			}
		}
	}
}