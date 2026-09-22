# PDACMOC 2.6.1

Documentation only. The package and the classification results are identical to 2.6.0.

* The README metric badges are grouped by source (GitHub, citations, Shiny app, Zenodo)
  and include the citations of the bioRxiv preprint.
* The example script (`inst/examples/example.R`), the README and the help pages show the
  same example again, pointing `reticulate` to the Python of the `pdacmoc` environment.

# PDACMOC 2.6.0

New Shiny app and simpler installation. Classification results are identical to 2.5.5.

* Redesigned Shiny app with five tabs (Classify, Results, Summary, Performance, Help),
  subtypes in the colours used in the paper, a High/Low confidence column, subtype
  summaries and the published balanced accuracy of every classifier. Light and dark mode.
* Classifications run in the background: the app stays responsive, shows the progress and
  the estimated time, and a classification can be cancelled at any time. At most two
  classifications run at the same time on a server (`options(PDACMOC.max_jobs = n)`);
  later ones wait for a free slot.
* Messages about the file and the classification (e.g. averaged or imputed genes) are
  listed in the Run card instead of pop-up notifications.
* 'Reset app' is now 'New classification', available once there are results.
* The app accepts .tsv, .csv and .txt files. `read.counts()` is exported to read files the
  same way from R (the former internal `read.expression.file()` still works).
* Installation: `install/install.sh` creates a conda environment with all the dependencies,
  including the pinned scikit-learn 1.3.1 and org.Hs.eg.db 3.18.0. The README explains it
  in three commands.
* The public server counts finished classifications and classified samples (totals only,
  nothing about the users or their data) for the usage badge of the README. The counter is
  off unless `options(PDACMOC.stats_dir)` is set.
* New dependencies: bslib and ggplot2.

# PDACMOC 2.5.5

Bug fixes. Classification results for inputs that already worked in 2.5.4 are unchanged.

* The Shiny app no longer freezes when an error happens: errors are shown as notifications
  and the app can be reset.
* Gene IDs that map to the same Ensembl ID (gene symbols or Entrez IDs) and repeated gene IDs
  in the input file are now averaged instead of stopping the classification.
* The Entrez ID option works (it used an invalid key type).
* Classification of a single sample no longer fails. Single-sample input is neither
  batch-corrected nor variance-stabilized; this will be revised in a future version.
* Non-integer counts (e.g. after averaging duplicated genes without batch correction) are
  rounded before the variance stabilizing transformation instead of failing.
* Clear error messages for non-numeric, negative or missing values, unrecognised gene IDs and
  missing training data.
* The Shiny app shows the package messages (imputed genes, averaged duplicates) as notifications.
* The Shiny app runs the classification only once per click, even if several files were uploaded.
* Downloading the stroma PDAConsensus classification returned an empty file.
* Virtual microdissection writes its temporary plots to a unique temporary file.
* Virtual microdissection is faster: it uses 4 parallel workers by default (set
  `options(PDACMOC.cores = n)` to change it). The bundled ADVOCATE package is updated to
  0.1.0.1, whose only change is reading the number of workers from the `ADVOCATE.cores` option.
* Virtual microdissection is retried once if a parallel worker fails.
