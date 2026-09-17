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
