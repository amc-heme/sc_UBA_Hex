# sc_UBA_Hex: On-the-fly hexbin computation

The package now computes hexagonal binning on-the-fly for any supported single-cell object using the new generic **get_hexbin_info()**. Highlights:

- **Generic dispatch**: `get_hexbin_info(object)` returns a list with `cID` and `hexbin.matrix`.
- **Supported classes**: Methods implemented for `SingleCellExperiment` and `Seurat` objects.
- **No metadata writes**: Internal extractors `.extract_cID` and `.extract_hexbin` compute bins at plotting time, so original objects remain unchanged.
- **SCUBA integration**: Uses SCUBA to read from additional object types.

You can now call any plot function directly on supported objects, for example:

```r
library(Seurat)
seurat_obj <- Load10X_Spatial("path/to/data")
seurat_obj <- RunUMAP(seurat_obj)
plot_hexbin_density(seurat_obj)
```

## Building and testing {#build-test}

To build and check the package locally from the terminal:
```bash
R CMD build .
R CMD check schex_*.tar.gz --as-cran
```

Alternatively, using **devtools** within R:
```r
# install.packages("devtools")
devtools::load_all()       # load all package functions
devtools::document()       # update documentation from roxygen comments
devtools::test()           # run unit tests (testthat)
devtools::check()          # run package checks
```

For Bioconductor compliance, you can also run:
```r
# install.packages("BiocManager")
BiocManager::install("BiocCheck")
BiocCheck::BiocCheck(".")
```

