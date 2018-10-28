# simpleSingleCell results

This repository contains compiled Markdown files for the [_simpleSingleCell_](http://bioconductor.org/packages/devel/workflows/html/simpleSingleCell.html) workflow.
Each file represents the compiled version of one of the vignettes.
The aim is to allow the maintainers to easily diagnose changes to the workflow results upon updates to dependent packages.
To check the results on a fresh system:

1. Run `setup.sh` to pull the [`simpleSingleCell`](https://github.com/MarioniLab/simpleSingleCell) repository adjacent to this one.
This assumes you do not have a `package/` directory next to the current directory.
2. Run `runner.R` to compile them while retaining the intermediate Markdown files.
This requires all packages dependent on _simpleSingleCell_ to be installed:

    ```
    BiocManager::install("simpleSingleCell")
    ```

3. A simple `git diff` will then highlight the changes to the compiled results.
Images are not diff'd but most code chunks will report numbers for checking.

Please report any discrepancies as issues at the [`simpleSingleCell`](https://github.com/MarioniLab/simpleSingleCell) repository, rather than this repository.
