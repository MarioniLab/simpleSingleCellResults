rmarkdown::render("work-1-reads.Rmd", output_dir=".", intermediates_dir=".", clean=FALSE)
rmarkdown::render("work-2-umis.Rmd", output_dir=".", intermediates_dir=".", clean=FALSE)
rmarkdown::render("work-3-tenx.Rmd", output_dir=".", intermediates_dir=".", clean=FALSE)
rmarkdown::render("work-5-mnn.Rmd", output_dir=".", intermediates_dir=".", clean=FALSE)

rmarkdown::render("xtra-1-qc.Rmd", output_dir=".", intermediates_dir=".", clean=FALSE)
rmarkdown::render("xtra-2-spike.Rmd", output_dir=".", intermediates_dir=".", clean=FALSE)
rmarkdown::render("xtra-3-var.Rmd", output_dir=".", intermediates_dir=".", clean=FALSE)
rmarkdown::render("xtra-4-misc.Rmd", output_dir=".", intermediates_dir=".", clean=FALSE)

