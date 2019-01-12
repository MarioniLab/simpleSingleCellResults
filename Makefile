targets: all

RCMD=R

# Obtaining the vignette sources.
src:
	 git clone https://github.com/MarioniLab/simpleSingleCell src

update: src
	cd src && git pull
	${RCMD} CMD INSTALL --preclean src/

# Moving the scripts to avoid directory chaos.
%.Rmd: src
	cp src/vignettes/$@ .

ref.bib: src
	cp src/vignettes/ref.bib .

# Defining the HTML outputs.
%.knit.md : %.Rmd ref.bib
	${RCMD} --no-save --slave -e "rmarkdown::render('$<', clean=FALSE)"

all: batch.knit.md bigdata.knit.md de.knit.md doublets.knit.md intro.knit.md misc.knit.md qc.knit.md reads.knit.md spike.knit.md tenx.knit.md umis.knit.md var.knit.md

# Cleaning commands.
uncache:
	rm -rf *_cache

clean:
	rm -rf *_cache *.html *_files *.Rmd raw_data
	git checkout *.knit.md
	cd src && git reset --hard HEAD
