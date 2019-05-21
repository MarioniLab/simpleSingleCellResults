targets: all

RCMD=Rdevel

# Obtaining the vignette sources.
src:
	 git clone https://github.com/MarioniLab/simpleSingleCell src

update: src
	cd src && git pull
	${RCMD} CMD INSTALL --preclean src/

# Moving the scripts to avoid directory chaos.
src/vignettes/%.Rmd: src

%.Rmd: src/vignettes/%.Rmd
	cp $< $@

ref.bib: src
	cp src/vignettes/ref.bib .

.PRECIOUS: %.Rmd # do not delete upon completion!

# Creating the *.knit.md file (destroying it if the command fails).
%.knit.md : %.Rmd ref.bib
	${RCMD} --no-save --slave -e "rmarkdown::render('$<', clean=FALSE)" || rm $@

bigdata.knit.md: tenx.knit.md

de.knit.md: reads.knit.md batch.knit.md

misc.knit.md: var.knit.md

qc.knit.md: reads.knit.md tenx.knit.md

var.knit.md: reads.knit.md

all: intro.knit.md \
        reads.knit.md \
        umis.knit.md \
        tenx.knit.md \
        batch.knit.md \
        doublets.knit.md \
        qc.knit.md \
        spike.knit.md \
        var.knit.md \
        de.knit.md \
        bigdata.knit.md \
        multibatch.knit.md \
        misc.knit.md

# Cleaning commands.
clean: 
	rm -rf *.Rmd *.html *_files *_cache *.utf8.md *.rds *.tsv *.xls
	rm -rf *.knit.md

distclean: clean
	rm -rf raw_data
	rm -rf src
