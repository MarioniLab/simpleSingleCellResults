targets: all

RCMD=R

# Obtaining the vignette sources.
src:
	 git clone https://github.com/MarioniLab/simpleSingleCell src

update: src
	cd src && git pull
	${RCMD} CMD INSTALL --preclean src/

# Moving the scripts to avoid directory chaos.
RMD=batch.Rmd bigdata.Rmd de.Rmd doublets.Rmd intro.Rmd misc.Rmd qc.Rmd reads.Rmd spike.Rmd tenx.Rmd umis.Rmd var.Rmd
$(RMD): src
	cp src/vignettes/*.Rmd .

# Defining the HTML outputs.
HTML=$(RMD:.Rmd=.html)
$(HTML): %.html: %.Rmd
	rm -f $@
	${RCMD} --no-save --slave -e "rmarkdown::render('$<', clean=FALSE)"

all: $(HTML)

# Cleaning commands.
uncache:
	rm -rf *_cache

clean:
	rm -rf *_cache *.html *_files *.Rmd raw_data
	git checkout *.knit.md
	cd src && git reset --hard HEAD
