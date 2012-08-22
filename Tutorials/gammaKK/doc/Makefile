
PROJECT = gammaKKdoc

TEX_FILES = $(shell ls *.tex *.bib)

.PHONY all: $(PROJECT)

$(PROJECT): $(TEX_FILES)
#	cp log log.old
	@rm -f log
	pdflatex $@.tex | tee log
	bibtex   $@ | tee log
	pdflatex $@.tex | tee log
	pdflatex $@.tex | tee log
#	dvipdfm $@.dvi
	open $@.pdf

clean:
	@rm -f *.aux
	@rm -f *.log
	@rm -f *.lof
	@rm -f *.lot
	@rm -f *.toc
	@rm -f *.bbl
	@rm -f *.blg
	@rm -f gammaKKdoc.pdf
