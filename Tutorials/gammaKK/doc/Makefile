
PROJECT = gammaKKdoc

TEX_FILES = $(shell ls *.tex)

.PHONY all: $(PROJECT)

$(PROJECT): $(TEX_FILES)
	pdflatex $@.tex
#	bibtex   $@
	pdflatex $@.tex
	pdflatex $@.tex
	open $@.pdf

clean:
	@rm -f *.aux
	@rm -f *.log
	@rm -f *.lof
	@rm -f *.lot
	@rm -f *.toc
	@rm -f *.bbl
	@rm -f *.blg
