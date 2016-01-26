default: hgroup clean

hgroup: hgroup.tex hgroup.bbl
	pdflatex hgroup.tex
	pdflatex hgroup.tex

hgroup.bbl:
	pdflatex hgroup
	bibtex hgroup

clean: 
	rm -f *.dvi *.aux *.log *.blg *~

realclean: clean
	rm -f hgroup.pdf hgroup.ps *.bbl
