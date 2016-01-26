default: multigrid clean

multigrid: multiGrid.tex

	pdflatex multiGrid.tex
	pdflatex multiGrid.tex


clean: 
	rm -f *.dvi *.aux *.log *.blg *~

realclean: clean
	rm -f multiGrid.pdf multiGrid.ps *.bbl

html:
	latex2html -split 3 -antialias -transparent -image_type gif multiGrid.tex
