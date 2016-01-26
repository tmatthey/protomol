default: userguide	

userguide: userguide.tex userguide.bbl
	pdflatex userguide
	pdflatex userguide

clean: 
	rm -f *.dvi *.aux *.log *.blg *.toc *.lof *~

realclean: clean
	rm -f userguide.pdf userguide.ps *.bbl

userguide.bbl:
	pdflatex userguide
	bibtex userguide
html:
	latex2html -split 3 -antialias -transparent -image_type gif userguide.tex
