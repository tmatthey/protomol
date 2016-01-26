default: ewald	

ewald: ewald.tex ewald.bbl
	pdflatex ewald
	pdflatex ewald

clean: 
	rm -f *.dvi *.aux *.log *.blg *.toc *.lof *~

realclean: clean
	rm -f ewald.pdf ewald.ps *.bbl

ewald.bbl:
	pdflatex ewald
	bibtex ewald
html:
	latex2html -split 3 -antialias -transparent -image_type gif ewald.tex
