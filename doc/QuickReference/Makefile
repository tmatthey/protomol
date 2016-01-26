default: quickref	

quickref: quickref.tex quickref.bbl
	pdflatex quickref
	pdflatex quickref

clean: 
	rm -f *.dvi *.aux *.log *.blg *.toc *.lof *~

realclean: clean
	rm -f quickref.pdf quickref.ps *.bbl

quickref.bbl:
	pdflatex quickref
	bibtex quickref
html:
	latex2html -split 3 -antialias -transparent -image_type gif quickref.tex
