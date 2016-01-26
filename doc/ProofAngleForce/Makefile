default: proof clean

proof: proof.tex

	pdflatex proof.tex
	pdflatex proof.tex


clean: 
	rm -f *.dvi *.aux *.log *.blg *~

realclean: clean
	rm -f proof.pdf proof.ps *.bbl
