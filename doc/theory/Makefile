all: paircomp-theory.pdf

clean:
	rm -f *.ps *.pdf *.dvi *.toc *.log *.aux

paircomp-theory.pdf: paircomp-theory.ps
	ps2pdf paircomp-theory.ps

paircomp-theory.ps: paircomp-theory.tex
	latex paircomp-theory.tex && latex paircomp-theory.tex
	dvips -Pcmz -Z9 -o paircomp-theory.ps paircomp-theory.dvi
