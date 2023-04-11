FILE = YarnCurve

all:
	pdflatex -shell-escape $(FILE)
	bibtex $(FILE)
	pdflatex -shell-escape $(FILE)

clean:
	rm -f YarnCurve.aux
	rm -f YarnCurve.bbl
	rm -f YarnCurve.blg
	rm -f YarnCurve.log
	rm -f YarnCurve.pdf
	rm -rf _minted-YarnCurve
