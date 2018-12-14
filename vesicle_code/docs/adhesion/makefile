filename=adhesion
compileLatexOptions=--quiet --c-style-errors
compilePdfOptions=

compile= ${compPdftex}
compLatex=latex ${compileLatexOptions}
compPdftex=pdflatex ${compilePdfOptions}

all:
	${compile}  ${filename}
	-bibtex ${filename}
	bibtex ${filename}
	${compile} ${filename}
	${compile} ${filename}
	@echo Done ....!

dvi:
	latex ${compileOptions} ${filename}

bib:
	${compile} ${filename}
	bibtex ${filename}

fast:
	${compile}  ${filename}

clean:
	rm *.aux *.log *.spl *.blg *.bbl 
	@echo Cleaned ....! 

%.pdf: %.tex
	${compile} $*
	${compile} $* 

