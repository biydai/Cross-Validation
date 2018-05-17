cross-validation.pdf: cross-validation.tex articles.bib main.tex
	pdflatex cross-validation
	bibtex cross-validation
	pdflatex cross-validation
	pdflatex cross-validation
