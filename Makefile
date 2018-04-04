cross-validation.pdf: cross-validation.tex articles.bib
	pdflatex cross-validation
	bibtex cross-validation
	pdflatex cross-validation
	pdflatex cross-validation
