cross-validation.pdf: cross-validation.tex articles.bib main.tex
	cleantex -be cross-validation

supplementary.pdf: supplementary.tex
	cleantex -be supplementary
