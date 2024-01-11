all: cross-validation.pdf supplementary.pdf

cross-validation.pdf: cross-validation.tex articles.bib main.tex
	cleantex --bbl -btq cross-validation

supplementary.pdf: supplementary.tex
	cleantex -beq supplementary
