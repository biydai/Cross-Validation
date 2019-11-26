coxcv.pdf: coxcv.tex ../../main.tex
	cp ../../main.tex ./
	cp ../../articles.bib ./
	sed -i '/eprint/d' articles.bib
	sed -i 's/\\citet{Verweij1993}/Verweij and Van Houwelingen\\citep{Verweij1993}/' main.tex
	sed -i 's/\\citet{shedden2008gene}/Shedden et al.\\citep{shedden2008gene}/' main.tex
	cleantex -be coxcv

# all: cover.pdf coxcv.pdf

# cover.pdf: cover.tex
# 	cleantex cover
