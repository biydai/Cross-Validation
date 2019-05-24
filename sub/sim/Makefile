all: cover.pdf coxcv.pdf

coxcv.pdf: coxcv.tex ../../main.tex
	cp ../../main.tex ./
	sed -i 's/\\citet{Verweij1993}/Verweij and Van Houwelingen\\citep{Verweij1993}/' main.tex
	sed -i 's/\\citet{shedden2008gene}/Shedden et al.\\citep{shedden2008gene}/' main.tex
	cleantex -be coxcv
	rm coxcv.pag

cover.pdf: cover.tex
	cleantex cover
