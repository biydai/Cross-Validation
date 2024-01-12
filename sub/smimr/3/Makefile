all: cover.pdf coxcv.pdf supp.pdf

coxcv.pdf: smimr.tex ../../../main.tex
	cp ../../../main.tex ./
	cp ../../../manuscript_figure/*.png ./
	Rscript ../scrub-url.r
	sed -i 's/\\citet{Verweij1993}/Verweij and Van Houwelingen\\citep{Verweij1993}/' main.tex
	sed -i 's/\\citet{shedden2008gene}/Shedden et al.\\citep{shedden2008gene}/' main.tex
	sed -i 's/15cm/14cm/' main.tex
	singletex smimr.tex coxcv.tex
	cleantex -eq coxcv
	grep --color includegraphics coxcv.tex

supp.pdf: ../../../supplementary.tex
	cp ../../../supplementary.tex ./
	singletex supplementary.tex supp.tex
	sed -i 's/\\bibliographystyle{ims-nourl}/\\bibliographystyle{sageV}/' supp.tex
	cleantex -eq supp
	grep --color includegraphics supp.tex

cover.pdf: cover.tex
	cleantex cover
