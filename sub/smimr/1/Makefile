coxcv.pdf: smimr.tex ../../../main.tex
	cp ../../../main.tex ./
	cp ../../../manuscript_figure/*.png ./
	Rscript scrub-url.R
	sed -i 's/\\citet{Verweij1993}/Verweij and Van Houwelingen\\citep{Verweij1993}/' main.tex
	sed -i 's/\\citet{shedden2008gene}/Shedden et al.\\citep{shedden2008gene}/' main.tex
	sed -i 's/15cm/14cm/' main.tex
	singletex smimr.tex coxcv.tex
	cleantex -beq coxcv
	grep --color includegraphics coxcv.tex

supp.pdf: ../../supplementary.tex
	cd ../../../;\
	singletex supplementary.tex sub/smimr/supp.tex
	sed -i 's/\\bibliographystyle{ims-nourl}/\\bibliographystyle{sageV}/' supp.tex
	cleantex -beq supp
	grep --color includegraphics supp.tex

# all: cover.pdf coxcv.pdf

# cover.pdf: cover.tex
# 	cleantex cover
