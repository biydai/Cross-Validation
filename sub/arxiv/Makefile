coxcv.pdf: ../../main.tex
	cd ../../;\
	singletex cross-validation.tex > sub/arxiv/coxcv.tex ;\
	cp manuscript_figure/*.png sub/arxiv/
	cleantex -be coxcv
	grep --color includegraphics coxcv.tex
