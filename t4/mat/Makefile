octave.log: t4.m
	octave $< > $@
	sed -i 's/$$/\\\\ \\hline/g' OP_tab_teo.tex
	sed -i 's/$$/\\\\ \\hline/g' data1_teo.tex
	sed -i 's/$$/\\\\ \\hline/g' data2_teo.tex
	sed -i 's/$$/\\\\ \\hline/g' freq_teo.tex
	sed -i 's/$$/\\\\ \\hline/g' merit.tex

clean:
	@rm -f octave.log octave-workspace *.eps *~ *_tab.tex octave-workspace  octave-workspace *.tex
.PHONY: clean
