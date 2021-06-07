octave.log: t5.m
	octave $< > $@
	sed -i 's/$$/\\\\ \\hline/g' data1_teo.tex
	sed -i 's/$$/\\\\ \\hline/g' data2_teo.tex
	sed -i 's/$$/\\\\ \\hline/g' data3_teo.tex

clean:
	@rm -f octave.log octave-workspace *.eps *~ *_tab.tex octave-workspace  octave-workspace *.tex

.PHONY: clean
