octave.log: t3.m
	octave $< > $@
	sed -i 's/$$/\\\\ \\hline/g' Envelope_tab.tex
	sed -i 's/$$/\\\\ \\hline/g' Regulator_tab.tex

clean:
	@rm -f octave.log octave-workspace *.eps *~ *_tab.tex octave-workspace  octave-workspace 
	@rm -rf data.txt
.PHONY: clean
