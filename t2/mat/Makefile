octave.log: t2.m
	octave $< > $@
	mv data_ng.txt ../sim
	sed -i 's/$$/\\\\ \\hline/g' Nodal1_tab.tex

clean:
	@rm -f octave.log octave-workspace *.eps *~ *_tab.tex octave-workspace 
	@rm -rf data.txt
.PHONY: clean
