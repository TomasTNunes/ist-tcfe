octave.log: t2.m
	octave $< > $@
	mv *_ng.txt ../sim
	sed -i 's/$$/\\\\ \\hline/g' Nodal1_tab.tex
	sed -i 's/tau/\\tau/g' Req_tau_tab.tex
	sed -i 's/$$/\\\\ \\hline/g' Req_tau_tab.tex
	sed -i 's/$$/\\\\ \\hline/g' phasors_tab.tex

clean:
	@rm -f octave.log octave-workspace *.eps *~ *_tab.tex octave-workspace  octave-workspace 
	@rm -rf data.txt
.PHONY: clean
