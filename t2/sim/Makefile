all: ngspice_1.log ngspice_2.log ngspice_3.log
	$(eval NGSFIGS=$(shell grep _FIG $< | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	$(eval NGSTABS=$(shell grep _TAB $< | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $< | grep -v $i_TAB | grep -v $i_END | grep -v '#' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(eval NGSFIGS=$(shell grep _FIG $(word 2,$^) | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	$(eval NGSTABS=$(shell grep _TAB $(word 2,$^) | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $(word 2,$^) | grep -v $i_TAB | grep -v $i_END | sed 's/vx#branch/Ix/g' | grep -v '#' | sed 's/v(5)-v(7)/Vx/g' | sed 's/(Vx)\/Ix/Req/g' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)
	$(eval NGSFIGS=$(shell grep _FIG $(word 3,$^) | sed 's/_FIG//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(eval NGSFIGPS=$(addsuffix .ps, $(NGSFIGS)))
	$(foreach i, $(NGSFIGPS), ps2pdf $i;)
	$(eval NGSTABS=$(shell grep _TAB $(word 3,$^) | sed 's/_TAB//g' | sed ':a;N;$!ba;s/\n/ /g'))
	$(foreach i, $(NGSTABS), sed -n '/^$i_TAB/,/^$i_END/{p;/^$i_END/q}' $(word 3,$^) | grep -v $i_TAB | grep -v $i_END | sed 's/vx#branch/Ix/g' | grep -v '#' | sed 's/v(5)-v(7)/Vx/g' | sed 's/(Vx)\/Ix/Req/g' | sed 's/\=/\&/g' | sed 's/$$/\\\\ \\hline/g' > $i_tab.tex;)

ngspice_1.log: ../sim/t2_1.net
	ngspice -b $< -o $@

ngspice_2.log: ../sim/t2_2.net
	ngspice -b $< -o $@

ngspice_3.log: ../sim/t2_3.net
	ngspice -b $< -o $@

clean:
	@rm -f ngspice_1.log *.tex *.pdf *.ps *~
	@rm -f ngspice_2.log
	@rm -f ngspice_3.log
	@rm -rf *_ng.txt

.PHONY: all clean
