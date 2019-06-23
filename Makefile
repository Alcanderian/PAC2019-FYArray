all:
	cd src&&make

run:
	cd src&&make
	cp FYArray.exe ~/preliminary/PAC2019/FYArrayLJX.exe
	cd ~/preliminary/PAC2019&&KMP_AFFINITY=compact ./FYArrayLJX.exe

vtune:export DATE=$(shell date '+%Y-%m-%d_%H-%M-%S')
vtune:
	cd src&&make
	cp FYArray.exe ~/preliminary/PAC2019/FYArrayLJX.exe
	cd ~/preliminary/PAC2019&&KMP_AFFINITY=compact amplxe-cl -collect hotspots -r hs_$(DATE) ./FYArrayLJX.exe
	cp -r src ~/preliminary/PAC2019/hs_$(DATE)
	cp -r include ~/preliminary/PAC2019/hs_$(DATE)
	tar czf hs_$(DATE).tar.gz ~/preliminary/PAC2019/hs_$(DATE)/

.PHONY:clean
clean:
	rm FYArray.exe -f
	rm src/*.o -f
	rm src/FYArray.exe -f
