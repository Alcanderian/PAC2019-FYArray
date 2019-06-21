DATE=`date '+%Y-%m-%d_%H-%M-%S'`

all:
	cd src&&make

run:
	cd src&&make
	cp FYArray.exe ~/preliminary/PAC2019/FYArrayLJX.exe
	cd ~/preliminary/PAC2019&&./FYArrayLJX.exe

vtune:
	cd src&&make
	cp FYArray.exe ~/preliminary/PAC2019/FYArrayLJX.exe
	cd ~/preliminary/PAC2019&&amplxe-cl -collect hotspots -r hs_$(DATE) ./FYArrayLJX.exe

.PHONY:clean
clean:
	rm FYArray.exe -f
	rm src/*.o -f
	rm src/FYArray.exe -f
