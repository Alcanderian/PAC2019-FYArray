CHECK_TXT_PATH=~/preliminary/PAC2019

all:
	cd src&&make

run:
	cd src&&make
	mv FYArray.exe $(CHECK_TXT_PATH)/FYArray.exe
	cd $(CHECK_TXT_PATH)&&KMP_AFFINITY=compact OMP_NUM_THREADS=34 ./FYArray.exe

vtune:export DATE=$(shell date '+%Y-%m-%d_%H-%M-%S')
vtune:
	cd src&&make
	mv FYArray.exe $(CHECK_TXT_PATH)/FYArray.exe
	cd $(CHECK_TXT_PATH)&&KMP_AFFINITY=compact OMP_NUM_THREADS=34 amplxe-cl -collect hotspots -r hs_$(DATE) ./FYArray.exe
	cp -r src $(CHECK_TXT_PATH)/hs_$(DATE)
	cp -r include $(CHECK_TXT_PATH)/hs_$(DATE)
	tar czf hs_$(DATE).tar.gz $(CHECK_TXT_PATH)/hs_$(DATE)/

.PHONY:clean
clean:
	rm FYArray.exe -f
	rm src/*.o -f
	rm src/FYArray.exe -f
