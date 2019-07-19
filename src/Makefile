TARGET=FYArray.exe

CXX=icpc
CXXFLAGS=-O3 -xHost -qopt-report-phase=vec -qopt-report=5 -qopt-zmm-usage=high -g -qopenmp -ltbbmalloc_proxy
LDFLAGS=-qopenmp -ltbbmalloc_proxy

INCLUDE=../include

OBJ=FYStorage.o main.o

$(TARGET):$(OBJ)
	$(CXX) $(LDFLAGS) $^ -o $@ 
	cp $(TARGET) ../

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) -c $< -o $@ 

.PHONY:clean
clean:
	rm *.o -f
