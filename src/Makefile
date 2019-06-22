TARGET=FYArray.exe

CXX=icpc
CXXFLAGS=-O2 -xHost -qopt-report-phase=vec -qopt-report=5 -qopt-zmm-usage=high -qopenmp
LDFLAGS=-O2 -xHost -qopt-report-phase=vec -qopt-report=5 -qopt-zmm-usage=high -qopenmp

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
