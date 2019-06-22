TARGET=FeatureTest.exe

CXX=icpc
CXXFLAGS=-O2 -xHost -qopt-report-phase=vec -qopt-report=5 -qopt-zmm-usage=high -g

INCLUDE=../include

OBJ=main.o

$(TARGET):$(OBJ)
	$(CXX) $^ -o $@ 
	cp $(TARGET) ../

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) -c $< -o $@ 

.PHONY:clean
clean:
	rm *.o -f
