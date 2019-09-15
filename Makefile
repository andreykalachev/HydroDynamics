TARGET = HydroDynamicsRun

BUILDDIR = build
SOURCEDIR = HydroDynamics
RESULTS_OTPUT = Results.txt
LIBSDIR = fade3d/lib_ubuntu17.04_x86_64
LIBS = fade3d

SRC = HydroDynamics.cpp
OBJECTS = HydroDynamics.o

CXX = g++

.PHONY: all, clean

all: $(TARGET)


$(TARGET): $(OBJECTS)
	$(CXX) -o $(TARGET) $(BUILDDIR)/$(OBJECTS) -L $(SOURCEDIR)/$(LIBSDIR) -l $(LIBS) -Wl,-R $(SOURCEDIR)/$(LIBSDIR)

$(OBJECTS): $(SOURCEDIR)/$(SRC)
	$(CXX) -c -o $(BUILDDIR)/$(OBJECTS) $(SOURCEDIR)/$(SRC)

clean: 
	rm -rf $(TARGET) $(BUILDDIR)/*.o

run:	
	./$(TARGET)
