CC	=	g++ -Wall -std=c++11 -O3 
UNAME = $(shell uname)
BINDIR = bin
OBJDIR = .obj
SRCDIR = src
SOURCESXX := $(wildcard $(SRCDIR)/*.cxx)
SOURCES  := $(wildcard $(SRCDIR)/*.cpp)
OBJECTS  := $(SOURCES:$(SRCDIR)/%.cpp=$(OBJDIR)/%.o)
EXECUTABLES := $(SOURCESXX:$(SRCDIR)/%.cxx=$(BINDIR)/%)
OPENCV = `pkg-config --libs opencv`

ifeq ($(UNAME), Darwin)
   LIBS = -L/opt/local/lib -L/usr/lib -lgsl -lblas -I/opt/local/include 
   FRAMEW = -framework OpenGL
else
   LIBS	=  -L/usr/lib -lgsl -lblas -l armadillo $(OPENCV) -fopenmp -lglut -lGL 
   INCS = -I /usr/local/include -I/usr/include  
endif

all: $(EXECUTABLES)

$(OBJECTS): $(OBJDIR)/%.o : $(SRCDIR)/%.cpp 
	$(CC) $(CFLAGS) -c $< -o $@ $(LIBS) $(FRAMEW) $(INCS) 
$(EXECUTABLES): $(BINDIR)/% : $(SRCDIR)/%.cxx $(OBJECTS) 
	$(CC) -o $@ $< $(OBJECTS) $(LIBS) $(FRAMEW) $(INCS) 

.PHONY: all clear

clean:
	rm $(EXECUTABLES) 
	rm $(OBJECTS)


