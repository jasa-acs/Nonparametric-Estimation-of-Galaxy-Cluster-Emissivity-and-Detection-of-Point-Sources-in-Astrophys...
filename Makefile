CC := g++
CFLAGS := -std=c++17 -fopenmp -pedantic -Wall -m64 -march=native -mno-sse5
COPTFLAGS := -O3
CLIB := -lstdc++fs -lcfitsio -lCCfits

PROJECTNAME := ALIAS

SOURCEDIR := src
SOURCES := $(shell find $(SOURCEDIR) -name '*.cpp')
HEADERDIR := include
OBJECTDIR := obj
OBJECTS := $(addprefix $(OBJECTDIR)/,$(SOURCES:%.cpp=%.o))

all: $(PROJECTNAME)

rebuild: clean $(PROJECTNAME)

$(PROJECTNAME): $(OBJECTS)
	$(CC) $(OBJECTS) $(CLIB) -fopenmp -o $(PROJECTNAME)

$(OBJECTDIR)/%.o: %.cpp
	mkdir -p $(OBJECTDIR)/$(dir $<) && $(CC) $(CFLAGS) $(COPTFLAGS) -I $(HEADERDIR) -c $< -o $@

.PHONY: clean

clean:
	rm -rf $(OBJECTDIR) $(PROJECTNAME)
