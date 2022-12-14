OUTPUTDIR := bin/

CFLAGS := -std=c++14 -fvisibility=hidden -lpthread -Wall -Wextra

ifeq (,$(CONFIGURATION))
	CONFIGURATION := release
endif

ifeq (debug,$(CONFIGURATION))
CFLAGS += -g
else
CFLAGS += -O2 -fopenmp
endif

HEADERS := src/*.h
SOURCES := src/*.cpp

TARGETBIN := VRP-$(CONFIGURATION)

CXX = mpic++

.SUFFIXES:
.PHONY: all clean

all: $(TARGETBIN)

$(TARGETBIN): $(SOURCES) $(HEADERS)
	$(CXX) -o $@ $(CFLAGS) $(SOURCES) 

format:
	clang-format -i src/*.cpp src/*.h

clean:
	rm -rf ./VRP-$(CONFIGURATION)*

FILES = src/*.cpp \
		src/*.h

handin.tar: $(FILES)
	tar cvf handin.tar $(FILES)
