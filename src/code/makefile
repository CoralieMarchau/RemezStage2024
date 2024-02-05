SRCS= $(wildcard *.cpp)

BINS=$(SRCS:.cpp=)
LIB_DIR     = /usr/local/lib
#LDLIBS          = -lglpk
O_FILES         = $(SRC_FILES:%.cc=%.o)


CXXFLAGS := $(shell pkg-config --cflags ibex)
LIBS	 := $(shell pkg-config --libs  ibex)


ifeq ($(DEBUG), yes)
CXXFLAGS := $(CXXFLAGS) -O0 -g -pg -Wall 
else
CXXFLAGS := $(CXXFLAGS) -O3 -DNDEBUG 
endif

CXXFLAGS := $(CXXFLAGS) -std=c++11
# -U__STRICT_ANSI__


all: $(BINS)

% :	%.cpp
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -o $@ $< $(LIBS) -lglpk -lsoplex -lz -lgmp -std=c++17

clean:
	rm -f $(BINS)
	$(RM) $(EXEC)
	
