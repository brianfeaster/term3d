CXXFLAGS = -std=c++14
CXX_ANALYZER?= c++-analyzer

all: asc++clock

clean::
	rm -f asc++clock

lint::
	$(CXX_ANALYZER) $(CXXFLAGS) -fsyntax-only asc++clock.cc
