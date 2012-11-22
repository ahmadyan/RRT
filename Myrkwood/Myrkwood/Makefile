CXXFLAGS =	-O2 -g -Wall -fmessage-length=0 -I/usr/local/include -I/usr/include -I/usr/include/gsl

OBJS =	utility.o sampling.o point.o config.o VoronoiDiagramGenerator.o plot.o main.o

LIBS =  -lgsl -lgslcblas -lm  -L/usr/lib

TARGET =	Reachability

$(TARGET):	$(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS) $(LIBS) 

all:	$(TARGET)

clean:
	rm -f $(OBJS) $(TARGET)
