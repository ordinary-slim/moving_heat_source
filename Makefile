CXX:=g++
MODULE:=MovingHeatSource.so
SRCDIR:=./src
OBJDIR=./obj

SRCS=$(shell find $(SRCDIR) -name *.cpp)
OBJS=$(SRCS:%=$(OBJDIR)/%.o) 
DEPS=$(OBJS:%.o=%.d)

INCLUDE_PATH:=external/pybind11/include/
CFLAGS=-c -fPIC -I$(INCLUDE_PATH) 
LFLAGS=-lpython3.8 -shared

.PHONY: default debug clean run

default: $(MODULE)

debug: CFLAGS += -g
debug: LFLAGS += -g
debug: $(MODULE)

run:
	cd run && python3 run.py

$(MODULE): $(OBJS)
	$(CXX) $(LFLAGS) -o $@ $^

$(OBJDIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) $< -o $@

clean:
	rm -r $(OBJDIR)
	rm $(MODULE)

-include $(DEPS)
