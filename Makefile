CXX:=g++
EXE:=main
SRCDIR:=./src
OBJDIR=./obj

SRCS=$(shell find $(SRCDIR) -name *.cpp)
OBJS=$(SRCS:%=$(OBJDIR)/%.o) 
DEPS=$(OBJS:%.o=%.d)

CFLAGS=-I$(INCLUDE_PATH) -M -MMD -MP
LFLAGS=-lpython3.8#matplotlib

.PHONY: default debug clean run

default: $(EXE)

debug: CFLAGS += -g
debug: LFLAGS += -g
debug: $(EXE)

run:
	./$(EXE)

$(EXE): $(OBJS)
	$(CXX) -o $@ $^ $(LFLAGS)

$(OBJDIR)/%.cpp.o: %.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CFLAGS) -c $< -o $@

clean:
	rm -r $(OBJDIR)
	rm $(EXE)

-include $(DEPS)
