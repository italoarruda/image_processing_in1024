
ifeq (, $(findstring CYGWIN, $(shell uname -s)))
  EXE =
  CORE = core*
else
  EXE = .exe
  CORE = *.stackdump
endif


CC        = gcc
CFLAGS    = -Wall -O2
LDFLAGS   = 
LDLIBS    = -lm


TARGET    = main
BIN      := $(addsuffix $(EXE), $(TARGET))

SRCS     := $(wildcard *.c)
OBJS     := $(patsubst %.c, %.o, $(SRCS))
OBJS     += fftsg.o fftsg2d.o


vpath %.c fft2d/


.PHONY: all clean


all: $(TARGET)

$(TARGET): $(OBJS)

clean:
	$(RM) $(BIN) *.o $(CORE)

