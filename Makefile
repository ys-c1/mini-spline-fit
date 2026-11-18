CC ?= gcc
CFLAGS ?= -std=c11 -Wall -Wextra -Wpedantic -O2
CPPFLAGS ?= -Isrc
LDLIBS ?= -lm

TARGET := test_spline

.PHONY: all test clean

all: $(TARGET)

$(TARGET): tests/test_spline.c src/spline.c src/spline.h
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ tests/test_spline.c src/spline.c $(LDFLAGS) $(LDLIBS)

test: $(TARGET)
	./$(TARGET)

clean:
	rm -f $(TARGET)
