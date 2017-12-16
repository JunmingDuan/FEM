BIN = $(patsubst %.eps, %.pdf, $(wildcard *.eps))

all: $(BIN)
$(BIN) : %.pdf : %.eps
	epstopdf $<
.PHONY: clean
clean:
	-rm $(BIN)
