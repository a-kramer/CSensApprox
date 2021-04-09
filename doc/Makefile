AUTHOR = "Andrei Kramer <andreikr@kth.se>"
TITLE = "CSensApprox Code Documentation"
CSS = default.css
CODEDOC_OPTIONS =  --author $(AUTHOR) --title $(TITLE)  --docversion 0.1 --copyright "CC BY 4.0" --css $(CSS)

.PHONY: all

all: main.html


main.html: ../main.c
	codedoc $(CODEDOC_OPTIONS) --section "main program" $< > $@

