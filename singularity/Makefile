.PHONY = run clean img shell

SINGULARITY = /usr/local/bin/singularity
IMG = trim-galore-0.6.5.img
DEF = image.def

run:
	sudo $(SINGULARITY) run $(IMG)

clean:
	rm -f $(IMG)

img: clean
	sudo $(SINGULARITY) build $(IMG) $(DEF)

shell:
	sudo $(SINGULARITY) shell --writable -B $(shell pwd):/tmp $(IMG)
