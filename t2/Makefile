# ------------------------------------------------------------
# type "make" command in the Linux terminal create report.pdf
#
# type "make clean" to delete all generated files
# ------------------------------------------------------------

all:
	make -C doc
	cp doc/report.pdf .

dataset:
	python2.7 t2_datagen.py > data.txt

clean:
	make -C doc clean

.PHONY: all clean