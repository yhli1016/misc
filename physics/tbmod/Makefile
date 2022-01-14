.PHONY: core

core: core.pyx
	python setup.py build_ext --inplace
	rm -rf build core.c
