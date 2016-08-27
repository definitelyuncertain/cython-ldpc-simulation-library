all : libLDPC

PY=python
CFLAGS=build_ext -i
CUSTOM_LIB_PATH=/usr/local/lib

libLDPC : setup_LDPC.py
	$(PY) setup_LDPC.py $(CFLAGS)