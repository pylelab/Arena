CC=clang++
CFLAGS=-O3
LDFLAGS=#-static

MissingRNAatom: MissingRNAatom.cpp MissingRNAatom.hpp Superpose.hpp IdealRNA.hpp PDBParser.hpp pstream.h GeometryTools.hpp cssr.hpp BondLengths.hpp BaseConformation.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
