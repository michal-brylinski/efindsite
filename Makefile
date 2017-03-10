#===============================================================================
#         ___________.__            .____________.__  __          
#     ____\_   _____/|__| ____    __| _/   _____/|__|/  |_  ____  
#   _/ __ \|    __)  |  |/    \  / __ |\_____  \ |  \   __\/ __ \ 
#   \  ___/|     \   |  |   |  \/ /_/ |/        \|  ||  | \  ___/ 
#    \___  >___  /   |__|___|  /\____ /_______  /|__||__|  \___  >
#        \/    \/            \/      \/       \/               \/ 
#
#                                                  
#   eFindSite - ligand-binding site prediction from meta-threading
#
#   Computational Systems Biology Group
#   Department of Biological Sciences
#   Center for Computation & Technology
#   Louisiana State University
#   407 Choppin Hall, Baton Rouge, LA 70803, USA
#
#   http://www.brylinski.org
#
#   Report bugs to michal@brylinski.org
#
#   Copyright 2013 Michal Brylinski
#
#   This file is part of eFindSite.
#
#   eFindSite is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   eFindSite is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with eFindSite. If not, see <http://www.gnu.org/licenses/>.
#
#===============================================================================

CXX = g++

CC = gcc

FC = gfortran

EXE = efindsite efindsite_screen

SRC_EFINDSITE = src-efindsite
SRC_GZSTREAM  = src-gzstream
SRC_DBSCAN    = src-dbscan
SRC_FRTMALIGN = src-frtmalign
SRC_LIBSVM    = src-libsvm
SRC_NWALIGN   = src-nwalign
SRC_QCPROT    = src-qcprot

CPPFLAGS = -O2 -Wall -Wno-write-strings -fPIC -fopenmp -I. -I$(SRC_EFINDSITE) -I$(SRC_DBSCAN) -I$(SRC_GZSTREAM) -I$(SRC_FRTMALIGN) -I$(SRC_LIBSVM) -I$(SRC_NWALIGN) -I$(SRC_QCPROT)

CCFLAGS = -O2 -Wall -ffast-math -pedantic -std=c99

FFLAGS = -O2

LDFLAGS = -lz -lgfortran -fopenmp -lm -L.

OBJ_EFINDSITE = $(SRC_EFINDSITE)/cluster.o \
                $(SRC_EFINDSITE)/cmps.o \
                $(SRC_EFINDSITE)/coords.o \
                $(SRC_EFINDSITE)/data.o \
                $(SRC_EFINDSITE)/efindsite.o \
                $(SRC_EFINDSITE)/list.o \
                $(SRC_EFINDSITE)/pocket.o \
                $(SRC_EFINDSITE)/runsvm.o \
                $(SRC_EFINDSITE)/tanimoto.o \
                $(SRC_EFINDSITE)/target.o \
                $(SRC_EFINDSITE)/template.o \
                $(SRC_EFINDSITE)/walltime.o \
                $(SRC_DBSCAN)/clusters.o \
                $(SRC_DBSCAN)/dbscan.o \
                $(SRC_DBSCAN)/kdtree2.o \
                $(SRC_DBSCAN)/utils.o \
                $(SRC_FRTMALIGN)/frtmalign.o \
                $(SRC_GZSTREAM)/gzstream.o \
                $(SRC_LIBSVM)/svm.o \
                $(SRC_NWALIGN)/nwalign.o \
                $(SRC_QCPROT)/qcprot.o \
                $(SRC_QCPROT)/rmsd_qcp.o

OBJ_EFINDSITE_SCREEN = $(SRC_EFINDSITE)/efindsite_screen.o \
                       $(SRC_EFINDSITE)/runsvm.o \
                       $(SRC_EFINDSITE)/tanimoto.o \
                       $(SRC_EFINDSITE)/walltime.o \
                       $(SRC_GZSTREAM)/gzstream.o \
                       $(SRC_LIBSVM)/svm.o \

default: $(EXE)

all: $(EXE)

efindsite: $(OBJ_EFINDSITE)
	$(CXX) -o $@ $(OBJ_EFINDSITE) $(LDFLAGS)
	@mkdir -p bin/
	@mv efindsite bin/

efindsite_screen: $(OBJ_EFINDSITE_SCREEN)
	$(CXX) -o $@ $(OBJ_EFINDSITE_SCREEN) $(LDFLAGS)
	@mkdir -p bin/
	@mv efindsite_screen bin/

#=== eFindSite =================================================================

cluster.o: cluster.C
	$(CXX) $(CPPFLAGS) -c -o cluster.o cluster.C

cmps.o: cmps.C
	$(CXX) $(CPPFLAGS) -c -o cmps.o cmps.C

coords.o: coords.C
	$(CXX) $(CPPFLAGS) -c -o coords.o coords.C

data.o: data.C
	$(CXX) $(CPPFLAGS) -c -o data.o data.C

efindsite.o: efindsite.C
	$(CXX) $(CPPFLAGS) -c -o efindsite.o efindsite.C

efindsite_screen.o: efindsite_screen.C
	$(CXX) $(CPPFLAGS) -c -o efindsite_screen.o efindsite_screen.C

list.o: list.C
	$(CXX) $(CPPFLAGS) -c -o list.o list.C

pocket.o: pocket.C
	$(CXX) $(CPPFLAGS) -c -o pocket.o pocket.C

runsvm.o: runsvm.C
	$(CXX) $(CPPFLAGS) -c -o runsvm.o runsvm.C

tanimoto.o: tanimoto.C
	$(CXX) $(CPPFLAGS) -c -o tanimoto.o tanimoto.C

target.o: target.C
	$(CXX) $(CPPFLAGS) -c -o target.o target.C

template.o: target.C
	$(CXX) $(CPPFLAGS) -c -o template.o template.C

walltime.o: walltime.C
	$(CXX) $(CPPFLAGS) -c -o walltime.o walltime.C

#=== DBSCAN ====================================================================

clusters.o: clusters.C
	$(CXX) $(CPPFLAGS) -c -o clusters.o clusters.C

dbscan.o: dbscan.C
	$(CXX) $(CPPFLAGS) -c -o dbscan.o dbscan.C

kdtree2.o: kdtree2.C
	$(CXX) $(CPPFLAGS) -c -o kdtree2.o kdtree2.C

utils.o: utils.C
	$(CXX) $(CPPFLAGS) -c -o utils.o utils.C

#=== fr-TM-align ===============================================================

frtmalign.o: frtmalign.f
	$(FC) $(FFLAGS) -c -o frtmalign.o frtmalign.f

#=== gzstream ==================================================================

gzstream.o: gzstream.C
	$(CXX) $(CPPFLAGS) -c -o gzstream.o gzstream.C

#=== libsvm ====================================================================

svm.o: svm.C
	$(CXX) $(CPPFLAGS) -c -o svm.o svm.C

#=== nwalign ===================================================================

nwalign.o: nwalign.f
	$(FC) $(FFLAGS) -c -o nwalign.o nwalign.f

#=== qcprot ====================================================================

qcprot.o: qcprot.c
	$(CC) $(CCFLAGS) -c -o qcprot.o qcprot.c

rmsd_qcp.o: rmsd_qcp.c
	$(CC) $(CCFLAGS) -c -o rmsd_qcp.o rmsd_qcp.c

#=== clean =====================================================================

clean:
	@(rm -f ${EXE} bin/efindsite bin/efindsite_screen ${OBJ_EFINDSITE} ${OBJ_EFINDSITE_SCREEN})
