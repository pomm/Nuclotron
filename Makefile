.DELETE_ON_ERROR:

#ifndef CLASTOOL
#    $(error "Please set the variable CLASTOOL")
#endif

ROOTCONFIG := root-config
ROOTCFLAGS := $(shell $(ROOTCONFIG) --cflags)
ROOTLDFLAGS := $(shell $(ROOTCONFIG) --ldflags)
ROOTGLIBS := $(shell $(ROOTCONFIG) --glibs)

CXX := c++
CXXFLAGS := -O2 -Wall -fPIC $(ROOTCFLAGS)
LD := c++
LDFLAGS := -O2 $(ROOTLDFLAGS) -L./

INCLUDES := 
# -I/u/home/orchen/CLAS/analyzer/analysis_lib/include \
               -I$(CLASTOOL)/include


#LIBS := $(ROOTGLIBS) -lPluto -lstlloader_C
LIBS := $(ROOTGLIBS) 
#\
#               -L$(CLASTOOL)/slib/${OS_NAME} -lClasTool \
#               -L/u/home/orchen/CLAS/analyzer/analysis_lib/slib/ -lTIdentificator
#LIBS := $(ROOTGLIBS) \
               -L/u/home/orchen/CLAS/analyzer/analysis_lib/slib/ -lTIdentificator
#LIBS := $(ROOTGLIBS) \
               -L$(CLASTOOL)/slib/${OS_NAME} -lClasTool \


#FILES := x2_Res
#FILES := CondT_calc
FILES := Generator
#FILES :=Aeep
#FILES :=pSRC
#FILES :=pSRC_E6Corrected

.PHONY: all clean

all: $(FILES)

#x2_Res: x2_Res.o
#CondT_calc: CondT_calc.o
Generator: Generator.o
#Aeep: Aeep.o
#pSRC: pSRC.o
#pSRC_E6Corrected: pSRC_E6Corrected.o



%: %.o
	$(LD) $(LDFLAGS) $(LIBS) -o $@ $^

%.o: %.cxx
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $< -o $@


clean:
	rm -f $(FILES:%=%.o) $(FILES)
