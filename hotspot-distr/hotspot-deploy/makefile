# File: makefile
# Author: Scott Kuehn
# Date Created: Wed Aug 19 22:22:27 PDT 2009

all = dist

CPP_SRCS += \
	./src/Cluster.cpp \
	./src/Hotspot.cpp \
	./src/HotspotDefaults.cpp \
	./src/InputDataReader.cpp \
	./src/MappableCountsDataReader.cpp 
OBJS += \
	./src/Cluster.o \
	./src/Hotspot.o \
	./src/HotspotDefaults.o \
	./src/InputDataReader.o \
	./src/MappableCountsDataReader.o 

GSL = `gsl-config --libs`
LIBS := ${GSL}
BUILDOPTS = -O3 -Wall

RM := rm -rf


dist: prep hotspot

prep:
	mkdir -p bin

hotspot: $(OBJS)
	@echo 'Building target: $@'
	@echo 'Invoking: GCC C++ Linker'
	g++  -o"bin/hotspot" $(OBJS) $(LIBS)
	@echo 'Finished building target: $@'
	@echo ' '

src/%.o: ./src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ ${BUILDOPTS} -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '

clean:
	${RM} bin/hotspot ${OBJS} src/*.d
	-@echo ' '

.PHONY: all clean dependents
.SECONDARY:
