################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../01_modScores.cpp \
../Chromosome.cpp \
../DictionaryK.cpp \
../FileReader.cpp \
../GeneticMap.cpp \
../ModScores.cpp \
../Partition.cpp \
../VCF.cpp 

OBJS += \
./01_modScores.o \
./Chromosome.o \
./DictionaryK.o \
./FileReader.o \
./GeneticMap.o \
./ModScores.o \
./Partition.o \
./VCF.o 

CPP_DEPS += \
./01_modScores.d \
./Chromosome.d \
./DictionaryK.d \
./FileReader.d \
./GeneticMap.d \
./ModScores.d \
./Partition.d \
./VCF.d 


# Each subdirectory must supply rules for building sources it contributes
%.o: ../%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


