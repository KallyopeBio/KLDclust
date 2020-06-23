################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../Penalty_Class/DictionaryPenalty.cpp \
../Penalty_Class/PwrChiSq.cpp \
../Penalty_Class/ThresholdPenalty.cpp 

OBJS += \
./Penalty_Class/DictionaryPenalty.o \
./Penalty_Class/PwrChiSq.o \
./Penalty_Class/ThresholdPenalty.o 

CPP_DEPS += \
./Penalty_Class/DictionaryPenalty.d \
./Penalty_Class/PwrChiSq.d \
./Penalty_Class/ThresholdPenalty.d 


# Each subdirectory must supply rules for building sources it contributes
Penalty_Class/%.o: ../Penalty_Class/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -std=c++0x -O3 -Wall -c -fmessage-length=0 -fopenmp -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


