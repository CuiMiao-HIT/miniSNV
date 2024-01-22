################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
./src/caller/callerstep1A3.cpp \
./src/caller/main.cpp 

OBJS += \
./src/caller/callerstep1A3.o \
./src/caller/main.o 

CPP_DEPS += \
./src/caller/callerstep1A3.d \
./src/caller/main.d 

# Each subdirectory must supply rules for building sources it contributes
src/caller/%.o: ../src/caller/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I../src/htslib -O3 -std=c++11 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


