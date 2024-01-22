################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
./src/cmtools/haplotype.cpp \
./src/cmtools/math_func.cpp \
./src/cmtools/load_bam.cpp \
./src/cmtools/mantaAssembler.cpp \
./src/cmtools/edlib.cpp \
./src/cmtools/form.cpp \

OBJS += \
./src/cmtools/haplotype.o \
./src/cmtools/math_func.o \
./src/cmtools/load_bam.o \
./src/cmtools/mantaAssembler.o \
./src/cmtools/edlib.o \
./src/cmtools/form.o \

CPP_DEPS += \
./src/cmtools/haplotype.d \
./src/cmtools/math_func.d \
./src/cmtools/load_bam.d \
./src/cmtools/mantaAssembler.d \
./src/cmtools/edlib.d \
./src/cmtools/form.d \

# Each subdirectory must supply rules for building sources it contributes
src/cmtools/%.o: ../src/cmtools/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -I../src/htslib -O3 -std=c++11 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


