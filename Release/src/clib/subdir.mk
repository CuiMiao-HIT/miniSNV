################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
C_SRCS += \
../src/clib/utils.c
# ../src/clib/bam_file.c \
# ../src/clib/vcf_file.c 

OBJS += \
./src/clib/utils.o 
# ./src/clib/bam_file.o \
# ./src/clib/vcf_file.o 

C_DEPS += \
./src/clib/utils.d 
# ./src/clib/bam_file.d \
# ./src/clib/vcf_file.d 

# Each subdirectory must supply rules for building sources it contributes
src/clib/%.o: ../src/clib/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -I../src/htslib -O3 -std=c99 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


