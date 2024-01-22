################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 

C_SRCS += \
./src/abpoa/abpoa_align.c \
./src/abpoa/abpoa_graph.c \
./src/abpoa/abpoa_plot.c \
./src/abpoa/abpoa_seq.c \
./src/abpoa/abpoa_seed.c \
./src/abpoa/abpoa_output.c \
./src/abpoa/kalloc.c \
./src/abpoa/kstring.c \
./src/abpoa/simd_abpoa_align.c \
./src/abpoa/simd_check.c \
./src/abpoa/utils.c \

OBJS += \
./src/abpoa/abpoa_align.o \
./src/abpoa/abpoa_graph.o \
./src/abpoa/abpoa_plot.o \
./src/abpoa/abpoa_seq.o \
./src/abpoa/abpoa_seed.o \
./src/abpoa/abpoa_output.o \
./src/abpoa/kalloc.o \
./src/abpoa/kstring.o \
./src/abpoa/simd_abpoa_align.o \
./src/abpoa/simd_check.o \
./src/abpoa/utils.o \

C_DEPS += \
./src/abpoa/abpoa_align.d \
./src/abpoa/abpoa_graph.d \
./src/abpoa/abpoa_plot.d \
./src/abpoa/abpoa_seq.d \
./src/abpoa/abpoa_seed.d \
./src/abpoa/abpoa_output.d \
./src/abpoa/kalloc.d \
./src/abpoa/kstring.d \
./src/abpoa/simd_abpoa_align.d \
./src/abpoa/simd_check.d \
./src/abpoa/utils.d \

# Each subdirectory must supply rules for building sources it contributes
src/abpoa/%.o: ../src/abpoa/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -I../src/htslib -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


