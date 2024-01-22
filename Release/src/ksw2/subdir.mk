################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 

C_SRCS += \
./src/ksw2/kalloc.c \
./src/ksw2/ksw2_extd2_sse.c \
./src/ksw2/ksw2_extd.c \
./src/ksw2/ksw2_extf2_sse.c \
./src/ksw2/ksw2_exts2_sse.c \
./src/ksw2/ksw2_extz2_sse.c \
./src/ksw2/ksw2_extz.c \
./src/ksw2/ksw2_gg2.c \
./src/ksw2/ksw2_gg2_sse.c \
./src/ksw2/ksw2_gg.c \

OBJS += \
./src/ksw2/kalloc.o \
./src/ksw2/ksw2_extd2_sse.o \
./src/ksw2/ksw2_extd.o \
./src/ksw2/ksw2_extf2_sse.o \
./src/ksw2/ksw2_exts2_sse.o \
./src/ksw2/ksw2_extz2_sse.o \
./src/ksw2/ksw2_extz.o \
./src/ksw2/ksw2_gg2.o \
./src/ksw2/ksw2_gg2_sse.o \
./src/ksw2/ksw2_gg.o \

C_DEPS += \
./src/ksw2/kalloc.d \
./src/ksw2/ksw2_extd2_sse.d \
./src/ksw2/ksw2_extd.d \
./src/ksw2/ksw2_extf2_sse.d \
./src/ksw2/ksw2_exts2_sse.d \
./src/ksw2/ksw2_extz2_sse.d \
./src/ksw2/ksw2_extz.d \
./src/ksw2/ksw2_gg2.d \
./src/ksw2/ksw2_gg2_sse.d \
./src/ksw2/ksw2_gg.d \

# Each subdirectory must supply rules for building sources it contributes
src/ksw2/%.o: ../src/ksw2/%.c
	@echo 'Building file: $<'
	@echo 'Invoking: Cross GCC Compiler'
	gcc -I../src/htslib -O3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


