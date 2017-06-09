################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../src/ClassicSVfit-master/bin/testClassicSVfit.cc 

OBJS += \
./src/ClassicSVfit-master/bin/testClassicSVfit.o 

CC_DEPS += \
./src/ClassicSVfit-master/bin/testClassicSVfit.d 


# Each subdirectory must supply rules for building sources it contributes
src/ClassicSVfit-master/bin/%.o: ../src/ClassicSVfit-master/bin/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/include/root -I/home/giles/Dropbox/Lisbon/LIP/project_work_environment/ -I/home/giles/programs/Delphes-3.3.2/external -I/usr/local/include -I/home//giles/programs/Delphes-3.3.2 -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


