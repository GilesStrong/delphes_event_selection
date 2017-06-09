################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../src/ClassicSVfit-master/src/ClassicSVfit.cc \
../src/ClassicSVfit-master/src/ClassicSVfitIntegrand.cc \
../src/ClassicSVfit-master/src/MeasuredTauLepton.cc \
../src/ClassicSVfit-master/src/SVfitIntegratorMarkovChain.cc \
../src/ClassicSVfit-master/src/svFitAuxFunctions.cc \
../src/ClassicSVfit-master/src/svFitHistogramAdapter.cc 

OBJS += \
./src/ClassicSVfit-master/src/ClassicSVfit.o \
./src/ClassicSVfit-master/src/ClassicSVfitIntegrand.o \
./src/ClassicSVfit-master/src/MeasuredTauLepton.o \
./src/ClassicSVfit-master/src/SVfitIntegratorMarkovChain.o \
./src/ClassicSVfit-master/src/svFitAuxFunctions.o \
./src/ClassicSVfit-master/src/svFitHistogramAdapter.o 

CC_DEPS += \
./src/ClassicSVfit-master/src/ClassicSVfit.d \
./src/ClassicSVfit-master/src/ClassicSVfitIntegrand.d \
./src/ClassicSVfit-master/src/MeasuredTauLepton.d \
./src/ClassicSVfit-master/src/SVfitIntegratorMarkovChain.d \
./src/ClassicSVfit-master/src/svFitAuxFunctions.d \
./src/ClassicSVfit-master/src/svFitHistogramAdapter.d 


# Each subdirectory must supply rules for building sources it contributes
src/ClassicSVfit-master/src/%.o: ../src/ClassicSVfit-master/src/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/usr/local/include/root -I/home/giles/Dropbox/Lisbon/LIP/project_work_environment/ -I/home/giles/programs/Delphes-3.3.2/external -I/usr/local/include -I/home//giles/programs/Delphes-3.3.2 -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


