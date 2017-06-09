################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CC_SRCS += \
../src/SVfit_standalone-master/src/LikelihoodFunctions.cc \
../src/SVfit_standalone-master/src/SVfitStandaloneAlgorithm.cc \
../src/SVfit_standalone-master/src/SVfitStandaloneLikelihood.cc \
../src/SVfit_standalone-master/src/SVfitStandaloneMarkovChainIntegrator.cc \
../src/SVfit_standalone-master/src/svFitStandaloneAuxFunctions.cc 

OBJS += \
./src/SVfit_standalone-master/src/LikelihoodFunctions.o \
./src/SVfit_standalone-master/src/SVfitStandaloneAlgorithm.o \
./src/SVfit_standalone-master/src/SVfitStandaloneLikelihood.o \
./src/SVfit_standalone-master/src/SVfitStandaloneMarkovChainIntegrator.o \
./src/SVfit_standalone-master/src/svFitStandaloneAuxFunctions.o 

CC_DEPS += \
./src/SVfit_standalone-master/src/LikelihoodFunctions.d \
./src/SVfit_standalone-master/src/SVfitStandaloneAlgorithm.d \
./src/SVfit_standalone-master/src/SVfitStandaloneLikelihood.d \
./src/SVfit_standalone-master/src/SVfitStandaloneMarkovChainIntegrator.d \
./src/SVfit_standalone-master/src/svFitStandaloneAuxFunctions.d 


# Each subdirectory must supply rules for building sources it contributes
src/SVfit_standalone-master/src/%.o: ../src/SVfit_standalone-master/src/%.cc
	@echo 'Building file: $<'
	@echo 'Invoking: Cross G++ Compiler'
	g++ -std=c++0x -I/home/cms/giles/programs/include/root -I/home/cms/giles/programs/Delphes-3.3.2/external -I/home/cms/giles/programs/include -I/home/cms/giles/programs/Delphes-3.3.2 -O0 -g3 -Wall -c -fmessage-length=0 -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


