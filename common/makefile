rwildcard=$(foreach d,$(wildcard $1*),$(call rwildcard,$d/,$2) $(filter $(subst *, %,$2),$d))
ifeq ($(OS),Windows_NT)
	CHK_DIR_EXISTS = if not exist "$(strip $1)" mkdir "$(strip $1)"
	NUKE = rmdir /s /q
	COPY_DIR = xcopy $1 $2 /E /H /Y
	COPY = xcopy $1 $2
	FIX_PATH = $(subst /,\,$1)
	LIB_SUFFIX = .dll
else
	CHK_DIR_EXISTS = @test -d $1 || mkdir -p $1
	NUKE = rm -r $1
	COPY_DIR = cp -rv $1 $2
	COPY = cp $1 $2
	FIX_PATH = $1
	LIB_SUFFIX = .so
	REPORT = @echo -e "\e[4;1;37m$1\033[0m"
endif
PROJECT_DIR :=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
SRC_DIR := $(PROJECT_DIR)src
LIB_DIR := $(PROJECT_DIR)lib
LIB_NAME := compPhys
LIB_FULL_NAME := lib$(LIB_NAME)$(LIB_SUFFIX)

C_FILES := $(call rwildcard, $(SRC_DIR), *.c)

C_FLAGS += --std=gnu99 -flto -lm -O3 -march=native -I$(PROJECT_DIR)headers

all: lib

lib: $(C_FILES)
	$(call CHK_DIR_EXISTS, $(LIB_DIR))
	$(call REPORT, Building $(LIB_FULL_NAME))
	gcc -fpic -shared $(C_FLAGS) -o $(LIB_DIR)/$(LIB_FULL_NAME) $<

clean:
	$(call REPORT,Cleaning...)
	-$(NUKE) "$(LIB_DIR)"
