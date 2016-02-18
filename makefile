PROJECT_DIR :=$(dir $(realpath $(lastword $(MAKEFILE_LIST))))
SUBDIRS :=$(wildcard $(PROJECT_DIR)/chapter_??)
REPORT = @echo -e "\e[4;1;37m$1\033[0m"
all:
	$(foreach subdir, $(SUBDIRS), $(MAKE) -C $(subdir);)

clean:
	$(foreach subdir, $(SUBDIRS), $(MAKE) -C $(subdir) clean;)
	make -C $(PROJECT_DIR)/common clean
