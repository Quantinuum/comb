# https on the safe side in case ssh not configured

M4RI_REPO = https://github.com/malb/m4ri.git
BRANCH = release-20240729
TMPDIR := $(shell mktemp -d)
ROOT_DIR := $(dir $(realpath $(lastword $(MAKEFILE_LIST))))
M4RI_CONFIGURE_OPTIONS = --with-pic --prefix=$(ROOT_DIR)/m4ri --disable-png --disable-dependency-tracking

m4ri:
	cd $(TMPDIR) && git clone $(M4RI_REPO)  -b $(BRANCH) &&	cd m4ri/ && autoreconf --install && ./configure	$(M4RI_CONFIGURE_OPTIONS) && $(MAKE) && $(MAKE) check && $(MAKE) install

m4ri-clean:
	rm -rf m4ri/
