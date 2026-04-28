# Standalone terra C++ library + small demo executable.
#   make          — build build/libterra.a and ./terra
#   make lib      — library only
#   make clean
#
# Requires: g++, gdal-config, pkg-config proj, same stack as the R package.

CXX      ?= g++
BUILD    := build
GDALDATADIR   := $(shell gdal-config --datadir)
PROJDATADIR   := $(shell pkg-config --variable=datadir proj 2>/dev/null)
GDALPLUGINDIR := $(or $(shell gdal-config --plugindir 2>/dev/null),$(shell gdal-config --prefix)/lib/gdalplugins)

# Quoted -D values so paths with spaces survive the preprocessor.
CPPFLAGS := -std=c++17 -Dstandalone $(shell gdal-config --cflags) -I../src/ \
	'-DTERRA_GDAL_DATA="$(GDALDATADIR)"' \
	'-DTERRA_PROJ_DATA="$(PROJDATADIR)"' \
	'-DTERRA_GDAL_PLUGINDIR="$(GDALPLUGINDIR)"'

# Every ../src/*.cpp except the R → Rcpp bridge (not used with -Dstandalone).
TERRA_EXCL := RcppExports.cpp RcppFunctions.cpp RcppModule.cpp
ALL_SRC_CPPS := $(notdir $(wildcard ../src/*.cpp))
TERRA_SRCS := $(filter-out $(TERRA_EXCL),$(sort $(ALL_SRC_CPPS)))

TERRA_OBJS := $(TERRA_SRCS:%.cpp=$(BUILD)/%.o)

APP_SRCS := gdal_init.cpp main.cpp show.cpp
APP_OBJS := $(APP_SRCS:%.cpp=$(BUILD)/%.o)

LIBA := $(BUILD)/libterra.a
BIN  := terra

.PHONY: all lib app clean

all: $(BIN)

lib: $(LIBA)

$(BIN): $(APP_OBJS) $(LIBA)
	$(CXX) -o $@ $(APP_OBJS) $(LIBA) \
		$(shell gdal-config --libs) $(shell gdal-config --dep-libs) -lgeos_c

$(LIBA): $(TERRA_OBJS)
	ar rcs $@ $(TERRA_OBJS)

app: $(BIN)

$(BUILD)/%.o: ../src/%.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) -c -o $@ $<

$(BUILD)/%.o: %.cpp | $(BUILD)
	$(CXX) $(CPPFLAGS) -I. -c -o $@ $<

$(BUILD):
	mkdir -p $(BUILD)

clean:
	rm -f $(BIN) $(LIBA) $(TERRA_OBJS) $(APP_OBJS)
	-rmdir $(BUILD) 2>/dev/null || true
