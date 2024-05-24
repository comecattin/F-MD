FC = gfortran
SRC_DIR = source
BIN_DIR = bin

# List of the sources files
SRCS = $(wildcard $(SRC_DIR)/*.f90)

# List of the object files
OBJS = $(patsubst $(SRC_DIR)/%.f90,$(BIN_DIR)/%.o,$(SRCS))

# Main target
md_simulation: $(OBJS)
	$(FC) -o $@ $^

# Rule to compile the object files
$(BIN_DIR)/%.o: $(SRC_DIR)/%.f90 | $(BIN_DIR)
	$(FC) -c -o $@ $<

# Create the bin directory if it doesn't exist
$(BIN_DIR):
	mkdir -p $(BIN_DIR)


.PHONY: clean
clean:
	rm -rf $(BIN_DIR) md_simulation
