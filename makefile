
CC = gcc

CFLAGS = -g -O0 -I./include
LFLAGS = -L./build/lib -lm

LIB = build/lib/libdg1d.a
EXE = build/example

OBJS := build/objs/dg_assemble.o \
	build/objs/dg_basis.o \
	build/objs/spmat.o

exe: lib $(EXE)
	@echo "\033[0;32m=== built exe ===\033[0;0m"

lib: $(LIB)
	@echo "\033[0;32m=== built library ===\033[0;0m"

$(EXE): $(LIB) build/objs/example.o
	@echo "\033[0;36mlinking $@ \033[0;0m"
	@$(CC) build/objs/example.o -o $(EXE) -ldg1d $(LFLAGS)

$(LIB): $(OBJS)
	@echo "\033[0;36mlinking $@ \033[0;0m"

	@mkdir -p build/lib
	@ar rcs $(LIB) $(OBJS)

build/objs/%.o: src/%.c

	@echo "\033[0;35mcompiling $< \033[0;0m"

	@mkdir -p build
	@mkdir -p build/objs
	@$(CC) $(CFLAGS) -c $< -o $@

clean:
	@rm -rf build
	@echo "\033[0;32m=== cleaned build directory ===\033[0;0m"

