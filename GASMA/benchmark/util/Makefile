FOLDER_ROOT=../..
FOLDER_BUILD=../../build
FOLDER_BIN=../../bin

UTILS=generate_dataset
UTILS_SRC=$(addsuffix .c, $(UTILS))

OBJS=$(FOLDER_BUILD)/*.o

all: LIBS+=$(LD_FLAGS)
all: FLAGS=$(CC_FLAGS)
all: $(UTILS)

generate_dataset: generate_dataset.c
	$(CC) $(CC_FLAGS) generate_dataset.c -o $(FOLDER_BIN)/generate_dataset $(LD_FLAGS)