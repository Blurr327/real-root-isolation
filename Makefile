NAME = real_root_isolation

SRC_DIR = src
OBJ_DIR = obj
INC_DIR = include

CC = gcc
CFLAGS = -I$(INC_DIR) -g
LDFLAGS = -lm -lgmp -lflint

SRC = $(wildcard $(SRC_DIR)/*.c)

OBJ = $(SRC:$(SRC_DIR)/%.c=$(OBJ_DIR)/%.o)

OBJ_NO_MAIN = $(filter-out obj/real_root_isolation.o, $(OBJ))

all: $(NAME)

$(NAME): $(OBJ)
	$(CC) $(CFLAGS) $(OBJ) $(LDFLAGS) -o $(NAME)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -c $< -o $@

test_%: $(OBJ_NO_MAIN)
	$(CC) $(CFLAGS) test/$@.c -o $@ $^ $(LDFLAGS)

vg_%: test_%
	valgrind --leak-check=full ./$<

clean:
	rm -rf $(OBJ_DIR) test_* $(NAME)

# Rebuild complet
re: clean all

.PHONY: all clean re run
