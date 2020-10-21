CC = g++
CFLAGS = -Wall -I/usr/X11R6/include -L/usr/X11R6/lib -lX11 -Wno-deprecated-declarations
PROG = opengl

SRCS = main.cpp

ifeq ($(shell uname),Darwin)
	LIBS = -framework OpenGL -framework GLUT
else
	LIBS = -lglut
endif

all: $(PROG)

$(PROG):	$(SRCS)
	$(CC) $(CFLAGS) -o $(PROG) $(SRCS) $(LIBS)

clean:
	rm -f $(PROG)
	rm *.o
