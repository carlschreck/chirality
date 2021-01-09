CC = g++
CFLAGS = -Wall -Wextra -g -O0 -I/usr/local/include -I/usr/include -L/usr/local/lib -lX11 -lm -lGL -lGLU -lglut -Wno-deprecated-declarations  -I/usr/X11R6/include -L/usr/X11R6/lib

CXXFLAGS = -Wall -Wextra -g -O0 

IFLAGS = -I/usr/local/include -I/usr/include

LFLAGS = -L/usr/local/lib -lm -lGL -lGLU -lglut

PROG = opengl

SRCS = main.cpp

#ifeq ($(shell uname),Darwin)
#	LIBS = -framework OpenGL -framework GLUT
#else
#	LIBS = -lglut
#endif

all: $(PROG)

$(PROG):	$(SRCS)
	$(CC) -o $(PROG) $(SRCS) $(CFLAGS)

clean:
	rm -f $(PROG)
	rm *.o
