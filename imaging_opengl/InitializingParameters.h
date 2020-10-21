// OpenGL parameters

const int glutInitWindowSizeX=640;//720;
const int glutInitWindowSizeY=640;//720;
const int glutInitWindowPositionX=1;
const int glutInitWindowPositionY=1;
static GLfloat spin = 1.0;
static GLint bold = 0;
static GLint pbc = 1;
static GLfloat spinstep = 1.0;
int spinstep_last;	


// Program-specific parameters

const int NMAX = 10000;
const int nframesMAX = 4001;

static GLfloat R[NMAX][nframesMAX];
static GLfloat G[NMAX][nframesMAX];
static GLfloat B[NMAX][nframesMAX];

static GLfloat rStep=1.0;
static GLfloat gStep=1.0;
static GLfloat bStep=1.0;

fstream file;
string filename;

int i, j, file_index, N[nframesMAX], nframes;
float phi, delry[nframesMAX];;
float xread[NMAX][nframesMAX];
float yread[NMAX][nframesMAX];
float thetaread[NMAX][nframesMAX];
float sig1read[NMAX][nframesMAX];
float sig2read[NMAX][nframesMAX];
float depth[NMAX][nframesMAX];
float c[NMAX][nframesMAX];

int movie;
double L, ymax;
