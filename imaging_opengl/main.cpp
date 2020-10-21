#include "Headers.h"
#include "convert.h"
#include "InitializingParameters.h"
#include "FunctionDeclarations.h"
#include "FunctionDefinitions.h"
#include "FunctionDefinitions2.h"

int main(int argc,char** argv) {

  // initialize opengl
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_DEPTH|GLUT_RGBA|GLUT_DOUBLE);
  glutInitWindowSize(glutInitWindowSizeX,glutInitWindowSizeY);
  glutCreateWindow("*** PRESS 's' TO GRAB SCREENSHOT ***");
  glutInitWindowPosition(glutInitWindowPositionX,glutInitWindowPositionY);
  
  // input file name & box length
  filename=argv[1];
  L=atof(argv[2]);  
  file.open(filename,ios::in);

  j=0;
  movie=0;
  while (true){
    j++;
    file >> N[j];    
    if( file.eof() ) break;
    ymax=0.0;
    for(i=1; i<=N[j]; i++){
      file >> xread[i][j] >> yread[i][j] >> thetaread[i][j]
	   >> sig1read[i][j] >> sig2read[i][j] >> depth[i][j] >> c[i][j];

      // Alter color of cell based
      R[i][j]=0.0;
      G[i][j]=1.0;
      B[i][j]=1.0;		  

      // sig = diameter in file, radius in imaging program
      sig1read[i][j] = sig1read[i][j]/2.0;
      sig2read[i][j] = sig2read[i][j]/2.0;

      // find maximum distance for output to screen
      if (yread[i][j]>ymax){
	ymax=yread[i][j];
      }      
    }
    cout << j << "  " << N[j] << " " << ymax << endl;
    for(i=1; i<=N[j]; i++){
      xread[i][j]=xread[i][j]/L;      // give box dimensions LxL
      yread[i][j]=yread[i][j]/L; 	  
      sig1read[i][j] = sig1read[i][j]/L;
      sig2read[i][j] = sig2read[i][j]/L;
    }
  }

  // Set # of frames based on input file
  nframes=j;
  
  glutDisplayFunc(drawEllipse2);
  glutReshapeFunc(reshape);
  glutKeyboardFunc(keyboard);
  glutMouseFunc(mouse);
  glutMainLoop();
  
  return 0;
}
