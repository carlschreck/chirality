//-------------------------------------------------------------
/// \file	tgasave.h
/// \date	9-feb-2005
/// \author	Rob Bateman
/// \brief	writes a 24 or 32 bit tga file to the disk.
//-------------------------------------------------------------

bool SaveScreenGrab(const char* filename);

#ifndef TGA_SAVE__H__
#define TGA_SAVE__H__

bool WriteTga(const char* filename,
	      const unsigned w,
	      const unsigned h,
	      const unsigned bpp,
	      const unsigned char* pixels);

#endif


/// the tga header
#pragma pack(push,1)
struct TGAHeader
{
  // sometimes the tga file has a field with some custom info in. This 
  // just identifies the size of that field. If it is anything other
  // than zero, forget it.
  unsigned char m_iIdentificationFieldSize;
  
  // This field specifies if a colour map is present, 0-no, 1 yes...
  unsigned char m_iColourMapType;
  
  // only going to support RGB/RGBA/8bit - 2, colour mapped - 1
  unsigned char m_iImageTypeCode;
  
  // ignore this field....0
  unsigned short m_iColorMapOrigin;
  
  // size of the colour map
  unsigned short m_iColorMapLength;
  
  // bits per pixel of the colour map entries...
  unsigned char m_iColourMapEntrySize;
  
  // ignore this field..... 0
  unsigned short m_iX_Origin;
  
  // ignore this field..... 0
  unsigned short m_iY_Origin;
  
  // the image width....
  unsigned short m_iWidth;
  
  // the image height.... 
  unsigned short m_iHeight;
  
  // the bits per pixel of the image, 8,16,24 or 32
  unsigned char m_iBPP;
  
  // ignore this field.... 0
  unsigned char m_ImageDescriptorByte;
};
#pragma pack(pop)

bool WriteTga(const char* filename,
	      const unsigned w,
	      const unsigned h,
	      const unsigned bpp,
	      const unsigned char* pixels) {
  
  // a flag to see if writing 32 bit image
  bool rgba=false;
  
  // make sure format is valid (allows bits or bytes per pixel)
  switch(bpp) {
  case 4: case 32:
    rgba = true;
    break;
    
  case 3: case 24:
    rgba = false;
    break;
    
  default:
    std::cerr << "[ERROR] Unsupported bit depth requested" << std::endl;
    return false;
  }
  
  FILE* fp = fopen(filename,"wb");
  
  if(!fp) {
    std::cerr << "[ERROR] could not save file \"" 
	      << filename << "\"" << std::endl;
    return false;
  }
  
  // fill the file header with correct info....
  TGAHeader header;
  
  // wipe to 0
  memset(&header,0,sizeof(TGAHeader));
  
  // rgb or rgba image
  header.m_iImageTypeCode = 2;
  
  // set image size
  header.m_iWidth = w;
  header.m_iHeight = h;
  
  // set bits per pixel
  header.m_iBPP = rgba ? 32 : 24;
  
  // write header as first 18 bytes of output file
  fwrite(&header,sizeof(TGAHeader),1,fp);
  
  // get num pixels
  unsigned int total_size = w * h * (3 + (rgba?1:0));
  unsigned int this_pixel_start = 0;
  
  // loop through each pixel
  for( ; this_pixel_start != total_size; this_pixel_start += 3 ) {
    
    // get address of pixel data
    const unsigned char* pixel = pixels + this_pixel_start;
    
    // write as BGR
    fputc(pixel[2],fp);
    fputc(pixel[1],fp);
    fputc(pixel[0],fp);
    
    // if RGBA also write alpha value
    if (rgba) {
      fputc(pixel[3],fp);
      ++this_pixel_start;
    }
  }
  
  fclose(fp);
  
  return true;		  
}


//------------------------------------------------------------	Reshape()
void Reshape(int w, int h)
{
  if (h==0) {
    h=1;
  }
  glViewport(0,0,w,h);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  gluPerspective(45,(float)w/h,0.1,100);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
}

//------------------------------------------------------------	Draw()
//
void Draw() {
  
  glClear(GL_DEPTH_BUFFER_BIT|GL_COLOR_BUFFER_BIT);
  
  glLoadIdentity();
  gluLookAt(1,2,6,0,0,0,0,1,0);
  
  glutSolidTeapot(1);
  
  glutSwapBuffers();
}


#ifndef SCREEN_GRAB__H__
#define SCREEN_GRAB__H__

bool SaveScreenGrab(const char* filename) {
  
  // get some info about the screen size
  unsigned sw      =  glutGet(GLUT_WINDOW_WIDTH);
  unsigned sh      =  glutGet(GLUT_WINDOW_HEIGHT);
  unsigned bpp     =  glutGet(GLUT_WINDOW_RGBA) ? 4 : 3;
  GLenum   format  =  (bpp==4) ? GL_RGBA : GL_RGB;
  
  
  // allocate memory to store image data
  unsigned char* pdata = new unsigned char[sw*sh*bpp];
  
  // read from front buffer
  glReadBuffer(GL_FRONT);
  
  // read pixel data
  glReadPixels(0,0,sw,sh,format,GL_UNSIGNED_BYTE,pdata);
  
  // write data as a tga file
  bool ret = WriteTga(filename,sw,sh,bpp,pdata);
  
  // clean up
  delete [] pdata;
  
  // done
  return ret;
}

#endif
