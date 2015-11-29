// File: NeheTGATexture.h
// Copyright 2006 GameDev.net LLC, Carsten Haubold
//
// Title: A class for loading TGA-textures
//
// Version History:
//		12 July 06	-	created
//
#ifndef _NeHe_TGA
#define _NeHe_TGA

//include Windows (for VC++)
#ifdef _WIN32
	#include <windows.h>
	#define WIN32_LEAN_AND_MEAN		// trim the excess fat from Windows
	#define WIN32_EXTRA_LEAN
#endif

// OpenGL headers.
#include <OpenGL/gl.h>
#include <OpenGL/glu.h>

#include <iostream>
#include <fstream>
#include <sstream>



using namespace std;


	//Class: imageTGA
	//a class which loads a TGA file and binds it to an OpenGL texture ID
	class imageTGA{
		public:
			//Function: imageTGA
			//empty constructor
			imageTGA();

			//Function: imageTGA
			//constructor which loads the specified TGA file
			//
			//Parameters:
			//	filename - full or relative path to the texture, has to be .tga
			//	mipmaps - whether to create Mipmaps(using gluBuild2DMipmaps) or not
			//	texLog - a logfile which contains the messages from texture loading


			//Function: ~imageTGA
			//empty destructor
			~imageTGA();

			//Function: getId
			//Returns the ID of the texture which is used to access it under OpenGL
			//
			//Returns:
			//	ID of the texture, 0 if an error occured
			GLuint getId();
			bool loadTGA(const string &filename);

			//Function: hasAlpha
			//Returns true if the texture image contains an alpha channel
			//
			//Returns:
			//	true if texture has alpha channel
			bool hasAlpha();
		private:
			//The width of the texture image
			int m_width;
			//height of the texture image
			int m_height;
			//bits per pixel
        		int m_bpp;
			//the ID
			GLuint m_id;
	};



#endif
