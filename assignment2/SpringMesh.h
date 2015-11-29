//
//  SpringMesh.hpp
//  assignment2
//
//  Created by yuyaolong on 15/10/22.
//  Copyright © 2015年 yuyaolong. All rights reserved.
//

#ifndef __SPRINGMESH_H__
#define __SPRINGMESH_H__

#include <vector>
#include "Vector.h"
#include "Strut.h"
#include "Face.h"
#include "Point.h"

#define EDGELENGTH 1.2
#define LINES 30
#define ROWS 30

const unsigned int STRUTNUMBER = (LINES-1)*ROWS + LINES*(ROWS-1) + (LINES-1)*(ROWS-1);
const unsigned int FACENUMBER = (LINES-1)*(ROWS-1)*2;
const unsigned int VERTEXNUMBER = LINES*ROWS;

const int MESHHEIGHT = -3;

class SpringMesh {
public:
    std::vector<Strut> struts;
    std::vector<Face> faces;
    std::vector<Point> points;
    SpringMesh();
    
};




#endif
