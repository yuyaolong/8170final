//
//  Face.hpp
//  assignment2
//
//  Created by yuyaolong on 15/10/22.
//  Copyright © 2015年 yuyaolong. All rights reserved.
//

#ifndef __FACE_H__
#define __FACE_H__

class Face {
public:
    int strutIndices[3];
    double vertexAngles[3];
    int pointIndices[3];
    Face(int s[3], double v[3], int p[3]);
private:
    Face();
    
};


#endif 
