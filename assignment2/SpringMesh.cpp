//
//  SpringMesh.cpp
//  assignment2
//
//  Created by yuyaolong on 15/10/22.
//  Copyright © 2015年 yuyaolong. All rights reserved.
//

#include "SpringMesh.h"



#define WEIGHT 10

#define Kij 100
#define KijShear 50
#define Dij 70
#define DijShear 20
#define Lij0 EDGELENGTH


SpringMesh::SpringMesh(){
    
    double ki[STRUTNUMBER] ;
    double di[STRUTNUMBER] ;
    double li[STRUTNUMBER] ;
    int vi[STRUTNUMBER][2] ;
    int fi[STRUTNUMBER][2] ;
    
    
    
//    int si[FACENUMBER][3] = {{2,3,7}, {0,4,3}, {4,5,8}, {1,6,5}};
//    double vai[FACENUMBER][3] = {{45,45,90}, {90,45,45}, {45,45,90}, {90,45,45}};
    //int pi[FACENUMBER][3] = {{0,4,3}, {1,4,0}, {1,5,4}, {2,5,1}};
    int si[FACENUMBER][3] = {{0,0,0}};
    double vai[FACENUMBER][3] = {{0,0,0}};
    int pi[FACENUMBER][3] = {{0,0,0}};
    int com1=0,com2=0,com3=0,topCom=0, downCom=0;
    //faces initioan
    for (unsigned int i=0; i<(LINES-1); i++) {
        for (unsigned int j=0; j<((ROWS-1)*2); j++) {
        int index = i*((ROWS-1)*2)+j;
            //std::cout<<index<<"\n";
            
            if (index%2 == 0) {
                com1 = ROWS*(3*i+1)-2*i-1 + j;
                com2 = com1+1;
                com3 = com2+1;
                topCom = (3*ROWS-2)*(i+1)+j/2;
                downCom = topCom - (3*ROWS-2);
                si[index][0] = com1;
                si[index][1] = com2;
                si[index][2] = topCom;
                
                //std::cout<<si[index][0]<<"  "<<si[index][1]<<"  "<<si[index][2]<<"\n";
                vai[index][0] = 45;
                vai[index][1] = 90;
                vai[index][2] = 45;
                
                pi[index][0] = i*ROWS+j/2;
                pi[index][1] = (i+1)*ROWS + j/2;
                pi[index][2] = (i+1)*ROWS+1 + j/2;
                //std::cout<<pi[index][0]<<"  "<<pi[index][1]<<"  "<<pi[index][2]<<"\n";
                
                ki[com1] = Kij;
                ki[com2] = KijShear;
                ki[com3] = Kij;
                ki[topCom] = Kij;
                ki[downCom] = Kij;
                di[com1] = Dij;
                di[com2] = DijShear;
                di[com3] = Dij;
                di[downCom] = Dij;
                di[topCom] = Dij;
                li[com1] = Lij0;
                li[com2] = Lij0*sqrt(2);
                li[com3] = Lij0;
                li[topCom] = Lij0;
                li[downCom] = Lij0;
                
                vi[com1][0] = pi[index][0];
                vi[com1][1] = pi[index][1];
                vi[com2][0] = pi[index][0];
                vi[com2][1] = pi[index][2];
                vi[topCom][0] = pi[index][1];
                vi[topCom][1] = pi[index][2];
                
                if ( index%((ROWS-1)*2) == 0) {
                    //std::cout<<"line index: "<<index<<"\n";
                    fi[com1][0] = -1;
                    fi[com1][1] = index;
                }
                else
                {
                    fi[com1][0] = index-1;
                    fi[com1][1] = index;
                }
                if ( index%((ROWS-1)*2) ==  ((ROWS-1)*2-2) )
                {
                    //std::cout<<"line index: "<<index<<"\n";
                    fi[com3][0] = index+1;
                    fi[com3][1] = -1;
                }

                
                if (i == 0) {
                    fi[downCom][0] = index+1;
                    fi[downCom][1] = -1;
                }
                else
                {
                    fi[downCom][0] = index+1;
                    fi[downCom][1] = index-(ROWS-1)*2;
                }
                
                if( i == (LINES-2) )
                {
                    //std::cout<<"line index: "<<index<<"\n";
                    fi[topCom][0] = -1;
                    fi[topCom][1] = index;
                }
                
                
                fi[com2][0] = index;
                fi[com2][1] = index+1;
                
                
            }
            else
            {
                si[index][0] = topCom - (3*ROWS-2);
                si[index][1] = com3;
                si[index][2] = com2;
                //std::cout<<si[index][0]<<"  "<<si[index][1]<<"  "<<si[index][2]<<"\n";
                
                vai[index][0] = 45;
                vai[index][1] = 45;
                vai[index][2] = 90;
                
                pi[index][0] = i*ROWS+j/2;
                pi[index][1] = (i+1)*ROWS+1 + j/2;
                pi[index][2] = i*ROWS+1 + j/2;
                //std::cout<<pi[index][0]<<"  "<<pi[index][1]<<"  "<<pi[index][2]<<"\n";
                
                vi[com3][0] = pi[index][2];
                vi[com3][1] = pi[index][1];
                vi[downCom][0] = pi[index][0];
                vi[downCom][1] = pi[index][2];

                
            }
            //std::cout<<si[index][0]<<"  "<<si[index][1]<<"  "<<si[index][2]<<"  "<<"\n";
        faces.push_back(Face(si[index], vai[index], pi[index]));
    }
}
    
    
    for (int i=0; i<STRUTNUMBER; i++) {
        //std::cout<<fi[i][0]<<" "<<fi[i][1]<<"\n";
        //std::cout<<li[i]<<"\n";
        
         struts.push_back(Strut(ki[i], di[i], li[i], vi[i], fi[i]));
        
    }
    
    
    double mi[VERTEXNUMBER];
    Vector3d Xi[VERTEXNUMBER];
    Vector3d Vi[VERTEXNUMBER];
    Vector3d Fi[VERTEXNUMBER];
    for (int i=0; i<LINES; i++) {
        for (int j=0; j<ROWS; j++) {
            //std::cout<<j*EDGELENGTH<<"  "<<i*EDGELENGTH<<"\n";
            int index = i*ROWS+j;
            if ((0==index)||(ROWS-1)==index||((VERTEXNUMBER-1)==index)||(VERTEXNUMBER-ROWS) == index) {
                mi[index] = WEIGHT*1.0/4;
            }else if ((index%ROWS == 0)|| (index%ROWS == (ROWS-1)) || (index<ROWS) || (index>(VERTEXNUMBER-ROWS)))
            {
                mi[index] = WEIGHT*1.0/2;
            }else
            {
                mi[index] = WEIGHT;
            }
            //std::cout<<mi[index]<<" ";
            
            
            Xi[index] = Vector3d(j*EDGELENGTH-((ROWS-1)*EDGELENGTH)/2, MESHHEIGHT, i*EDGELENGTH-((LINES-1)*EDGELENGTH)/2);
            Vi[index] = Vector3d(0,0,0);
            Fi[index] = Vector3d(0,0,0);
            points.push_back(Point(mi[index], Xi[index], Vi[index], Fi[index]));
            //std::cout<<Xi[index]<<"\n";
        }
        
    }
    
    
    
   
    
    
}



