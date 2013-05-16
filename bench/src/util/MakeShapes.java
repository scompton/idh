package util;

import static edu.mines.jtk.util.ArrayMath.*;

public class MakeShapes {
  
  public static float[] cube() {
    float[] v = new float[] {
        -50.0f,-50.0f, 50.0f, // front
        -50.0f, 50.0f, 50.0f,
         50.0f, 50.0f, 50.0f,
         50.0f,-50.0f, 50.0f,
        -50.0f,-50.0f,-50.0f, // back
        -50.0f, 50.0f,-50.0f,
         50.0f, 50.0f,-50.0f,
         50.0f,-50.0f,-50.0f,
        -50.0f, 50.0f, 50.0f, // top
        -50.0f, 50.0f,-50.0f,
         50.0f, 50.0f,-50.0f,
         50.0f, 50.0f, 50.0f,
        -50.0f,-50.0f, 50.0f, // bottom
        -50.0f,-50.0f,-50.0f,
         50.0f,-50.0f,-50.0f,
         50.0f,-50.0f, 50.0f,
         50.0f,-50.0f, 50.0f, // right
         50.0f, 50.0f, 50.0f,
         50.0f, 50.0f,-50.0f,
         50.0f,-50.0f,-50.0f,
        -50.0f,-50.0f, 50.0f, //left
        -50.0f, 50.0f, 50.0f,
        -50.0f, 50.0f,-50.0f,
        -50.0f,-50.0f,-50.0f,
    };
    return v;
  }
  
  public static float[] icosahedron() {
    float[] v = new float[3*3*20];
    int j = 0;
    for (int i=0; i<20; i++) {
      float[] v1 = V_ICOSAHEDRON[T_ICOSAHEDRON[i][0]];
      v[j  ] = v1[0]; 
      v[j+1] = v1[1];
      v[j+2] = v1[2];
      j+=3;
      float[] v2 = V_ICOSAHEDRON[T_ICOSAHEDRON[i][1]];
      v[j  ] = v2[0]; 
      v[j+1] = v2[1];
      v[j+2] = v2[2];
      j+=3;
      float[] v3 = V_ICOSAHEDRON[T_ICOSAHEDRON[i][2]]; 
      v[j  ] = v3[0]; 
      v[j+1] = v3[1];
      v[j+2] = v3[2];
      j+=3;
    }
    return v;
  }
  
  public static float[] sphere() {
    int nSubdivide = 5;
    int nt = 20*(int)pow(4.0f,(float)nSubdivide);
    float[] v = new float[nt*3*3];
    print("nt="+nt+", nv="+v.length);
    int[] j = new int[]{0};
    for (int i=0; i<20; i++) {
      float[] v1 = V_ICOSAHEDRON[T_ICOSAHEDRON[i][0]];
      float[] v2 = V_ICOSAHEDRON[T_ICOSAHEDRON[i][1]];
      float[] v3 = V_ICOSAHEDRON[T_ICOSAHEDRON[i][2]];
      subdivide(v1,v2,v3,0,nSubdivide,v,j);
    }
    return v;
  }

  private static void subdivide(
      float[] v1, float[] v2, float[] v3, int c, int n, float[] v, int[] j) 
  {
    float[] vs = new float[4*3*3];
    float[] v12 = new float[3];
    float[] v23 = new float[3];
    float[] v31 = new float[3];
    for (int i=0; i<3; i++) {
      v12[i] = v1[i]+v2[i];
      v23[i] = v2[i]+v3[i];
      v31[i] = v3[i]+v1[i];
    }
    normalize(v12);
    normalize(v23);
    normalize(v31);
    c++;
    if (c<n) {
      subdivide(v1,v12,v31,c,n,v,j);
      subdivide(v2,v23,v12,c,n,v,j);
      subdivide(v3,v31,v23,c,n,v,j);
      subdivide(v12,v23,v31,c,n,v,j);
    } else {
      int k = 0;
      vs[k  ] = v1[0];
      vs[k+1] = v1[1];
      vs[k+2] = v1[2];
      k += 3;
      vs[k  ] = v12[0];
      vs[k+1] = v12[1];
      vs[k+2] = v12[2];
      k += 3;
      vs[k  ] = v31[0];
      vs[k+1] = v31[1];
      vs[k+2] = v31[2];
      k += 3;
      vs[k  ] = v2[0];
      vs[k+1] = v2[1];
      vs[k+2] = v2[2];
      k += 3;
      vs[k  ] = v23[0];
      vs[k+1] = v23[1];
      vs[k+2] = v23[2];
      k += 3;
      vs[k  ] = v12[0];
      vs[k+1] = v12[1];
      vs[k+2] = v12[2];
      k += 3;
      vs[k  ] = v3[0];
      vs[k+1] = v3[1];
      vs[k+2] = v3[2];
      k += 3;
      vs[k  ] = v31[0];
      vs[k+1] = v31[1];
      vs[k+2] = v31[2];
      k += 3;
      vs[k  ] = v23[0];
      vs[k+1] = v23[1];
      vs[k+2] = v23[2];
      k += 3;
      vs[k  ] = v12[0];
      vs[k+1] = v12[1];
      vs[k+2] = v12[2];
      k += 3;
      vs[k  ] = v23[0];
      vs[k+1] = v23[1];
      vs[k+2] = v23[2];
      k += 3;
      vs[k  ] = v31[0];
      vs[k+1] = v31[1];
      vs[k+2] = v31[2];
      k += 3;
      int nv = vs.length;
      copy(nv,0,vs,j[0],v);
      j[0]+=nv;
    }
  }

  private static void normalize(float[] v) {
    float d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]); 
    if (d == 0.0)
      return;
    v[0] /= d; v[1] /= d; v[2] /= d; 
  }
  
  private static final float X = .525731112119133606f;
  private static final float Z = .850650808352039932f;
  private static final float[][] V_ICOSAHEDRON = new float[][]{
      {-X,0.0f,Z},{X,0.0f,Z}, {-X,0.0f,-Z},{X,0.0f,-Z},    
      {0.0f,Z,X}, {0.0f,Z,-X},{0.0f,-Z,X},{0.0f,-Z,-X},    
      {Z,X,0.0f}, {-Z,X,0.0f},{Z,-X,0.0f},{-Z,-X,0.0f}
  };
  private static final int[][] T_ICOSAHEDRON = new int[][]{ 
      {0, 4, 1}, {0,9, 4}, {9, 5,4}, { 4,5,8}, {4,8, 1},    
      {8,10, 1}, {8,3,10}, {5, 3,8}, { 5,2,3}, {2,7, 3},    
      {7,10, 3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1, 6}, 
      {6, 1,10}, {9,0,11}, {9,11,2}, { 9,2,5}, {7,2,11}
  };

  private static void print(String s) {
    System.out.println(s);
  }
}
