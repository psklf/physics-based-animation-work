using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using System;
using System.IO;

public class FVM : MonoBehaviour
{
  float dt          = 0.003f;
  float mass        = 1;
  float stiffness_0	= 20000.0f;
  float stiffness_1	= 5000.0f;
  float damp        = 0.996f;

  Vector3 kGravity = new Vector3(0, -9.8F, 0);
  Vector3 kFloorNormal = new Vector3(0, 1.0F, 0);

  int[] Tet;
  int tet_number;			// The number of tetrahedra

  Vector3[] 	Force;
  Vector3[] 	V;
  Vector3[] 	X;
  int number;				// The number of vertices

  Matrix4x4[] inv_Dm;

  //For Laplacian smoothing.
  Vector3[] V_sum;
  int[]     V_num;

  SVD svd = new SVD();

  // Helper Function
  Matrix4x4 Matrix4_Sub(Matrix4x4 l, Matrix4x4 rh)
  {
    Matrix4x4 m = Matrix4x4.zero;
    for (int r = 0; r < 4; ++r)
    {
      for (int c = 0; c < 4; ++c)
      {
        m[r, c] = l[r, c] - rh[r, c];
      }
    }
    return m;
  }

  Matrix4x4 Matrix4_Add(Matrix4x4 l, Matrix4x4 rh)
  {
    Matrix4x4 m = Matrix4x4.zero;
    for (int r = 0; r < 4; ++r)
    {
      for (int c = 0; c < 4; ++c)
      {
        m[r, c] = l[r, c] + rh[r, c];
      }
    }
    return m;
  }

  Matrix4x4 Matrix4_Multiply(Matrix4x4 m, float v)
  {
    Matrix4x4 ret = Matrix4x4.zero;
    for (int r = 0; r < 4; ++r)
    {
      for (int c = 0; c < 4; ++c)
      {
        ret[r, c] = m[r, c] * v;
      }
    }
    return ret;
  }

  // Start is called before the first frame update
  void Start()
  {
    // FILO IO: Read the house model from files.
    // The model is from Jonathan Schewchuk's Stellar lib.
    {
      string fileContent = File.ReadAllText("Assets/house2.ele");
      string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);

      tet_number=int.Parse(Strings[0]);
      Tet = new int[tet_number*4];

      for(int tet=0; tet<tet_number; tet++)
      {
        Tet[tet*4+0]=int.Parse(Strings[tet*5+4])-1;
        Tet[tet*4+1]=int.Parse(Strings[tet*5+5])-1;
        Tet[tet*4+2]=int.Parse(Strings[tet*5+6])-1;
        Tet[tet*4+3]=int.Parse(Strings[tet*5+7])-1;
      }
    }
    {
      string fileContent = File.ReadAllText("Assets/house2.node");
      string[] Strings = fileContent.Split(new char[]{' ', '\t', '\r', '\n'}, StringSplitOptions.RemoveEmptyEntries);
      number = int.Parse(Strings[0]);
      X = new Vector3[number];
      for(int i=0; i<number; i++)
      {
        X[i].x=float.Parse(Strings[i*5+5])*0.4f;
        X[i].y=float.Parse(Strings[i*5+6])*0.4f;
        X[i].z=float.Parse(Strings[i*5+7])*0.4f;
      }
      //Centralize the model.
      Vector3 center=Vector3.zero;
      for(int i=0; i<number; i++)		center+=X[i];
      center=center/number;
      for(int i=0; i<number; i++)
      {
        X[i]-=center;
        float temp=X[i].y;
        X[i].y=X[i].z;
        X[i].z=temp;
      }
    }
    /*tet_number=1;
      Tet = new int[tet_number*4];
      Tet[0]=0;
      Tet[1]=1;
      Tet[2]=2;
      Tet[3]=3;

      number=4;
      X = new Vector3[number];
      V = new Vector3[number];
      Force = new Vector3[number];
      X[0]= new Vector3(0, 0, 0);
      X[1]= new Vector3(1, 0, 0);
      X[2]= new Vector3(0, 1, 0);
      X[3]= new Vector3(0, 0, 1);*/


    //Create triangle mesh.
    Vector3[] vertices = new Vector3[tet_number*12];
    int vertex_number=0;
    for(int tet=0; tet<tet_number; tet++)
    {
      vertices[vertex_number++]=X[Tet[tet*4+0]];
      vertices[vertex_number++]=X[Tet[tet*4+2]];
      vertices[vertex_number++]=X[Tet[tet*4+1]];

      vertices[vertex_number++]=X[Tet[tet*4+0]];
      vertices[vertex_number++]=X[Tet[tet*4+3]];
      vertices[vertex_number++]=X[Tet[tet*4+2]];

      vertices[vertex_number++]=X[Tet[tet*4+0]];
      vertices[vertex_number++]=X[Tet[tet*4+1]];
      vertices[vertex_number++]=X[Tet[tet*4+3]];

      vertices[vertex_number++]=X[Tet[tet*4+1]];
      vertices[vertex_number++]=X[Tet[tet*4+2]];
      vertices[vertex_number++]=X[Tet[tet*4+3]];
    }

    int[] triangles = new int[tet_number*12];
    for(int t=0; t<tet_number*4; t++)
    {
      triangles[t*3+0]=t*3+0;
      triangles[t*3+1]=t*3+1;
      triangles[t*3+2]=t*3+2;
    }
    Mesh mesh = GetComponent<MeshFilter> ().mesh;
    mesh.vertices  = vertices;
    mesh.triangles = triangles;
    mesh.RecalculateNormals ();


    V 	  = new Vector3[number];
    Force = new Vector3[number];
    V_sum = new Vector3[number];
    V_num = new int[number];

    // Need to allocate and assign inv_Dm
    inv_Dm = new Matrix4x4[tet_number];
    for (int tet=0; tet<tet_number; tet++)
    {
      inv_Dm[tet] = Build_Edge_Matrix(tet).inverse;
    }
  }

  Matrix4x4 Build_Edge_Matrix(int tet)
  {
    Matrix4x4 ret = Matrix4x4.zero;
    // Need to build edge matrix here.
    Vector3 v0 = X[Tet[tet * 4 + 0]];
    Vector3 v1 = X[Tet[tet * 4 + 1]];
    Vector3 v2 = X[Tet[tet * 4 + 2]];
    Vector3 v3 = X[Tet[tet * 4 + 3]];
    ret.SetColumn(0, v1 - v0);
    ret.SetColumn(1, v2 - v0);
    ret.SetColumn(2, v3 - v0);
    ret.SetColumn(3, new Vector4(0, 0, 0, 1));
    return ret;
  }


  void _Update()
  {
    // Jump up.
    if(Input.GetKeyDown(KeyCode.Space))
    {
      for(int i=0; i<number; i++)
        V[i].y+=0.2f;
    }

    for(int i=0; i<number; i++)
    {
      // Add gravity to Force.
      Force[i] = kGravity * mass;
    }

    //
    for(int tet=0; tet<tet_number; tet++)
    {
      // Deformation Gradient
      Matrix4x4 X_cur = Build_Edge_Matrix(tet);
      Matrix4x4 F = X_cur * inv_Dm[tet];

      // Green Strain
      Matrix4x4 GS = Matrix4_Multiply(Matrix4_Sub(F.transpose * F, Matrix4x4.identity),
                                      0.5F);
      // Second PK Stress
      float trace = GS[0, 0] + GS[1, 1] + GS[2, 2];
      Matrix4x4 S = Matrix4_Add(Matrix4_Multiply(GS, 2.0F * stiffness_1),
                                Matrix4_Multiply(Matrix4x4.identity, stiffness_0 * trace));


      // Elastic Force
      Matrix4x4 P = F * S;
      Matrix4x4 force_123 = Matrix4_Multiply(P * inv_Dm[tet].transpose,
                                             -1.0F / inv_Dm[tet].determinant / 6.0F);

      // Apply Force to vertex
      Vector3 sum = Vector3.zero;
      Force[Tet[tet * 4 + 1]] += (Vector3)force_123.GetColumn(0);
      sum += (Vector3)force_123.GetColumn(0);
      Force[Tet[tet * 4 + 2]] += (Vector3)force_123.GetColumn(1);
      sum += (Vector3)force_123.GetColumn(1);
      Force[Tet[tet * 4 + 3]] += (Vector3)force_123.GetColumn(2);
      sum += (Vector3)force_123.GetColumn(2);
      Force[Tet[tet * 4 + 0]] += -sum;
    }
    //

    for(int i=0; i<number; i++)
    {
      //  Update X and V here.
      V[i] += Force[i] / mass * dt;
      V[i] *= damp;
      X[i] += V[i] * dt;


      //  (Particle) collision with floor.
      // impulse
      //
      float floor_y = -3.0F;
      float dist = (X[i].y - floor_y);
      if (dist < 0)
      {
        V[i] += (Math.Abs(dist) * kFloorNormal) / dt;
        X[i] += Math.Abs(dist) * kFloorNormal;
      }
      //
    }
  }

  // Update is called once per frame
  void Update()
  {
    for(int l=0; l<10; l++)
    {
      _Update();
    }
    // Dump the vertex array for rendering.
    Vector3[] vertices = new Vector3[tet_number*12];
    int vertex_number=0;
    for(int tet=0; tet<tet_number; tet++)
    {
      vertices[vertex_number++]=X[Tet[tet*4+0]];
      vertices[vertex_number++]=X[Tet[tet*4+2]];
      vertices[vertex_number++]=X[Tet[tet*4+1]];
      vertices[vertex_number++]=X[Tet[tet*4+0]];
      vertices[vertex_number++]=X[Tet[tet*4+3]];
      vertices[vertex_number++]=X[Tet[tet*4+2]];
      vertices[vertex_number++]=X[Tet[tet*4+0]];
      vertices[vertex_number++]=X[Tet[tet*4+1]];
      vertices[vertex_number++]=X[Tet[tet*4+3]];
      vertices[vertex_number++]=X[Tet[tet*4+1]];
      vertices[vertex_number++]=X[Tet[tet*4+2]];
      vertices[vertex_number++]=X[Tet[tet*4+3]];
    }
    Mesh mesh = GetComponent<MeshFilter> ().mesh;
    mesh.vertices  = vertices;
    mesh.RecalculateNormals ();
  }
}
