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
  float damp        = 0.999f;

  float restitution   = 6.0f;					// for collision
  float restitution_t = 0.3f;					// for collision

  Vector3 kGravity = new Vector3(0, -9.8F, 0);
  Vector3 kFloorNormal = new Vector3(0, 1.0F, 0);

  // use this switch between green strain method and SVD method
  bool use_green_strain_pk = false;

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

  float GetdWdLambda(float I_c, float II_c, float III_c, float lambda)
  {
    // neo-Hookean
    // float ret = (-1.0F/3.0F * stiffness_0 * I_c * Mathf.Pow(III_c, -4.0F / 3.0F) -
    //            0.5F * stiffness_1 * Mathf.Pow(III_c, -1.5F)) * 4.0F * Mathf.Pow(lambda, 3) +
    // stiffness_0 * Mathf.Pow(III_c, -1.0F / 3.0F) * 2.0F * lambda;


    // W = s0/2[I_c - 3]^2 + s1/4[III_c - 2I_c + 3]
    float ret = stiffness_1 * 0.25F * 4.0F * Mathf.Pow(lambda, 3) +
      (stiffness_0 * (I_c - 3) - 0.5F * stiffness_1) * 2.0F * lambda;

    return ret;
  }

  Matrix4x4 GetFirstPKStress(Matrix4x4 F)
  {
    if (use_green_strain_pk)
    {
      // Green Strain
      Matrix4x4 GS = Matrix4_Multiply(Matrix4_Sub(F.transpose * F, Matrix4x4.identity),
          0.5F);
      // Second PK Stress
      float trace = GS[0, 0] + GS[1, 1] + GS[2, 2];
      Matrix4x4 S = Matrix4_Add(Matrix4_Multiply(GS, 2.0F * stiffness_1),
          Matrix4_Multiply(Matrix4x4.identity, stiffness_0 * trace));

      // Elastic Force
      Matrix4x4 P = F * S;
      return P;
    }
    else
    {
      // SVD
      Matrix4x4 U = Matrix4x4.zero;
      Matrix4x4 S = Matrix4x4.zero;
      Matrix4x4 V_t = Matrix4x4.zero;
      svd.svd(F, ref U, ref S, ref V_t);

      float I_c = Mathf.Pow(S[0, 0], 2) + Mathf.Pow(S[1, 1], 2) +
        Mathf.Pow(S[2, 2], 2);
      float III_c = Mathf.Pow(S[0, 0], 4) + Mathf.Pow(S[1, 1], 4) +
        Mathf.Pow(S[2, 2], 4);

      Matrix4x4 diag = Matrix4x4.zero;
      for (int i = 0; i < 3; ++i) { diag[i, i] = GetdWdLambda(I_c, 0, III_c, S[i, i]); }

      // StVk
      /*
      float l0l0 = Mathf.Pow(S[0, 0], 2);
      float l1l1 = Mathf.Pow(S[1, 1], 2);
      float l2l2 = Mathf.Pow(S[2, 2], 2);
      float I_c = Mathf.Pow(S[0, 0], 2) + Mathf.Pow(S[1, 1], 2) +
                  Mathf.Pow(S[2, 2], 2);
      float II_c = l0l0 * l1l1 + l0l0 * l2l2 + l1l1 * l2l2;
      float d00 = stiffness_1 / 2.0F * (l1l1 * S[0, 0] + l2l2 * S[0, 0]) +
        2.0F * S[0, 0] * (stiffness_0 * I_c - 3.0F * stiffness_0 - 0.5F * stiffness_1 * I_c);
      float d11 = stiffness_1 / 2.0F * (l0l0 * S[1, 1] + l2l2 * S[1, 1]) +
        2.0F * S[1, 1] * (stiffness_0 * I_c - 3.0F * stiffness_0 - 0.5F * stiffness_1 * I_c);
      float d22 = stiffness_1 / 2.0F * (l1l1 * S[2, 2] + l0l0 * S[2, 2]) +
        2.0F * S[2, 2] * (stiffness_0 * I_c - 3.0F * stiffness_0 - 0.5F * stiffness_1 * I_c);

      diag[0, 0] = d00;
      diag[1, 1] = d11;
      diag[2, 2] = d22;
      */

      Matrix4x4 P = U * diag * V_t.transpose;
      return P;
    }
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

      // Reset
      V_sum[i] = Vector3.zero;
      V_num[i] = 0;
    }

    for(int tet=0; tet<tet_number; tet++)
    {
      // Deformation Gradient
      Matrix4x4 X_cur = Build_Edge_Matrix(tet);
      Matrix4x4 F = X_cur * inv_Dm[tet];

      Matrix4x4 P = GetFirstPKStress(F);

      /*
      // Green Strain
      Matrix4x4 GS = Matrix4_Multiply(Matrix4_Sub(F.transpose * F, Matrix4x4.identity),
                                      0.5F);
      // Second PK Stress
      float trace = GS[0, 0] + GS[1, 1] + GS[2, 2];
      Matrix4x4 S = Matrix4_Add(Matrix4_Multiply(GS, 2.0F * stiffness_1),
                                Matrix4_Multiply(Matrix4x4.identity, stiffness_0 * trace));


      // Elastic Force
      Matrix4x4 P = F * S;
      */
      Matrix4x4 force_123 = Matrix4_Multiply(P * inv_Dm[tet].transpose,
                                             -1.0F / inv_Dm[tet].determinant / 6.0F);

      // Apply Force to vertex
      Vector3 minus_f0 = Vector3.zero;
      Force[Tet[tet * 4 + 1]] += (Vector3)force_123.GetColumn(0);
      minus_f0 += (Vector3)force_123.GetColumn(0);
      Force[Tet[tet * 4 + 2]] += (Vector3)force_123.GetColumn(1);
      minus_f0 += (Vector3)force_123.GetColumn(1);
      Force[Tet[tet * 4 + 3]] += (Vector3)force_123.GetColumn(2);
      minus_f0 += (Vector3)force_123.GetColumn(2);
      Force[Tet[tet * 4 + 0]] += -minus_f0;
    }

    for (int i = 0; i < number; i++)
    {
      V[i] += Force[i] / mass * dt;
    }
    for (int tet = 0; tet < tet_number; ++tet)
    {
      int[] ids = new int[4];
      for (int j = 0; j < 4; ++j) { ids[j] = Tet[tet * 4 + j]; }
      Vector3 sum_v = V[ids[0]] + V[ids[1]] + V[ids[2]] + V[ids[3]];

      for (int j = 0; j < 4; ++j)
      {
        V_sum[ids[j]] += sum_v;
        V_num[ids[j]] += 4;
      }
    }

    for(int i=0; i<number; i++)
    {
      //  Update X and V here.
      V[i] = V_sum[i] / V_num[i];
      V[i] *= damp;
      X[i] += V[i] * dt;

      //  (Particle) collision with floor.
      // impulse
      float floor_y = -3.0F;
      float dist = (X[i].y - floor_y);
      if (dist < 0)
      {
        if (V[i].y >= 0) { continue; }

        Vector3 V_n0 = kFloorNormal * V[i].y;
        Vector3 V_t0 = V[i] - V_n0;
        float a =
          Math.Max(1.0F - restitution_t * (1 + restitution) * V_n0.magnitude / V_t0.magnitude, 0);
        V_t0 *= a;
        V_n0 = -restitution * V_n0;
        V[i] = V_n0 + V_t0;
        // Debug.Log("i: after " + i + " Vn: " + V_n0.ToString("F4") + " Vt: " + V_t0.ToString("F4"));
      }
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
