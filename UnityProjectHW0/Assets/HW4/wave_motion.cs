using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class wave_motion : MonoBehaviour
{
  int size 		= 100;
  float rate 		= 0.005f;
  float gamma		= 0.004f;
  float damping 	= 0.996f;
  float[,] 	old_h;
  float[,]	low_h;
  float[,]	vh;
  float[,]	b;

  bool [,]	cg_mask;
  float[,]	cg_p;
  float[,]	cg_r;
  float[,]	cg_Ap;
  bool 	tag = true;

  Vector3[] cube_v = new Vector3[2];
  Vector3[] cube_w = new Vector3[2];

  Vector3 kGravity    = new Vector3(0, -9.8F, 0);
  float dt            = 0.003F;
  float mass          = 500.0F;
  float m;
  float angular_decay = 0.8F;
  Matrix4x4 I_ref = Matrix4x4.identity;

  // Get the cross product matrix of vector a
	Matrix4x4 Get_Cross_Matrix(Vector3 a)
	{
		Matrix4x4 A = Matrix4x4.zero;
		A [0, 0] = 0;
		A [0, 1] = -a [2];
		A [0, 2] = a [1];
		A [1, 0] = a [2];
		A [1, 1] = 0;
		A [1, 2] = -a [0];
		A [2, 0] = -a [1];
		A [2, 1] = a [0];
		A [2, 2] = 0;
		A [3, 3] = 1;
		return A;
	}

  Quaternion Quaternion_Add(Quaternion l, Quaternion r)
  {
    return new Quaternion(l.x + r.x, l.y + r.y,
        l.z + r.z, l.w + r.w);
  }

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

  // Use this for initialization
  void Start ()
  {
    Mesh mesh = GetComponent<MeshFilter> ().mesh;
    mesh.Clear ();

    Vector3[] X=new Vector3[size*size];

    for (int i=0; i<size; i++)
      for (int j=0; j<size; j++)
      {
        X[i*size+j].x=i*0.1f-size*0.05f;
        X[i*size+j].y=0;
        X[i*size+j].z=j*0.1f-size*0.05f;
      }

    int[] T = new int[(size - 1) * (size - 1) * 6];
    int index = 0;
    for (int i=0; i<size-1; i++)
      for (int j=0; j<size-1; j++)
      {
        T[index*6+0]=(i+0)*size+(j+0);
        T[index*6+1]=(i+0)*size+(j+1);
        T[index*6+2]=(i+1)*size+(j+1);
        T[index*6+3]=(i+0)*size+(j+0);
        T[index*6+4]=(i+1)*size+(j+1);
        T[index*6+5]=(i+1)*size+(j+0);
        index++;
      }
    mesh.vertices  = X;
    mesh.triangles = T;
    mesh.RecalculateNormals ();

    low_h 	= new float[size,size];
    old_h 	= new float[size,size];
    vh 	  	= new float[size,size];
    b 	  	= new float[size,size];

    cg_mask	= new bool [size,size];
    cg_p 	= new float[size,size];
    cg_r 	= new float[size,size];
    cg_Ap 	= new float[size,size];

    for (int i=0; i<size; i++)
      for (int j=0; j<size; j++)
      {
        low_h[i,j]=99999;
        old_h[i,j]=0;
        vh[i,j]=0;
      }

    cube_v[0] = Vector3.zero;
    cube_v[1] = Vector3.zero;
    cube_w[0] = Vector3.zero;
    cube_w[1] = Vector3.zero;

    GameObject cube1 = GameObject.Find("Cube");
    GameObject cube2 = GameObject.Find("Block");
    Mesh cube_mesh = cube1.GetComponent<MeshFilter>().mesh;
    m = mass / cube_mesh.vertices.Length;
    InitCubeIref(cube_mesh.vertices);
  }

  void A_Times(bool[,] mask, float[,] x, float[,] Ax, int li, int ui, int lj, int uj)
  {
    for(int i=li; i<=ui; i++)
      for(int j=lj; j<=uj; j++)
        if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
        {
          Ax[i,j]=0;
          if(i!=0)		Ax[i,j]-=x[i-1,j]-x[i,j];
          if(i!=size-1)	Ax[i,j]-=x[i+1,j]-x[i,j];
          if(j!=0)		Ax[i,j]-=x[i,j-1]-x[i,j];
          if(j!=size-1)	Ax[i,j]-=x[i,j+1]-x[i,j];
        }
  }

  float Dot(bool[,] mask, float[,] x, float[,] y, int li, int ui, int lj, int uj)
  {
    float ret=0;
    for(int i=li; i<=ui; i++)
      for(int j=lj; j<=uj; j++)
        if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
        {
          ret+=x[i,j]*y[i,j];
        }
    return ret;
  }

  void Conjugate_Gradient(bool[,] mask, float[,] b, float[,] x, int li, int ui, int lj, int uj)
  {
    //Solve the Laplacian problem by CG.
    A_Times(mask, x, cg_r, li, ui, lj, uj);

    for(int i=li; i<=ui; i++)
      for(int j=lj; j<=uj; j++)
        if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
        {
          cg_p[i,j]=cg_r[i,j]=b[i,j]-cg_r[i,j];
        }

    float rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);

    for(int k=0; k<128; k++)
    {
      if(rk_norm<1e-10f)	break;
      A_Times(mask, cg_p, cg_Ap, li, ui, lj, uj);
      float alpha=rk_norm/Dot(mask, cg_p, cg_Ap, li, ui, lj, uj);

      for(int i=li; i<=ui; i++)
        for(int j=lj; j<=uj; j++)
          if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
          {
            x[i,j]   +=alpha*cg_p[i,j];
            cg_r[i,j]-=alpha*cg_Ap[i,j];
          }

      float _rk_norm=Dot(mask, cg_r, cg_r, li, ui, lj, uj);
      float beta=_rk_norm/rk_norm;
      rk_norm=_rk_norm;

      for(int i=li; i<=ui; i++)
        for(int j=lj; j<=uj; j++)
          if(i>=0 && j>=0 && i<size && j<size && mask[i,j])
          {
            cg_p[i,j]=cg_r[i,j]+beta*cg_p[i,j];
          }
    }
  }

  void Shallow_Wave(float[,] old_h, float[,] h, float [,] new_h)
  {
    // Step 1:
    // Compute new_h based on the shallow wave model.
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        new_h[i, j] = h[i, j] + damping * (h[i, j] - old_h[i, j]);
      }
    }
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        float h_ij = h[i, j];
        if (i > 0) { new_h[i, j] += rate * (h[i - 1, j] - h_ij); }
        if (i < (size - 1)) { new_h[i, j] += rate * (h[i + 1, j] - h_ij); }
        if (j > 0) { new_h[i, j] += rate * (h[i, j - 1] - h_ij); }
        if (j < (size - 1)) { new_h[i, j] += rate * (h[i, j + 1] - h_ij); }
      }
    }

    // Step 2: Block->Water coupling

    // for block 1, calculate low_h.
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        low_h[i, j] = new_h[i, j];
      }
    }

    // Cube related constants

    float cube_size = 0.5F;
    int cube_step = 9;

    // then set up b and cg_mask for conjugate gradient.
    GameObject cube = GameObject.Find("Cube");
    Vector3 c1 = cube.transform.position;

    int li1 = (int) ((c1.x - cube_size + 5.0F) * 0.1F * size);
    int ui1 = li1 + cube_step;
    int lj1 = (int) ((c1.z - cube_size + 5.0F) * 0.1F * size);
    int uj1 = lj1 + cube_step;

    Collider collider1 = cube.GetComponent<Collider>();

    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        cg_mask[i, j] = false;
        b[i, j] = 0.0F;
      }
    }  // end of water grid loop

    List<Vector3> hit_list_1 = new List<Vector3>();
    List<int[]> index_list_1 = new List<int[]>();
    for (int i = li1; i <= ui1; ++i)
    {
      for (int j = lj1; j <= uj1; ++j)
      {
        if (i < 0 || i >= size || j < 0 || j >= size) { continue; }
        // Whether this cube is collide w/ current water grid
        Vector3 cur_p = new Vector3(i * 10.0F / size - 5.0F, new_h[i, j],
            j * 10.0F / size - 5.0F);

        Vector3 ray_s = cur_p;
        ray_s.y -= 10.0F;
        Ray ray = new Ray(ray_s, new Vector3(0, 1, 0));

        float distance;
        if (collider1.bounds.IntersectRay(ray, out distance))
        {
          Vector3 hit_point = ray_s + new Vector3(0, distance, 0);

          if (hit_point.y < new_h[i, j])
          {
             // low_h[i, j] = c1.y - 0.5F;
            low_h[i, j] = hit_point.y;
            cg_mask[i, j] = true;
            b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;

            // save point
            hit_list_1.Add(hit_point - c1);
            index_list_1.Add(new int[2]{i, j});
          }
        }  // end of intersect ray
      }
    }  // end of cube1 low h loop


    // Solve the Poisson equation to obtain vh (virtual height).
    //     save result to x1
    float[,] x1 = new float[size, size];
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        x1[i, j] = 0;
      }
    }
    Conjugate_Gradient(cg_mask, b, x1, li1, ui1, lj1, uj1);

    //TODO: for block 2, calculate low_h.
    //
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        low_h[i, j] = new_h[i, j];
      }
    }

    //TODO: then set up b and cg_mask for conjugate gradient.
    //TODO: Solve the Poisson equation to obtain vh (virtual height).
    GameObject cube2 = GameObject.Find("Block");
    Vector3 c2 = cube2.transform.position;

    int li = (int) ((c2.x - 0.5F + 5.0F) * 0.1F * size);
    int ui = li + cube_step;
    int lj = (int) ((c2.z - 0.5F + 5.0F) * 0.1F * size);
    int uj = lj + cube_step;

    Collider collider2 = cube2.GetComponent<Collider>();

    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        cg_mask[i, j] = false;
        b[i, j] = 0.0F;
      }
    }

    List<Vector3> hit_list_2 = new List<Vector3>();
    List<int[]> index_list_2 = new List<int[]>();
    for (int i = li; i <= ui; ++i)
    {
      for (int j = lj; j <= uj; ++j)
      {
        if (i < 0 || i >= size || j < 0 || j >= size) {continue;}
        // Whether this cube is collide w/ current water grid
        Vector3 cur_p = new Vector3(i * 10.0F / size - 5.0F, new_h[i, j],
            j * 10.0F / size - 5.0F);

        Vector3 ray_s = cur_p;
        ray_s.y -= 1.0F;
        Ray ray = new Ray(ray_s, new Vector3(0, 1, 0));

        float distance;
        if (collider2.bounds.IntersectRay(ray, out distance))
        {
          Vector3 hit_point = ray_s + new Vector3(0, distance, 0);
          if (hit_point.y < new_h[i, j])
          {
            //low_h[i, j] = c2.y - 0.5F;
            low_h[i, j] = hit_point.y;
            cg_mask[i, j] = true;
            b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
            hit_list_2.Add(hit_point - c2);
            index_list_2.Add(new int[2]{i, j});
          }
        }  // end of intersect ray
      }
    }
    float[,] x2 = new float[size, size];
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        x2[i, j] = 0;
      }
    }
    Conjugate_Gradient(cg_mask, b, x2, li, ui, lj, uj);

    // Diminish vh by x1+x2
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        vh[i, j] = (x1[i, j] + x2[i, j]) * gamma;
      }
    }

    // Update new_h by vh.
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        float v_ij = vh[i, j];
        if (i != 0) { new_h[i, j] += rate * (vh[i - 1, j] - v_ij); }
        if (i != (size - 1)) { new_h[i, j] += rate * (vh[i + 1, j] - v_ij); }
        if (j != 0) { new_h[i, j] += rate * (vh[i, j - 1] - v_ij); }
        if (j != (size - 1)) { new_h[i, j] += rate * (vh[i, j + 1] - v_ij); }
      }
    }

    // Step 3
    // old_h <- h; h <- new_h;
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        old_h[i, j] = h[i, j];
        h[i, j] = new_h[i, j];
      }
    }

    // Step 4: Water->Block coupling.
    // More TODO here.

    // Dynamics for cubes

    // deltaA
    float area = 10.0F * 10.0F / (size * size);

    // cube1
    Quaternion q1 = cube.transform.rotation;
    Matrix4x4 R1 = Matrix4x4.Rotate(q1);

    Vector3 water_f1 = Vector3.zero;
    Vector3 water_f2 = Vector3.zero;

    // cube1
    // Get rotation by impuse
    Vector3 sum_tao = Vector3.zero;
    for (int i = 0; i < hit_list_1.Count; ++i)
    {
      Vector3 force = new Vector3(0,
          vh[index_list_1[i][0], index_list_1[i][1]] * 9.8F * 997.0F * area,
          0);
      // don't why but should give a scale
      sum_tao += Vector3.Cross(hit_list_1[i], force * 10);
  //    Debug.Log("hit force" + force);
 //     Debug.Log("point" + hit_list_1[i]);
      water_f1 += force;
    }
 //   Debug.Log("sum tao" + sum_tao.ToString("F4"));
    Matrix4x4 I = R1 * I_ref * R1.transpose;

    cube_w[0] += (I.inverse.MultiplyPoint3x4(sum_tao) * dt);
    cube_w[0] *= angular_decay;

    cube_v[0] += (kGravity + water_f1 / mass) * dt;
    cube_v[0] *= damping;

    // cube2
    Quaternion q2 = cube2.transform.rotation;
    Matrix4x4 R2 = Matrix4x4.Rotate(q2);

    Vector3 sum_tao_2 = Vector3.zero;
    for (int i = 0; i < hit_list_2.Count; ++i)
    {
      Vector3 force = new Vector3(0,
          vh[index_list_2[i][0], index_list_2[i][1]] * 9.8F * 997.0F * area,
          0);
      sum_tao_2 += Vector3.Cross(hit_list_2[i], force * 1);
      water_f2 += force;
    }

    Matrix4x4 I2 = R2 * I_ref * R2.transpose;
    cube_w[1] += I2.inverse.MultiplyPoint3x4(sum_tao_2) * dt;
    cube_w[1] *= angular_decay;

    cube_v[1] += (kGravity + water_f2 / mass) * dt;
    cube_v[1] *= damping;

    // Update position and rotation
    q1 = Update_Rotation_w(q1, cube_w[0], dt);
    c1 += cube_v[0] * dt;
    cube.transform.position = c1;
    cube.transform.rotation = q1;

    q2 = Update_Rotation_w(q2, cube_w[1], dt);
    c2 += cube_v[1] * dt;
    cube2.transform.position = c2;
    cube2.transform.rotation = q2;
  }

  // Update is called once per frame
  void Update ()
  {
    Mesh mesh = GetComponent<MeshFilter> ().mesh;
    Vector3[] X    = mesh.vertices;
    float[,] new_h = new float[size, size];
    float[,] h     = new float[size, size];

    //  Load X.y into h.
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        h[i, j] = X[i * size + j].y;
      }
    }

    // Add random water.
    if (Input.GetKeyDown ("r"))
    {
      float random_water = Random.value - 0.5F;
      int col = (int)(Random.value * (size - 3) + 1);
      int row = (int)(Random.value * (size - 3) + 1);
      h[col, row] += random_water;
      h[col - 1, row] -= random_water * 0.25F;
      h[col + 1, row] -= random_water * 0.25F;
      h[col, row - 1] -= random_water * 0.25F;
      h[col, row + 1] -= random_water * 0.25F;
    }

    for(int l=0; l<8; l++)
    {
      Shallow_Wave(old_h, h, new_h);
    }

    // Store h back into X.y and recalculate normal.
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        X[i * size + j].y = h[i, j];
      }
    }

    mesh.vertices = X;
    mesh.RecalculateNormals();
  }

  Quaternion Update_Rotation_w(Quaternion q, Vector3 w, float dt)
  {
    Vector3 wdt = w * dt * 0.5F;
    Quaternion omega_q = new Quaternion(wdt.x, wdt.y, wdt.z, 0);
    Quaternion dq = omega_q * q;
    return Quaternion_Add(q, dq).normalized;
  }

  Matrix4x4 InitCubeIref(Vector3[] vertices)
  {
    float m = mass / vertices.Length;
    for (int i=0; i<vertices.Length; i++)
    {
      float diag=m*vertices[i].sqrMagnitude;
      I_ref[0, 0]+=diag;
      I_ref[1, 1]+=diag;
      I_ref[2, 2]+=diag;
      I_ref[0, 0]-=m*vertices[i][0]*vertices[i][0];
      I_ref[0, 1]-=m*vertices[i][0]*vertices[i][1];
      I_ref[0, 2]-=m*vertices[i][0]*vertices[i][2];
      I_ref[1, 0]-=m*vertices[i][1]*vertices[i][0];
      I_ref[1, 1]-=m*vertices[i][1]*vertices[i][1];
      I_ref[1, 2]-=m*vertices[i][1]*vertices[i][2];
      I_ref[2, 0]-=m*vertices[i][2]*vertices[i][0];
      I_ref[2, 1]-=m*vertices[i][2]*vertices[i][1];
      I_ref[2, 2]-=m*vertices[i][2]*vertices[i][2];
    }
    I_ref [3, 3] = 1;
    return I_ref;
  }

  bool IsPointofCube(Vector3 center, Vector3 p, float l)
  {
    return (Mathf.Abs(p.x - center.x) <= l && Mathf.Abs(p.y - center.y) <= l &&
        Mathf.Abs(p.z - center.z) <= l);
  }
}
