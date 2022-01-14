using UnityEngine;
using System.Collections;

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
  float mass          = 300.0F;
	float angular_decay = 0.98F;
	float restitution   = 0.5f;					// for collision
	float restitution_t = 0.3f;					// for collision

	Matrix4x4 Get_Cross_Matrix(Vector3 a)
	{
		//Get the cross product matrix of vector a
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

  Quaternion Update_Rotation_w(Quaternion q, Vector3 w, float dt)
  {
    Vector3 wdt = w * dt * 0.5F;
    Quaternion omega_q = new Quaternion(wdt.x, wdt.y, wdt.z, 0);
    Quaternion dq = omega_q * q;
    return Quaternion_Add(q, dq).normalized;
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

    // then set up b and cg_mask for conjugate gradient.
    GameObject cube = GameObject.Find("Cube");
    Vector3 c1 = cube.transform.position;

    Collider collider1 = cube.GetComponent<Collider>();

    float cube_size = 0.5F;

    int li = size - 1;
    int ui = 0;
    int lj = size - 1;
    int uj = 0;
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        // Whether this cube is collide w/ current water grid
        Vector3 cur_p = new Vector3(i * 10.0F / size - 5.0F, new_h[i, j],
            j * 10.0F / size - 5.0F);

        Vector3 ray_s = cur_p;
        ray_s.y -= 1.0F;
        Ray ray = new Ray(ray_s, new Vector3(0, 1, 0));
        float distance;

        if (Mathf.Abs(cur_p.x - c1.x) < 0.5F && Mathf.Abs(cur_p.z - c1.z) < 0.5F &&
            Mathf.Abs(cur_p.y - c1.y) < 0.5F)
        {
          /*
          if (collider1.bounds.IntersectRay(ray, out distance))
          {
            Debug.Log("Raycast " + distance);
          }*/
          low_h[i, j] = c1.y - 0.5F;
          cg_mask[i, j] = true;
          b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
          // mark bound
          if (i < li) { li = i; }
          if (i > ui) { ui = i; }
          if (j < lj) { lj = j; }
          if (j > uj) { uj = j; }
        }
        else
        {
          cg_mask[i, j] = false;
          b[i, j] = 0.0F;
        }
      }
    }

    // Mesh cube_mesh = cube.GetComponent<MeshFilter>().mesh;

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
    // (bool[,] mask, float[,] b, float[,] x, int li, int ui, int lj, int uj)
    Conjugate_Gradient(cg_mask, b, x1, li, ui, lj, uj);

    //TODO: for block 2, calculate low_h.

    li = size - 1;
    ui = 0;
    lj = size - 1;
    uj = 0;

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
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        Vector3 cur_p =
          new Vector3(i * 10.0F / size - 5.0F, new_h[i, j], j * 10.0F / size - 5.0F);
        if (Mathf.Abs(cur_p.x - c2.x) < 0.5F && Mathf.Abs(cur_p.z - c2.z) < 0.5F &&
            Mathf.Abs(cur_p.y - c2.y) < 0.5F)
        {
          low_h[i, j] = c2.y - 0.5F;
          cg_mask[i, j] = true;
          b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
          // mark bound
          if (i < li) { li = i; }
          if (i > ui) { ui = i; }
          if (j < lj) { lj = j; }
          if (j > uj) { uj = j; }
        }
        else
        {
          cg_mask[i, j] = false;
          b[i, j] = 0.0F;
        }
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

    //Step 4: Water->Block coupling.
    //More TODO here.

    // deltaA
    float area = 10.0F * 10.0F / (size * size);
    Vector3 water_f1 = Vector3.zero;
    Vector3 water_f2 = Vector3.zero;
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        Vector3 cur_p =
          new Vector3(i * 10.0F / size - 5.0F, old_h[i, j], j * 10.0F / size - 5.0F);
        if (Mathf.Abs(cur_p.x - c1.x) < 0.5F && Mathf.Abs(cur_p.z - c1.z) < 0.5F &&
            Mathf.Abs(cur_p.y - c1.y) < 0.5F)
        {
          // water_f1.y += (old_h[i, j] - h[i, j]) * 9.8F * 997.0F * area;
          // water_f1.y += (vh[i, j]) * rate / dt / dt;
          water_f1.y += (vh[i, j]) * 9.8F * 997.0F * area;
        }
        if (Mathf.Abs(cur_p.x - c2.x) < 0.5F && Mathf.Abs(cur_p.z - c2.z) < 0.5F &&
            Mathf.Abs(cur_p.y - c2.y) < 0.5F)
        {
          water_f2.y += (vh[i, j]) * 9.8F * 997.0F * area;
        }
      }
    }

    Debug.Log("water force 1" + water_f1);
    cube_v[0] += (kGravity + water_f1 / mass) * dt;
    cube_v[0] *= damping;


    // Calculate w
    cube_w[0] *= angular_decay;

    // Update position and rotation

		Quaternion q1 = cube.transform.rotation;
    q1 = Update_Rotation_w(q1, cube_w[0], dt);

    c1 += cube_v[0] * dt;
    Debug.Log("pos 1" + c1);

    cube.transform.position = c1;
    cube.transform.rotation = q1;

    cube_v[1] += (kGravity + water_f2 / mass) * dt;
    cube_v[1] *= damping;
    c2 += cube_v[1] * dt;
    cube2.transform.position = c2;
  }

  // Update is called once per frame
  void Update ()
  {
    Mesh mesh = GetComponent<MeshFilter> ().mesh;
    Vector3[] X    = mesh.vertices;
    float[,] new_h = new float[size, size];
    float[,] h     = new float[size, size];

    // : Load X.y into h.
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

    Debug.Log(" start Shallow_Wave" );
    for(int l=0; l<8; l++)
    {
      Shallow_Wave(old_h, h, new_h);
    }
    Debug.Log(" stop Shallow_Wave" );

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
}
