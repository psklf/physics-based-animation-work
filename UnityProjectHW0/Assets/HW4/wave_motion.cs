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

  Vector3 	cube_v = Vector3.zero;
  Vector3 	cube_w = Vector3.zero;


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

    //Step 2: Block->Water coupling
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
    float cube_size = 0.5F;
    for (int i = 0; i < size; i++)
    {
      for (int j = 0; j < size; j++)
      {
        Vector3 cur_p = new Vector3(i * 10.0F / size - 5.0F, 0.0F, j * 10.0F / size - 5.0F);
        if (Mathf.Abs(cur_p.x - c1.x) < 0.5F && Mathf.Abs(cur_p.z - c1.z) < 0.5F)
        {
          low_h[i, j] -= 0.1F;
          cg_mask[i, j] = true;
          b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
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
    Conjugate_Gradient(cg_mask, b, x1, 0, size - 1, 0, size - 1);

    //TODO: for block 2, calculate low_h.
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
        Vector3 cur_p = new Vector3(i * 10.0F / size - 5.0F, 0.0F, j * 10.0F / size - 5.0F);
        if (Mathf.Abs(cur_p.x - c2.x) < 0.5F && Mathf.Abs(cur_p.z - c2.z) < 0.5F)
        {
          low_h[i, j] -= 0.1F;
          cg_mask[i, j] = true;
          b[i, j] = (new_h[i, j] - low_h[i, j]) / rate;
        }
        else
        {
          cg_mask[i, j] = false;
          b[i, j] = 0.0F;
        }
      }
    }
    float[,] x2 = new float[size, size];
    Conjugate_Gradient(cg_mask, b, x2, 0, size - 1, 0, size - 1);

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
}
