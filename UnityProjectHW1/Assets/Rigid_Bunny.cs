using UnityEngine;
using System.Collections;
using System.Collections.Generic;

public class Rigid_Bunny : MonoBehaviour
{
	bool launched 		= false;
	float dt 			= 0.015f;
	Vector3 v 			= new Vector3(0, 0, 0);	// velocity
	Vector3 w 			= new Vector3(0, 0, 0);	// angular velocity

	float mass;									// mass
	Matrix4x4 I_ref;							// reference inertia

	float linear_decay	= 0.999f;				// for velocity decay
	float angular_decay	= 0.98f;
	float restitution 	= 0.3f;					// for collision

  Vector3 gravity_a = new Vector3(0, -9.8F, 0);

	// Use this for initialization
	void Start ()
	{
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;

		float m=1;
		mass=0;
		for (int i=0; i<vertices.Length; i++)
		{
			mass += m;
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
	}

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

  float Phi(Vector3 point, Vector3 P, Vector3 N)
  {
    return Vector3.Dot((point - P), N);
  }

  // In this function, update v and w by the impulse due to the collision with
	//a plane <P, N>
	void Collision_Impulse(Vector3 P, Vector3 N)
	{
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;

    Matrix4x4 R = Matrix4x4.Rotate(transform.rotation);

    // Get all collide vertices
		List<Vector3> collide_vertices = new List<Vector3>();
		for (int i = 0; i < vertices.Length; i++)
		{
      Vector3 x_i = transform.position + R.MultiplyPoint3x4(vertices[i]);
      if (Phi(x_i, P, N) < 0)
      {
        collide_vertices.Add(vertices[i]);
      }
    }
    if (collide_vertices.Count == 0) { return; }

    Vector3 normal_v = Vector3.zero;
    for (int i = 0 ; i < collide_vertices.Count; ++i)
    {
      normal_v += collide_vertices[i];
    }
    normal_v /= collide_vertices.Count;
    Debug.Log("Get normal vertex: " + normal_v);

    Vector3 Rxr_i = R.MultiplyPoint3x4(normal_v);
    Vector3 sumv = v + Vector3.Cross(w, Rxr_i);
    if (Vector3.Dot(sumv, N) < 0)
    {
      Vector3 v_n = N * Vector3.Dot(sumv, N);
      Vector3 v_t = sumv - v_n;
      float a = Mathf.Max(
          1 - 0.5F * (1 + restitution) * v_n.magnitude / v_t.magnitude,
          0);
      v_n = -v_n * restitution;
      v_t = v_t * a;
      Vector3 sum_v_new = v_n + v_t;

      Matrix4x4 I = R * I_ref * R.transpose;
      Matrix4x4 Rri_cross_mat = Get_Cross_Matrix(Rxr_i);
      Matrix4x4 inv_mass_mat = Matrix4x4.identity;
      float inv_mass = 1.0F / mass;
      inv_mass_mat[0, 0] = inv_mass;
      inv_mass_mat[1, 1] = inv_mass;
      inv_mass_mat[2, 2] = inv_mass;
      inv_mass_mat[3, 3] = inv_mass;
      Matrix4x4 K =  Matrix4_Sub(inv_mass_mat, Rri_cross_mat * I.inverse * Rri_cross_mat);
      Vector3 j = K.inverse.MultiplyPoint3x4(sum_v_new - sumv);

      // Update real v and w
      v += j / mass;
      Vector3 add_w = I.inverse * (Rri_cross_mat.MultiplyPoint3x4(j));
      w += add_w;
      Debug.Log("The : Collision " + v + " w: " + w);
    }
  }

	// Update is called once per frame
	void Update ()
	{
		//Game Control
		if(Input.GetKey("r"))
		{
			transform.position = new Vector3 (0, 0.6f, 0);
      transform.rotation = new Quaternion(0, 0, 0, 1);
			restitution = 0.5f;
			launched=false;
		}
		if(Input.GetKey("l"))
		{
			v = new Vector3 (1, 2, 0);
      w = new Vector3 (3, 0, 3);
			launched=true;
		}

		// Part I: Update velocities
    if (launched)
    {
      // Calculate v
      v += gravity_a * mass * dt / mass;
      v *= linear_decay;

      // Calculate w
      w *= angular_decay;
    }

		// Part II: Collision Impulse
		Collision_Impulse(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0));
		// Collision_Impulse(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));

		// Part III: Update position & orientation
		//Update linear status
		Vector3 x    = transform.position;
		//Update angular status
		Quaternion q = transform.rotation;

    if (launched) {
      x += v * dt;

      Vector3 wdt = w * dt * 0.5F;
      Quaternion omega_q = new Quaternion(wdt.x, wdt.y, wdt.z, 0);
      Quaternion dq = omega_q * q;
      q = Quaternion_Add(q,dq).normalized;
      // Debug.Log("w: " + w + " q: " + q);
    }

		// Part IV: Assign to the object
		transform.position = x;
		transform.rotation = q;
	}
}
