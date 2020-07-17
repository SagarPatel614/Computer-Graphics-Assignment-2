package MainPackage;

public class Sphere {
	
	float x,y,z;	/* Center of the sphere */
	float radius;	/* Radius of the sphere */
	float Kd;		/* Diffuse Reflection Coefficient */
	
	public Sphere(float x, float y, float z, float radius, float Kd) {
		this.x = x;
		this.y = y;
		this.z = z;
		this.radius = radius;
		this.Kd = Kd;
	}

}
