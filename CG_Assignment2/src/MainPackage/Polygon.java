package MainPackage;

public class Polygon {
	
	Vector v[] = new Vector[4];	/* List of all the vertices */
	float N[] = new float[3];		/* Normal to the Polygon */
	float Kd;						/* Diffuse Reflection Coefficient */
	
	public Polygon(Vector v0, Vector v1, Vector v2, Vector v3, float Kd) {
		this.v[0] = v0;
		this.v[1] = v1;
		this.v[2] = v2;
		this.v[3] = v3;
		//this.N = N;
		this.Kd = Kd;
	}

}
