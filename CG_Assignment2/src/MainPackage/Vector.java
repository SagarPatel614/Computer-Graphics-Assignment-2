package MainPackage;

public class Vector {
	/* Matrix declaration */
	Matrix R = new Matrix();
	Matrix T = new Matrix();
	Matrix transR = new Matrix();
	Matrix invT = new Matrix();
	static Matrix MCW = new Matrix();
	static Matrix MWC = new Matrix();
	
	
	/* Vector Coordinates */
	float x, y, z;
	
	/* Unit Vectors */
	Vector u, n, v;
	
	/* Empty Constructor */
	public Vector() {
		//empty
	}
	
	/* Basic Constructor */
	public Vector(float x, float y, float z) {
		this.x = x;
		this.y = y;
		this.z = z;
	}
	
	public float getX() {
		return this.x;
	}
	
	public float getY() {
		return this.y;
	}
	
	public float getZ() {
		return this.z;
	}
	
	public void setX(float x) {
		this.x = x;
	}
	
	public void setY(float y) {
		this.y = y;
	}
	
	public void setZ(float z) {
		this.z = z;
	}
	
	

	static float focal = 0.05f;	/* focal length simulating 50 mm lens */
	
	/* Method for transforming from camera coordinates to world coordinates */
	public Vector transformC2W() {
		Vector V1;
		
		float x = MCW.value[0][0]*this.getX() + MCW.value[0][1]*this.getY() + MCW.value[0][2]*this.getZ() + MCW.value[0][3]*1;
		float y = MCW.value[1][0]*this.getX() + MCW.value[1][1]*this.getY() + MCW.value[1][2]*this.getZ() + MCW.value[1][3]*1;
		float z = MCW.value[2][0]*this.getX() + MCW.value[2][1]*this.getY() + MCW.value[2][2]*this.getZ() + MCW.value[2][3]*1;
		float w = MCW.value[3][0]*this.getX() + MCW.value[3][1]*this.getY() + MCW.value[3][2]*this.getZ() + MCW.value[3][3]*1;
		/*
		//System.out.println("\n");
		//MCW.print();
		float x = MCW.value[0][0]*0 + MCW.value[0][1]*0 + MCW.value[0][2]*0 + MCW.value[0][3]*1;
		float y = MCW.value[1][0]*0 + MCW.value[1][1]*0 + MCW.value[1][2]*0 + MCW.value[1][3]*1;
		float z = MCW.value[2][0]*0 + MCW.value[2][1]*0 + MCW.value[2][2]*0 + MCW.value[2][3]*1;
		float w = MCW.value[3][0]*0 + MCW.value[3][1]*0 + MCW.value[3][2]*0 + MCW.value[3][3]*1;
		*/
		if(w!=1 && w!=0) {
			x = x/w;
			y = y/w;
			z = z/w;
		}
		//System.out.println("\n\n value of w is:" + w + "\n" + "The value of x,y,z for provided vector are: " + MCW.value[0][3] +" "+this.getY()+" "+this.getZ());
		return V1 = new Vector(x,y,z);
	}
	
	/* Method for transforming from world coordinates to camera coordinates */
	public Vector transformW2C() {
		Vector V;
		float x = MWC.value[0][0]*this.getX() + MWC.value[0][1]*this.getY() + MWC.value[0][2]*this.getZ();
		float y = MWC.value[1][0]*this.getX() + MWC.value[1][1]*this.getY() + MWC.value[1][2]*this.getZ();
		float z = MWC.value[2][0]*this.getX() + MWC.value[2][1]*this.getY() + MWC.value[2][2]*this.getZ();
		float w = MWC.value[0][3]*this.getX() + MWC.value[1][3]*this.getY() + MWC.value[2][3]*this.getZ();
		if(w!=1 && w!=0) {
			x = x/w;
			y = y/w;
			z = z/w;
		}
		return V = new Vector(x,y,z);
	}
	
	/* Method for initializing Matrix for C2W transformation */
	public void initializeMCW(Vector VPN, Vector VUP, Vector VRP) {
		calculateUnitVectors(VPN, VUP);
		generateMatrixR();
		generateMatrixT(VRP);
		transposeR();
		inverseT();
		calculateMCW();
		calculateMWC();
	}

	private void calculateUnitVectors(Vector VPN_, Vector VUP_) {
		// Calculate n
		n = new Vector().calculateN(VPN_);
		
		
		// Calculate u
		u = VUP_.calculateU(VPN_);
		//u = VUP_.calculateU(n);
		
		
		// Calculate v
		v = new Vector().calculateV(n, u);
		
	}
	
	private Vector calculateN(Vector v1_) {
		float magnitude = v1_.module();
		this.x = (v1_.x / magnitude);
		this.y = (v1_.y / magnitude);
		this.z = (v1_.z / magnitude);
		
		return this;
	}
	
	private Vector calculateU(Vector v2_) {
		Vector v3 = this.crossProduct(v2_);
		float mag = (float) v3.module();
		//System.out.println(mag);
		Vector u = new Vector(v3.getX()/mag, v3.getY()/mag, v3.getZ()/mag);
		return u;
	}
	
	private Vector calculateV(Vector n, Vector u) {
		Vector v = n.crossProduct(u);
		return v;
	}
	
	float dotProduct(Vector v1_) {
		return ((this.x * v1_.x) + (this.y * v1_.y) + (this.z * v1_.z));
	}
	
	Vector crossProduct(Vector v2_) {
		Vector v3 = new Vector((this.getY()*v2_.getZ())-(this.getZ()*v2_.getY()) , (this.getX()*v2_.getZ())-(this.getZ()*v2_.getX()) , (this.getX()*v2_.getY())-(this.getY()*v2_.getX()));
		return v3;
	}

	private float module() {
		float mod = (float) Math.sqrt(Math.pow(this.x, 2) + Math.pow(this.y, 2) + Math.pow(this.z, 2));
		return mod;
	}
	
	void print() {
		System.out.println("("+ this.getX() + ","+ this.getY() +","+ this.getZ() + ")");
	}
	
	private void generateMatrixR() {
		R.MatrixR(u.getX(), u.getY(), u.getZ(), v.getX(), v.getY(), v.getZ(), n.getX(), n.getY(), n.getZ());
		
	}
	
	private void generateMatrixT(Vector VRP) {
		T.MatrixT(VRP.getX(),VRP.getY(),VRP.getZ());
		
	}
	
	private void transposeR() {
		transR.transposeR(R);
		
	}
	
	private void inverseT() {
		invT.inverseT(T);
		
	}
	
	private void calculateMCW() {
		MCW = Matrix.calculateMCW(invT, transR);
		
	}
	
	private void calculateMWC() {
		MWC = Matrix.calculateMWC(R, T);
		
	}
	
	public void copy(Vector v) {
		this.setX(v.getX());
		this.setY(v.getY());
		this.setZ(v.getZ());
	}
	
	public Vector createVector(Vector v1, Vector v2) {
		Vector nVector = new Vector();
		nVector.setX(v1.getX() - v2.getX());
		nVector.setY(v1.getY() - v2.getY());
		nVector.setZ(v1.getZ() - v2.getZ());
		//return new Vector((v1.getX() - v2.getX()) , (v1.getY() - v2.getY()) , (v1.getZ() - v2.getZ()));
		return nVector;
	}
	
	public void normalize() {
		float magnitude = this.module();
		this.x = (float) (this.x / magnitude);
		this.y = (float) (this.y / magnitude);
		this.z = (float) (this.z / magnitude);
	}
	
	public Vector multiply(float t) {
		this.x *= t;
		this.y *= t;
		this.z *= t;
		return this;
	}
	
	public Vector addVector(Vector v) {
		Vector result = new Vector();
		result.x = this.x + v.x;
		result.y = this.y + v.y;
		result.z = this.z + v.z;
		return result;
	}
	
	public Vector create2DVector(Vector v1, Vector v2) {
		Vector nVector = new Vector();
		nVector.setX(v1.getX() - v2.getX());
		nVector.setY(v1.getY() - v2.getY());
		nVector.setZ(v1.getZ() - v2.getZ());
		return nVector;
	}
	public Vector convert2D(int index) {
		if(index == 0) {
			Vector v2D = new Vector(0.0f, this.getY(), this.getZ());
			return v2D;
		} else if (index == 1) {
			Vector v2D = new Vector(this.getX(), 0.0f, this.getZ());
			return v2D;
		} else {
			Vector v2D = new Vector(this.getX(), this.getY(), 0.0f);
			return v2D;
		}
		
	}

	public Vector subtractVector(Vector newPoint) {
		Vector result = new Vector();
		result.x = this.x - newPoint.x;
		result.y = this.y - newPoint.y;
		result.z = this.z - newPoint.z;
		return result;
	}

	/*
	/* Transformation from the world to the camera coordinates 
	double Mwc[][] =
		{1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};

	/* Transformation from the camera to the world coordinates 
	double Mcw[][] = {1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 1.0};
	*/

}
