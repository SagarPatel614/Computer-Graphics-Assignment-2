package MainPackage;

public class Matrix {
	
	float[][] value = {{0,0,0,0},{0,0,0,0},{0,0,0,0},{0,0,0,0}};

	Matrix(){
		
	}
	
	//create rotational matrix
	void MatrixR(float Ux, float Uy, float Uz, float Vx, float Vy, float Vz, float Nx, float Ny, float Nz){
		this.value[0][0] = Ux;
		this.value[0][1] = Uy;
		this.value[0][2] = Uz;
		this.value[0][3] = 0;
		
		this.value[1][0] = Vx;
		this.value[1][1] = Vy;
		this.value[1][2] = Vz;
		this.value[1][3] = 0;
		
		this.value[2][0] = Nx;
		this.value[2][1] = Ny;
		this.value[2][2] = Nz;
		this.value[2][3] = 0;
		
		this.value[3][0] = 0;
		this.value[3][1] = 0;
		this.value[3][2] = 0;
		this.value[3][3] = 1;
	}
	
	//calculate translation matrix
	void MatrixT(float VRPx, float VRPy, float VRPz){
		this.value[0][0] = 1;
		this.value[0][1] = 0;
		this.value[0][2] = 0;
		this.value[0][3] = 0-VRPx;
		
		this.value[1][0] = 0;
		this.value[1][1] = 1;
		this.value[1][2] = 0;
		this.value[1][3] = 0-VRPy;
		
		this.value[2][0] = 0;
		this.value[2][1] = 0;
		this.value[2][2] = 1;
		this.value[2][3] = 0-VRPz;
		
		this.value[3][0] = 0;
		this.value[3][1] = 0;
		this.value[3][2] = 0;
		this.value[3][3] = 1;
	}
	
	//calculate matrix multiplication
	public Matrix multiply(Matrix T) {
		Matrix result = new Matrix();
		for (int i = 0; i<4; i++)
		{
			for (int j = 0; j<4; j++)
			{
				for (int k = 0; k<4; k++)
				{
					result.value[i][j] += this.value[i][k] * T.value[k][j];
				}
			}
		}
		
		return result;
	}
	
	//print the matrix
	public void print() {
		for (int i = 0; i<4 ; i++) {
			for (int j = 0; j<4 ; j++) {
				System.out.print(this.value[i][j] + " ");
			}
			System.out.print("\n");
		}
	}
	
	//transpose matrix R
		public void transposeR(Matrix R) {
			this.value[0][0] = R.value[0][0];
			this.value[1][0] = R.value[0][1];
			this.value[2][0] = R.value[0][2];
			this.value[3][0] = R.value[0][3];
			
			this.value[0][1] = R.value[1][0];
			this.value[1][1] = R.value[1][1];
			this.value[2][1] = R.value[1][2];
			this.value[3][1] = R.value[1][3];
			
			this.value[0][2] = R.value[2][0];
			this.value[1][2] = R.value[2][1];
			this.value[2][2] = R.value[2][2];
			this.value[3][2] = R.value[2][3];
			
			this.value[0][3] = R.value[3][0];
			this.value[1][3] = R.value[3][1];
			this.value[2][3] = R.value[3][2];
			this.value[3][3] = R.value[3][3];
		}
		
		//calculate inverse of T
		public void inverseT(Matrix T) {
			this.value[0][0] = T.value[0][0];
			this.value[0][1] = T.value[0][1];
			this.value[0][2] = T.value[0][2];
			if(T.value[0][3] != 0.0) {
				this.value[0][3] = T.value[0][3]*(-1);
			} else {
				this.value[0][3] = T.value[0][3];
			}
			
			this.value[1][0] = T.value[1][0];
			this.value[1][1] = T.value[1][1];
			this.value[1][2] = T.value[1][2];
			if(T.value[1][3] != 0.0) {
				this.value[1][3] = T.value[1][3]*(-1);
			} else {
				this.value[1][3] = T.value[1][3];
			}
			
			this.value[2][0] = T.value[2][0];
			this.value[2][1] = T.value[2][1];
			this.value[2][2] = T.value[2][2];
			if(T.value[2][3] != 0.0) {
				this.value[2][3] = T.value[2][3]*(-1);
			} else {
				this.value[2][3] = T.value[2][3];
			}
			
			this.value[3][0] = T.value[3][0];
			this.value[3][1] = T.value[3][1];
			this.value[3][2] = T.value[3][2];
			this.value[3][3] = T.value[3][3];
		}
		
		//calculate camera to world matrix
		static public Matrix calculateMCW(Matrix R, Matrix T) {
			Matrix MCW = R.multiply(T);
			return MCW;
		}
		
		//calculate world to camera matrix
		static public Matrix calculateMWC(Matrix R, Matrix T) {
			Matrix MWC = R.multiply(T);
			return MWC;
		}
		
		//calculate world to light matrix
		static public Matrix calculateMWL(Matrix R, Matrix T) {
			Matrix MWL = R.multiply(T);
			return MWL;
		}

		//calculate light to world matrix
		static public Matrix calculateMLW(Matrix R, Matrix T) {
			Matrix MLW = R.multiply(T);
			return MLW;
		}

		//calculate camera to light matrix
		static public Matrix calculateMCL(Matrix MCW, Matrix MWL) {
			//MCW*MWL
			Matrix MCL = MCW.multiply(MWL);
			return MCL;
		}

		//calculate light to camera matrix
		static public Matrix calculateMLC(Matrix MLW, Matrix MWC) {
			//MLW*MWC
			Matrix MLC = MLW.multiply(MWC);
			return MLC;
		}
		
		//check the test vector
		public Vector multiplyVector(Vector test) {
			float x = this.value[0][0]*test.getX() + this.value[0][1]*test.getY() + this.value[0][2]*test.getZ() + this.value[0][3]*1;
			float y = this.value[1][0]*test.getX() + this.value[1][1]*test.getY() + this.value[1][2]*test.getZ() + this.value[1][3]*1;
			float z = this.value[2][0]*test.getX() + this.value[2][1]*test.getY() + this.value[2][2]*test.getZ() + this.value[2][3]*1;
			float w = this.value[0][3]*test.getX() + this.value[1][3]*test.getY() + this.value[2][3]*test.getZ() + this.value[3][3]*1;
			if(w!=1 && w!=0) {
				x = x/w;
				y = y/w;
				z = z/w;
			}
			Vector v = new Vector(x,y,z);
			return v;
			
			//System.out.println("The result of Transformation is : " + x + " " + y + " " + z);
		}

}
