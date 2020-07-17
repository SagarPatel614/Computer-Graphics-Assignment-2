package MainPackage;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.awt.image.RenderedImage;
import java.io.File;
import java.io.IOException;

import javax.imageio.ImageIO;

public class MainClass {
	// Declaration of Output Image files
	private final static String filenamePolygon = "Polygon.jpeg";
    private final static String filenameSphere = "Sphere.jpeg";
    private static BufferedImage img = null;
    
    // Controller for selecting which object is in the space. Uncomment any one to obtain the result for the specified geometry
    static final String object = "Sphere";
    //static final String object = "Polygon";
    
	// Definition of Image buffer
	final static int ROWS = 512;
	final static int COLS = 512;
	static int[][] image = new int[ROWS][COLS];
	
	/* Definition of the window on the image plane in the camera coordinates.
	 * For mapping (j,i) in the screen coordinate to (x,y) on the image plane.
	 * The specified screen size stimulates the 35mm film.
	 */
	static float Xmin = 0.0175f;
	static float Ymin = -0.0175f;
	static float Xmax = -0.0175f;
	static float Ymax = 0.0175f;
	
	// Initializing Sphere and Polygon Objects
	static Sphere S = new Sphere(1.0f, 1.0f, 1.0f, 0.5f, 0.75f);	//Vector[] vertices = new Vector[] { new Vector(0.0f,0.0f,0.0f), new Vector(0.0f,0.0f,0.0f), new Vector(0.0f,0.0f,0.0f), new Vector(0.0f,0.0f,0.0f)};
	static Polygon PL = new Polygon(new Vector(0.0f,0.0f,0.0f), new Vector(0.0f,0.0f,2.0f), new Vector(2.0f,0.0f,2.0f), new Vector(2.0f,0.0f,0.0f), 0.8f);
	
	// Declaring Vector Variables
	static Vector P0 = new Vector();
	static Vector dirV = new Vector();
	static Vector P = new Vector();
	static Vector N = new Vector();
	static Vector colorVector = new Vector();
	static Vector polyNormal = new Vector();
	
	/* definition of the camera parameters */
	static Vector VRP = new Vector(1.0f, 2.0f, 3.5f);
	static Vector VPN = new Vector(0.0f, -1.0f, -2.5f);
	static Vector VUP = new Vector(0.0f, 1.0f, 0.0f);
	
	static float Kd;
	
	static Matrix polyTranslate = new Matrix();
	
	// Setting Light Source position and intensity parameters
	static Vector LRP = new Vector(-10.0f, 10.0f, 2.0f);
	static float Ip = 200.0f;
	
	//Global variables
	static Boolean found;
	static Boolean flag = true;
	
	// Main function
	public static void main(String[] args) throws IOException {
		//local variables
		int i, j, c = 0, hits = 0;
		
		img = map(512, 512);
		
		// Initialize Required Vector Data
		Vector ini = new Vector();
		System.out.println("Vector VRP :");
		VRP.print();
		System.out.println("Vector VPN :");
		VPN.print();
		System.out.println("Vector VUP :");
		VUP.print();
		System.out.println("Vector LRP :");
		LRP.print();
		ini.initializeMCW(VPN, VUP, VRP);
		
		for(i = 0; i < ROWS; i++) {
			for(j = 0; j < COLS; j++) {
				//construct the ray, V, starting from the CenterOfProjection, P0, and passing through the pixel (i,j)
				ray_construction(i,j);
				
				// if the ray intersects with an object save the returned shading value into the image buffer
				if(ray_tracing(P0, dirV)) { 
					float C = shading(P, N, Kd);
					image[i][j] = (int) C;
					colorVector = new Vector(1.0f,0.5f,1.5f).multiply(N.dotProduct(new Vector().createVector(LRP, P)) * Kd);
					//value of r, g, b for color
					int r = (int) Math.min(255 ,255 * Math.max(0, Math.min(1, ((colorVector.getX() < 0.0f)? colorVector.getX()*-1.0f: colorVector.getX()))) + 10);
					int g = (int) Math.min(255 ,255 * Math.max(0, Math.min(1, ((colorVector.getY() < 0.0f)? colorVector.getY()*-1.0f: colorVector.getY()))) + 10);
					int b = (int) Math.min(255 ,255 * Math.max(0, Math.min(1, ((colorVector.getZ() < 0.0f)? colorVector.getZ()*-1.0f: colorVector.getZ()))) + 10);
					Color color = new Color(r,g,b);
					img.setRGB(j, i, color.getRGB());
					hits++;
				}
				else { 
					//for no hits
					image[i][j] = (char) 0;
				}
			}
		}
		
		System.out.println("\n\n Image Created!!!!\nNumber of Hits are: " + hits);
		if(object.equals("Sphere")) { // for working on sphere
			saveImage(img, filenameSphere);
		} else if(object.equals("Polygon")) { // for working on polygon
			saveImage(img, filenamePolygon);
		}
		
	}

	// Ray Construction method
	static void ray_construction(int i, int j) {
		/* Step 1:
		 * Map (j,i) in the screen coordinates to (xc, yc) in the camera coordinates.
		 */
		float xc = (Xmax - Xmin) * j / (COLS - 1) + Xmin;
		float yc = (Ymax - Ymin) * i / (ROWS - 1) + Ymin;
		
		/* Step 2:
		 * Transform the origin (0.0, 0.0, 0.0) of the camera coordinate to P0 in the world coordinate using the transformation
		 * matrix Mcw. Note that the transformation result should be VRP.
		 */
		Vector Origin = new Vector(0f,0f,0f);
		P0 = Origin.transformC2W();

		/* Step 3:
		 * Transform the point (xc, yc, f) on the image plane in the camera coordinates to P1 in the world coordinates using the
		 * transformation matrix Mcw.
		 */
		Vector P1 = new Vector(xc, yc, Vector.focal).transformC2W();
		
		/* Step 4:
		 * Calculate dirV and normalize it to unit length.
		 */
		Vector dir = new Vector().createVector(P0, P1);
		dirV.copy(dir);
		//normalize dirV
		dirV.normalize();
		
		return;
	}
	
	static Boolean ray_tracing(Vector P0_, Vector dirV_) {
		/* If the ray intersects with any object, this function call will return a true value and the data of the intersection are
		 * returned in P, N and Kd.
		 */
		
		found = ray_object_intersection(P0_, dirV_);
		return found;
	}
	
	static Boolean ray_object_intersection(Vector P0_, Vector dirV_) {
		float t1 = 0.0f, t2 = 0.0f;
		
		if(object.equals("Sphere")) {
			//check intersection with sphere
			t1 = ray_sphere_intersection(P0_, dirV_, S);
		} else if(object.equals("Polygon")) {
			//check intersection with polygon
			t2 = ray_polygon_intersection(P0_, dirV_, PL);
		}
		
		if(t1 == 0.0f && t2 == 0.0f) {
			//no intersection with both the objects
			return false;
		} else if(t2 == 0.0f) {
			//intersection with only the sphere
			
			P = P0.addVector(dirV.multiply(t1));
			N = generateSphereSurfaceNormal(P);
			Kd = S.Kd;
			return true;
			
		} else if(t1 == 0) {
			//intersection with only polygon
			
			P = P0.addVector(dirV.multiply(t2));
			N.copy(polyNormal);;
			Kd = PL.Kd;
			return true;
		} else if (t1 < t2) {
			//intersection with both, but sphere is closer then polygon
			P = P0.addVector(dirV.multiply(t1));
			
			N = generateSphereSurfaceNormal(P);
			Kd = S.Kd;
			return true;
		} else {
			//intersection with both, but polygon is closer then sphere
			P = P0.addVector(dirV.multiply(t2));
			
			N.copy(polyNormal);;
			Kd = PL.Kd;
			return true;
		}
	}
	
	private static float ray_polygon_intersection(Vector p0_, Vector dirV_, Polygon pL2) {
		//generatePlaneEquation();
		generatePolygonNormal(pL2);
		float D = calculateD(pL2.v[0]);
		float t = computeIntersection(dirV_, p0_, D);
		if(flag) {
			
			Vector point = new Vector();
			point = p0_.addVector(dirV_.multiply(t));
			if(checkIsInside(pL2, point)) {
				return t;
			} else {
				return 0.0f;
			}
		} else {
			
			return 0.0f;
		}
	}

	private static Boolean checkIsInside(Polygon pL2_, Vector point) {
		int index = 0;
		if(Math.abs(polyNormal.getX()) >= Math.abs(polyNormal.getY()) && Math.abs(polyNormal.getX()) >= Math.abs(polyNormal.getZ())) {
			index = 0;
		} else if(Math.abs(polyNormal.getY()) >= Math.abs(polyNormal.getX()) && Math.abs(polyNormal.getY()) >= Math.abs(polyNormal.getZ())) {
			index = 1;
		} else {
			index = 2;
		}
		
		Vector[] v2D = new Vector[4];
		for(int i = 0; i < 4; i++) {
			v2D[i] = pL2_.v[i].convert2D(index);
		}
		Vector newPoint = new Vector();
		newPoint = point.convert2D(index);
		polyTranslate.MatrixT(newPoint.getX(), newPoint.getY(), newPoint.getZ());
		
		for(int i = 0; i < 4; i++) {
			v2D[i] = v2D[i].subtractVector(newPoint);//polyTranslate.multiplyVector(v2D[i]);
		}
		int intersection = 0;
		//calculate number of intersections with edges
		switch(index) {
		//Working on the Y-Z plane
		case 0:
			for(int i = 1; i < 4; i++) {
				if(v2D[i].getZ() > 0 && v2D[i-1].getZ() > 0) {
					//both are greater then zero i.e. on the positive side of Z-axis
					if(v2D[i].getY() < 0 && v2D[i-1].getY() < 0) {
						//both are on the positive side of Y-axis. thus the edge does not cross Z-axis
					} else if(v2D[i].getY() < 0 && v2D[i-1].getY() < 0) {
						//both are on the negative side of Y-axis. thus the edge does not cross Z-axis
					} else {
						//one is positive and one is negative. thus, the edge will cross Z-axis
						intersection++;
					}
				} else if(v2D[i].getZ() < 0 && v2D[i-1].getZ() < 0) {
					//both are on negative Z-axis. this is not our concern.
				} else {
					// one is on positive side, other is on negative side of Z-axis.
					intersection++;
				}
			}
			if(v2D[3].getZ() > 0 && v2D[0].getZ() > 0) {
				//both are greater then zero i.e. on the positive side of Z-axis
				if(v2D[3].getY() < 0 && v2D[0].getY() < 0) {
					//both are on the positive side of Y-axis. thus the edge does not cross Z-axis
				} else if(v2D[3].getY() < 0 && v2D[0].getY() < 0) {
					//both are on the negative side of Y-axis. thus the edge does not cross Z-axis
				} else {
					//one is positive and one is negative. thus, the edge will cross Z-axis
					intersection++;
				}
			} else if(v2D[3].getZ() < 0 && v2D[0].getZ() < 0) {
				//both are on negative Z-axis. this is not our concern.
			} else {
				// one is on positive side, other is on negative side of Z-axis.
				intersection++;
			}
			break;
		//Working on the X-Z plane
		case 1:
			for(int i = 1; i < 4; i++) {
				if(v2D[i].getX() > 0 && v2D[i-1].getX() > 0) {
					//both are greater then zero i.e. on the positive side of X-axis
					if(v2D[i].getZ() < 0 && v2D[i-1].getZ() < 0) {
						//both are on the positive side of Y-axis. thus the edge does not cross X-axis
					} else if(v2D[i].getZ() < 0 && v2D[i-1].getZ() < 0) {
						//both are on the negative side of Y-axis. thus the edge does not cross X-axis
					} else {
						//one is positive and one is negative. thus, the edge will cross X-axis
						intersection++;
					}
				}
				else if(v2D[i].getX() < 0 && v2D[i-1].getX() < 0) {
					//both are on negative X-axis. this is not our concern.
				} else {
					// one is on positive side, other is on negative side of X-axis.
					intersection++;
				}
			}
			if(v2D[3].getX() > 0 && v2D[0].getX() > 0) {
				//both are greater then zero i.e. on the positive side of X-axis
				if(v2D[3].getZ() < 0 && v2D[0].getZ() < 0) {
					//both are on the positive side of Y-axis. thus the edge does not cross X-axis
				} else if(v2D[3].getZ() < 0 && v2D[0].getZ() < 0) {
					//both are on the negative side of Y-axis. thus the edge does not cross X-axis
				} else {
					//one is positive and one is negative. thus, the edge will cross X-axis
					intersection++;
				}
			} else if(v2D[3].getX() < 0 && v2D[0].getX() < 0) {
				//both are on negative X-axis. this is not our concern.
			} else {
				// one is on positive side, other is on negative side of X-axis.
				intersection++;
			}
			break;
		//Working on the X-Y plane
		case 2:
			for(int i = 1; i < 4; i++) {
				if(v2D[i].getY() > 0 && v2D[i-1].getY() > 0) {
					//both are greater then zero i.e. on the positive side of Y-axis
					if(v2D[i].getX() < 0 && v2D[i-1].getX() < 0) {
						//both are on the positive side of Y-axis. thus the edge does not cross Y-axis
					} else if(v2D[i].getX() < 0 && v2D[i-1].getX() < 0) {
						//both are on the negative side of Y-axis. thus the edge does not cross Y-axis
					} else {
						//one is positive and one is negative. thus, the edge will cross Y-axis
						intersection++;
					}
				} else if(v2D[i].getY() < 0 && v2D[i-1].getY() < 0) {
					//both are on negative Y-axis. this is not our concern.
				} else {
					// one is on positive side, other is on negative side of Y-axis.
					intersection++;
				}
			}
			if(v2D[3].getY() > 0 && v2D[0].getY() > 0) {
				//both are greater then zero i.e. on the positive side of Y-axis
				if(v2D[3].getX() < 0 && v2D[0].getX() < 0) {
					//both are on the positive side of X-axis. thus the edge does not cross Y-axis
				} else if(v2D[3].getX() < 0 && v2D[0].getX() < 0) {
					//both are on the negative side of X-axis. thus the edge does not cross Y-axis
				} else {
					//one is positive and one is negative. thus, the edge will cross Y-axis
					intersection++;
				}
			} else if(v2D[3].getY() < 0 && v2D[0].getY() < 0) {
				//both are on negative Y-axis. this is not our concern.
			} else {
				// one is on positive side, other is on negative side of Y-axis.
				intersection++;
			}
			break;
		}
		
		if((intersection % 2) == 0) {
			return false;
		} else {
			return true;
		}
	}

	private static float computeIntersection(Vector dir_, Vector p0_, float D) {
		if(polyNormal.dotProduct(dir_) == 0.0f) {
			flag = false;
			return 0;
		} else {
			flag = true;
			float t = -1 * ((polyNormal.dotProduct(p0_) + D) / polyNormal.dotProduct(dir_));
			
			return t;
		}
	}

	private static float calculateD(Vector v0) {
		float D = -1 * polyNormal.dotProduct(v0);// -(polyNormal.getX()*v0.getX() + polyNormal.getY()*v0.getY() + polyNormal.getZ()*v0.getZ());
		
		return D;
	}

	private static void generatePolygonNormal(Polygon pL2) {
		Vector v1 = new Vector().createVector(pL2.v[1], pL2.v[0]);
		Vector v2 = new Vector().createVector(pL2.v[2], pL2.v[0]);
		
		Vector normal = new Vector();
		normal = v1.crossProduct(v2);
		
		normal.normalize();
		polyNormal.copy(normal);
		
	}

	private static Vector generateSphereSurfaceNormal(Vector p2) {
		Vector N = new Vector().createVector(p2, new Vector(S.x, S.y, S.z));
		N.normalize();
		return N;
	}

	private static float ray_sphere_intersection(Vector p0_, Vector dirV_, Sphere s2) {
		Vector L = new Vector().createVector(p0_, new Vector(s2.x, s2.y, s2.z));
		
		float Tca = L.dotProduct(dirV_);
		
		if(Tca < 0) {
			return 0;
		}
		float d = (float) Math.sqrt(L.dotProduct(L) - Math.pow(Tca, 2));
		
		if(d > s2.radius) {
			return 0;
		} else if(d == s2.radius) {
			return Tca;
		} else {
			
			float Thc = ((s2.radius * s2.radius) - (d * d));
			
			float t0 = Tca - Thc;
			float t1 = Tca + Thc;
			if(t0 > t1) {
				float temp = t0;
				t0 = t1;
				t1 = temp;
			}
			if(t0 < 0.0f) {
				t0 = t1;
				if(t0 < 0.0f) {
					return 0;
				}
			}
			return t0;
		}
	}

	static float shading(Vector P_, Vector N_, float Kd_) {
		/* Obtain the light-ray hitting the point.
		 * LRP is the light source point
		 */
		Vector L = new Vector().createVector(P_, LRP); //LRP, P_);
		
		//normalize L to unit length
		L.normalize();
		
		//Calculate the shading value
		float C = (float) (Ip * Kd * (N.dotProduct(L)));	// here (N * L) is a dot product
		
		return C;
	}
	
	// Adjusting Color value for the image pixel
	private static Color brightness(Color c, float scale) {
		
	    int r = Math.min(255, Math.max(0, (int) (c.getRed() +(0.5 * scale))));
	    int g = Math.min(255, Math.max(0, (int) (c.getGreen() + (0.5 * scale))));
	    int b = Math.min(255, Math.max(0, (int) (c.getBlue() + (0.5 * scale))));
	    return new Color(r,g,b);
	}
	
	/* Methods for image creating */
	private static BufferedImage map( int sizeX, int sizeY ){
        final BufferedImage res = new BufferedImage( sizeX, sizeY, BufferedImage.TYPE_INT_RGB );
        for (int x = 0; x < sizeX; x++){
            for (int y = 0; y < sizeY; y++){
                res.setRGB(x, y, Color.BLACK.getRGB() );
            }
        }
        return res;
    }

    private static void saveImage( final BufferedImage bi, final String path ){
        try {
            RenderedImage rendImage = bi;
            //ImageIO.write(rendImage, "bmp", new File(path));
            //ImageIO.write(rendImage, "PNG", new File(path));
            ImageIO.write(rendImage, "jpeg", new File(path));
        } catch ( IOException e) {
            e.printStackTrace();
        }
    }

}
