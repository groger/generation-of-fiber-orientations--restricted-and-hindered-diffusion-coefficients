#ifndef SphereIkosohedron_h
#define SphereIkosohedron_h

/*
 * Notice that this is a non-recursive icosahedron subdivision.
 * Adapted from the NeuroLib project @ UNC
 */

#include <iostream>
#include <itkDataObject.h>
#include <itkObjectFactory.h>
#include <vector>
#include <math.h>
#include <sstream>
#include <vtkPolyData.h>
#include <vtkPolyDataWriter.h>
#include <vtkSmartPointer.h>
#include <vtkPolygon.h>
#include <vtkCellArray.h>

typedef std::vector<double>	VectorType ;
typedef std::vector<int>	IndexList ;

namespace itk {

  /*
   * 
   * The Icosahedron created is centered at (0,0,0), and has radius 1
   * 
   */
  
template < typename T >
class ITK_EXPORT SphereIkosahedron : public DataObject
{
	public:

		typedef SphereIkosahedron 		Self;
		typedef DataObject 			Superclass;
		typedef SmartPointer<Self> 		Pointer;
		typedef SmartPointer<const Self> 	ConstPointer;
		
		itkNewMacro(Self);
		itkTypeMacro(SphereIkosahedron, DataObject);
		itkSetMacro(SubdivisionLevel, int);
		itkGetMacro(SubdivisionLevel, int);
		
		void 					Initialize();
		void 					WriteIcosahedronToVTKFile(const char *fname = "Icosahedron.vtk");
		vtkSmartPointer<vtkPolyData>		CreateVTKPolyData() ;
		int 					GetNumberOfVertices() ;
		VectorType 				GetCoordinateTableatIndex(int) ;	//Vertex cartesian coordinates
		VectorType 				GetPhiThetaTableatIndex(int) ;
        std::vector< VectorType > GetCoordinateTable(); //Vertex cartesian coordinates
        std::vector< VectorType > GetPhiThetaTable(); 
		int 					PhiThetaToIndex(double, double) ;
		IndexList				GetTriangleVertices(int) ;
		IndexList				GetSurroundingTriangles(int) ;		//To get the surrounding triangles of a vertex
		const std::vector<VectorType> &		GetTriangle(int) ;
		
	protected:
		void 					ComputeSubdivision();
		SphereIkosahedron(){}
		~SphereIkosahedron(){}
		
		std::vector< VectorType > 		m_CoordinateTable;
		std::vector< VectorType > 		m_PhiThetaTable;			//Phi, first coordinate, theta second one
		std::vector< IndexList >		m_OrdinatedTriangles;
		std::vector<std::vector< VectorType > >	m_all_triangs;
		int 					m_SubdivisionLevel;
		int 					m_NumberOfVertices;
		int 					m_NumberOfTriangles;
		T *					m_FeatureTable;
		std::vector<IndexList>			m_VertexTriangleMap ;
		
	private:
		SphereIkosahedron(const Self&); //purposely not implemented
		void 	operator=(const Self&); //purposely not implemented
};

#define X .525731112119133606
#define Z .850650808352039932

#ifndef M_PI
#define M_PI 3.141592653589793
#endif

//Vertices, triangles, edges of a single icosahedron
static double vert[12][3] = {
	{-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z}, {X, 0.0, -Z},
	{0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
	{Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}
};
static int edge[30][2] = {
	{0,1}, {0,4}, {0,6}, {0,9}, {0,11}, {1,4}, {1,6}, {1,8}, {1,10}, {2,3},
	{2,5}, {2,7}, {2,9}, {2,11}, {3,5}, {3,7}, {3,8}, {3,10}, {4,5}, {4,8},
	{4,9}, {5,8}, {5,9}, {6,7}, {6,10}, {6,11}, {7,10}, {7,11}, {8,10}, {9,11}
};
static int triang[20][3] = {
	{0,4,1}, {0,9,4}, {9,5,4}, {4,5,8}, {4,8,1},
	{8,10,1}, {8,3,10}, {5,3,8}, {5,2,3}, {2,7,3},
	{7,10,3}, {7,6,10}, {7,11,6}, {11,0,6}, {0,1,6},
	{6,1,10}, {9,0,11}, {9,11,2}, {9,2,5}, {7,2,11}
};

//************************************************************************************************

template < typename T >
int SphereIkosahedron< T >::GetNumberOfVertices()
{
	return m_NumberOfVertices;
}

//************************************************************************************************

template < typename T >
void SphereIkosahedron< T >::ComputeSubdivision()
{
	int i;
	
	m_NumberOfVertices = 12 + 30*m_SubdivisionLevel + 20*m_SubdivisionLevel*(m_SubdivisionLevel-1)/2;//12 is the number of vertices of the simple iko, 30: number of edge, 20: number of triangles
	m_NumberOfTriangles = 20*(m_SubdivisionLevel+1)*(m_SubdivisionLevel+1);

	int subdiv = m_SubdivisionLevel + 1;
	int m,n;
	double x1, y1, z1, x2, y2, z2, x3, y3, z3, dx12, dy12, dz12, dx23, dy23, dz23, length;
	for(i=0; i<30; i++) 
	{
		x1 = vert[edge[i][0] ][0];
		y1 = vert[edge[i][0] ][1];
		z1 = vert[edge[i][0] ][2];
		x2 = vert[edge[i][1] ][0];
		y2 = vert[edge[i][1] ][1];
		z2 = vert[edge[i][1] ][2];
		dx12 = (x2 - x1)/subdiv;
		dy12 = (y2 - y1)/subdiv;
		dz12 = (z2 - z1)/subdiv;
		for(n=1; n<subdiv; n++)		//Note that nothing is done when subdiv == 1
		{
			VectorType Point, PhiTheta;
			Point.push_back(x1 + n*dx12);
			Point.push_back(y1 + n*dy12);
			Point.push_back(z1 + n*dz12);
			
			length = sqrt((double) Point[0]*Point[0]+
					Point[1]*Point[1]+
					Point[2]*Point[2]);
			//Normalisation (to make sure each new point is at a distance of 1 from the center so that we approximate a sphere)
			Point[0] /= length;
			Point[1] /= length;
			Point[2] /= length;
			m_CoordinateTable.push_back(Point);
			//From cartesian to polar coordinates
			double temp_phi = atan2(Point[1] , Point[0]);
			while(temp_phi >= M_PI)
			{
				temp_phi -= 2*M_PI;
			}
			while(temp_phi < -M_PI)
			{
				temp_phi += 2*M_PI;
			}
			PhiTheta.push_back(temp_phi);//phi
			PhiTheta.push_back(acos(Point[2]));//theta
		
			m_PhiThetaTable.push_back(PhiTheta);
		}
	}

	if(subdiv > 2) 
	{
		//If the subdivision level is > 2, we have to create points inside each basic triangle
		for(i=0; i<20; i++) 
		{
			//For each triangle, we subdivide each edge.
			x1 = vert[triang[i][0] ][0];
			y1 = vert[triang[i][0] ][1];
			z1 = vert[triang[i][0] ][2];
			x2 = vert[triang[i][1] ][0];
			y2 = vert[triang[i][1] ][1];
			z2 = vert[triang[i][1] ][2];
			x3 = vert[triang[i][2] ][0];
			y3 = vert[triang[i][2] ][1];
			z3 = vert[triang[i][2] ][2];
			dx12 = (x2 - x1)/subdiv;
			dy12 = (y2 - y1)/subdiv;
			dz12 = (z2 - z1)/subdiv;
			dx23 = (x3 - x2)/subdiv;
			dy23 = (y3 - y2)/subdiv;
			dz23 = (z3 - z2)/subdiv;

			n = 1;
			do 
			{
				for(m=1; m<=n; m++) 
				{
					//Then, we subdivide each distance between 2 opposite points from 2 edges of the same triangle
					VectorType Point, PhiTheta;
					Point.push_back(x1 + (n+1)*dx12 + m*dx23);
					Point.push_back(y1 + (n+1)*dy12 + m*dy23);
					Point.push_back(z1 + (n+1)*dz12 + m*dz23);
					//Normalisation
					length = sqrt((double) Point[0]*Point[0]+
							Point[1]*Point[1]+
							Point[2]*Point[2]);
					Point[0] /= length;
					Point[1] /= length;
					Point[2] /= length;
					m_CoordinateTable.push_back(Point);
					//From cartesian to polar coordinates
					double temp_phi = atan2(Point[1] , Point[0]);
					while(temp_phi >= M_PI)
					{
						temp_phi -= 2*M_PI;
					}
					while(temp_phi < -M_PI)
					{
						temp_phi += 2*M_PI;
					}
					PhiTheta.push_back(temp_phi);//phi
					PhiTheta.push_back(acos(Point[2]));//theta
		
					m_PhiThetaTable.push_back(PhiTheta);
				}
				n++;
			}while( n<=(subdiv-2) );
		}
	}
	//Create the triangle from the new made vertices
	if (subdiv > 1) 
	{
		for(i=0; i<20; i++) 
		{
			x1 = vert[triang[i][0] ][0];
			y1 = vert[triang[i][0] ][1];
			z1 = vert[triang[i][0] ][2];
			x2 = vert[triang[i][1] ][0];
			y2 = vert[triang[i][1] ][1];
			z2 = vert[triang[i][1] ][2];
			x3 = vert[triang[i][2] ][0];
			y3 = vert[triang[i][2] ][1];
			z3 = vert[triang[i][2] ][2];
			dx12 = (x2 - x1)/subdiv;
			dy12 = (y2 - y1)/subdiv;
			dz12 = (z2 - z1)/subdiv;
			dx23 = (x3 - x2)/subdiv;
			dy23 = (y3 - y2)/subdiv;
			dz23 = (z3 - z2)/subdiv;

			n = 1;
			do 
			{
				for(m=1; m<=n; m++) 
				{
					// Draw lower triangle
					VectorType triangle1, triangle2, triangle3;
					std::vector< VectorType > zero_triangs;
					
					triangle1.push_back(x1 + n*dx12 + m*dx23);
					triangle1.push_back(y1 + n*dy12 + m*dy23);
					triangle1.push_back(z1 + n*dz12 + m*dz23);
					length = sqrt( triangle1[0]*triangle1[0] + triangle1[1]*triangle1[1] + triangle1[2]*triangle1[2]);
					triangle1[0] /= length;
					triangle1[1] /= length;
					triangle1[2] /= length;
					zero_triangs.push_back(triangle1);
					
					triangle2.push_back(x1 + (n-1)*dx12 + (m-1)*dx23);
					triangle2.push_back(y1 + (n-1)*dy12 + (m-1)*dy23);
					triangle2.push_back(z1 + (n-1)*dz12 + (m-1)*dz23);
					length = sqrt( triangle2[0]*triangle2[0] + triangle2[1]*triangle2[1] + triangle2[2]*triangle2[2]);
					triangle2[0] /= length;
					triangle2[1] /= length;
					triangle2[2] /= length;
					zero_triangs.push_back(triangle2);
					
					triangle3.push_back(x1 + n*dx12 + (m-1)*dx23);
					triangle3.push_back(y1 + n*dy12 + (m-1)*dy23);
					triangle3.push_back(z1 + n*dz12 + (m-1)*dz23);
					length = sqrt( triangle3[0]*triangle3[0] + triangle3[1]*triangle3[1] + triangle3[2]*triangle3[2]);
					triangle3[0] /= length;
					triangle3[1] /= length;
					triangle3[2] /= length;
					zero_triangs.push_back(triangle3);
					
					m_all_triangs.push_back(zero_triangs);
					
					if ( m != n ) 
					{
						// Draw lower left triangle
						VectorType triangle4, triangle5, triangle6;
						std::vector< VectorType > one_triangs;
						
						triangle4.push_back(x1 + n*dx12 + m*dx23);
						triangle4.push_back(y1 + n*dy12 + m*dy23);
						triangle4.push_back(z1 + n*dz12 + m*dz23);
						length = sqrt( triangle4[0]*triangle4[0] + triangle4[1]*triangle4[1] + triangle4[2]*triangle4[2]);
						triangle4[0] /= length;
						triangle4[1] /= length;
						triangle4[2] /= length;
						one_triangs.push_back(triangle4);
						
						triangle5.push_back(x1 + (n-1)*dx12 + m*dx23);
						triangle5.push_back(y1 + (n-1)*dy12 + m*dy23);
						triangle5.push_back(z1 + (n-1)*dz12 + m*dz23);
						length = sqrt( triangle5[0]*triangle5[0] + triangle5[1]*triangle5[1] + triangle5[2]*triangle5[2]);
						triangle5[0] /= length;
						triangle5[1] /= length;
						triangle5[2] /= length;
						one_triangs.push_back(triangle5);
						
						triangle6.push_back(x1 + (n-1)*dx12 + (m-1)*dx23);
						triangle6.push_back(y1 + (n-1)*dy12 + (m-1)*dy23);
						triangle6.push_back(z1 + (n-1)*dz12 + (m-1)*dz23);
						length = sqrt( triangle6[0]*triangle6[0] + triangle6[1]*triangle6[1] + triangle6[2]*triangle6[2]);
						triangle6[0] /= length;
						triangle6[1] /= length;
						triangle6[2] /= length;
						one_triangs.push_back(triangle6);
						
						m_all_triangs.push_back(one_triangs);
					}
				}
				n++;
			} while( n<=subdiv );
		}
	}
	
	if (m_SubdivisionLevel == 0){	//0-level subdivision
		for (int i = 0; i < m_NumberOfTriangles; i++){
			std::vector<VectorType> triangle(3) ;
			
			for (int j = 0; j < 3; j++){
				triangle[j].resize(3) ;	//The point
				for (int k = 0; k < 3; k++){
					triangle[j][k] = vert[triang[i][j]][k] ;
				}
			}
			
			m_all_triangs.push_back(triangle) ;
		}
	}
	
	double epsilon = 0.00001;
	
	for (i = 0 ; i < m_NumberOfTriangles ; i++)
	{
		std::vector<int> triangs;
		for( int l = 0; l < 3 ; l++ )
		{
			triangs.push_back(-1);
		}
		m_OrdinatedTriangles.push_back(triangs);
	}

	// find indexes
	if (subdiv > 1){
		for(i = 0; i < m_NumberOfVertices; i++) 
		{
			for (int j = 0; j < m_NumberOfTriangles ; j++) 
			{
				for(int k = 0 ; k < 3 ; k++)
				{
					if (m_OrdinatedTriangles[j][k] < 0)
					{
						
						if ( (fabs(m_CoordinateTable[i][0] - m_all_triangs[j][k][0]) < epsilon) && 
											(fabs(m_CoordinateTable[i][1] - m_all_triangs[j][k][1]) < epsilon) && 
											(fabs(m_CoordinateTable[i][2] - m_all_triangs[j][k][2]) < epsilon ) ) 
						{
							m_OrdinatedTriangles[j][k] = i;
						}
					}
				}
			}
		}
	}
	else{	//Zero-level subdivision
		for (int i = 0; i < m_NumberOfTriangles; i++){
			for (int j = 0; j < 3; j++){
				m_OrdinatedTriangles[i][j] = triang[i][j] ;
			}
		}
	}

	// Compute the vertex-triangle map
	m_VertexTriangleMap.resize(m_NumberOfVertices) ;
	for (int i = 0; i < m_NumberOfTriangles; i++){
		for (int j = 0; j < 3; j++){
			m_VertexTriangleMap[m_OrdinatedTriangles[i][j]].push_back(i) ;
		}
	}
}

//************************************************************************************************

template < typename T >
void SphereIkosahedron< T >::Initialize()
{
	int i = 0;
	
	for( i = 0 ; i < 12 ; i++ )
	{
		VectorType PhiTheta, Coordinate;
		
		Coordinate.push_back(vert[i][0]);
		Coordinate.push_back(vert[i][1]);
		Coordinate.push_back(vert[i][2]);
		m_CoordinateTable.push_back(Coordinate);
		
		double temp_phi = atan2(Coordinate[1] , Coordinate[0]);
		while(temp_phi >= M_PI)
		{
			temp_phi -= 2*M_PI;
		}
		while(temp_phi < -M_PI)
		{
			temp_phi += 2*M_PI;
		}
		PhiTheta.push_back(temp_phi);//phi
		PhiTheta.push_back(acos(Coordinate[2]));//theta
		m_PhiThetaTable.push_back(PhiTheta);
	}

	ComputeSubdivision();
	WriteIcosahedronToVTKFile();

}

//************************************************************************************************


template< typename T >
int SphereIkosahedron< T >::PhiThetaToIndex(double Phi, double Theta)
{
	int indice = -1;
	int i,j;
	
	for(i = 0 ; i < m_NumberOfVertices ; i++)
	{
		if( Phi == m_PhiThetaTable[i][0] && Theta == m_PhiThetaTable[i][1])
		{
			indice = i;
		}
	}
	if(indice == -1)
	{
		std::cout<<"Values of Phi and/or Theta not correct, index not found."<<std::endl;
	}
	return indice;
}

//************************************************************************************************

template < typename T >
std::vector<double> SphereIkosahedron< T >::GetPhiThetaTableatIndex(int index)
{
	return m_PhiThetaTable[index];
}

//************************************************************************************************

template < typename T >
std::vector< VectorType > SphereIkosahedron< T >::GetPhiThetaTable() 
{
	return m_PhiThetaTable;
}


//************************************************************************************************

template<typename T>
vtkSmartPointer<vtkPolyData> SphereIkosahedron<T>::CreateVTKPolyData(){
  
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New() ;
	points->SetDataTypeToDouble() ;
	
	for (int i = 0; i < m_NumberOfVertices; i++){
		points->InsertNextPoint(m_CoordinateTable[i][0], m_CoordinateTable[i][1], m_CoordinateTable[i][2]) ;
	}

	vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New() ;

	for (int i = 0; i < m_NumberOfTriangles; i++){
		vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New() ;
		
		for (int j = 0; j < 3; j++){
			polygon->GetPointIds()->InsertNextId(m_OrdinatedTriangles[i][j]) ;
		}
		
		polygons->InsertNextCell(polygon) ;
	}

	vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New() ;
	data->SetPoints(points) ;
	data->SetPolys(polygons) ;

	return data ;
}

//************************************************************************************************

template < typename T >
void SphereIkosahedron< T >::WriteIcosahedronToVTKFile(const char *fname){

	vtkSmartPointer<vtkPolyData> data = CreateVTKPolyData() ;
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New() ;
	writer->SetInput(data) ;
	writer->SetFileName(fname) ;
	std::stringstream ss ;
	ss << "Icosahedron subdivision level: " << m_SubdivisionLevel ;
	writer->SetHeader(ss.str().c_str()) ;
	writer->Write() ;
	
}

//************************************************************************************************

template < typename T >
std::vector<double> SphereIkosahedron< T >::GetCoordinateTableatIndex(int index)
{
	return m_CoordinateTable[index];
}

//************************************************************************************************

template < typename T >
std::vector< VectorType > SphereIkosahedron< T >::GetCoordinateTable() 
{
	return m_CoordinateTable;
}

//************************************************************************************************

template <typename T>
IndexList SphereIkosahedron<T>::GetTriangleVertices(int triangle){
	return m_OrdinatedTriangles[triangle] ;
}

//************************************************************************************************

template <typename T>
IndexList SphereIkosahedron<T>::GetSurroundingTriangles(int vertex){
	return m_VertexTriangleMap[vertex] ;
}

//************************************************************************************************

template <typename T>
const std::vector<VectorType> &SphereIkosahedron<T>::GetTriangle(int triangle){
	return m_all_triangs[triangle] ;
}

//************************************************************************************************

}//end of namespace itk
#endif
