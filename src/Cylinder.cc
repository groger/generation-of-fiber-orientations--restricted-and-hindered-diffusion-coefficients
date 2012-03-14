/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <Cylinder.h>
#include <FiberODF_Common.h>
#include <Geometry_Common.h>
#include <Point.h>
#include <cmath>
#include <vtkPolygon.h>
#include <vtkCellArray.h>

Cylinder::Cylinder(const double height, const double radius){
	if (radius < 0)
		throw InvalidDataRange("Cylinder::Cylinder : radius") ;
	if (height <= 0)
		throw InvalidDataRange("Cylinder::Cylinder : height") ;
	
	this->height = height ;
	if (radius == 0)
		this->radius = geometry::EPSILON ;
	else
		this->radius = radius ;
}

vtkSmartPointer<vtkPolyData> Cylinder::GetVTKPolyData(const geometry::Vector &axis, const geometry::Matrix &m) const{

	//Compute the transform matrix
	geometry::Matrix transform ;
	
	geometry::Vector ax = axis ;
	ax.normalize() ;
	if (geometry::crossProduct(ax, geometry::Vector(0, 0, 1)).magnitude() < geometry::EPSILON){
		//The input cylinder axis is aligned with z-axis
		transform = m ;
	}
	else{
		//Compute the axis vector around which to rotate (0, 0, 1) to the cylinder's axis
		geometry::Vector rotateAxis = geometry::crossProduct(geometry::Vector(0, 0, 1), ax) ;
		rotateAxis.normalize() ;
		double rotateAngle = std::acos(geometry::Vector(0, 0, 1) * ax) / geometry::PI * 180 ;	//In degrees
		
		//Compute the matrix to rotate the rotateAxis to xz plane about z-axis
		const double rotAxisProjOnXY = std::sqrt(rotateAxis[0] * rotateAxis[0] + rotateAxis[1] * rotateAxis[1]) ;
		geometry::Matrix rotXZ(
			rotateAxis[0] / rotAxisProjOnXY,	rotateAxis[1] / rotAxisProjOnXY,	0,	0,
			-rotateAxis[1] / rotAxisProjOnXY,	rotateAxis[0] / rotAxisProjOnXY,	0,	0,
			0,					0,					1,	0,
			0,					0,					0,	1
		) ;
		
		//Compute the matrix to rotate the rotateAxis already in xz plane to align with z-axis around y-axis
		//Note that the rotateAxis is already normalized to unit vector
		geometry::Matrix rotZAxis(
			rotateAxis[2],		0,		-rotAxisProjOnXY,		0,
			0,			1,		0,				0,
			rotAxisProjOnXY,	0,		rotateAxis[2],			0,
			0,			0,		0,				1
		) ;
		
		transform = m * geometry::inverse(rotXZ) * geometry::inverse(rotZAxis) * geometry::rotateAroundZ(rotateAngle) * rotZAxis * rotXZ ;
	}
  
//DEBUG
	//std::cout << geometry::rotateAroundZ(rotateAngle) << std::endl ;
	//std::cout << transform << std::endl ;
//
  
	//The cylinder is approximated by being divided into rectangles
  
	const int NUM_POINTS_ON_CIRCLE = 60 ;
  
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New() ;
	points->SetDataTypeToDouble() ;
	
	for (int i = 0; i < NUM_POINTS_ON_CIRCLE; i++){
		double angle = i * geometry::PI * 2 / NUM_POINTS_ON_CIRCLE ;
		geometry::Point p(std::cos(angle) * radius, std::sin(angle) * radius, - height / 2) ;
		p = transform * p ;
		points->InsertNextPoint(p[0], p[1], p[2]) ;
	}

	for (int i = 0; i < NUM_POINTS_ON_CIRCLE; i++){
		double angle = i * geometry::PI * 2 / NUM_POINTS_ON_CIRCLE ;
		geometry::Point p(std::cos(angle) * radius, std::sin(angle) * radius, height / 2) ;
		p = transform * p ;
		points->InsertNextPoint(p[0], p[1], p[2]) ;
	}

	vtkSmartPointer<vtkCellArray> polygons = vtkSmartPointer<vtkCellArray>::New() ;

	for (int i = 0; i < NUM_POINTS_ON_CIRCLE - 1; i++){
		vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New() ;
		
		polygon->GetPointIds()->InsertNextId(i) ;
		polygon->GetPointIds()->InsertNextId(i + 1) ;
		polygon->GetPointIds()->InsertNextId(NUM_POINTS_ON_CIRCLE + i + 1) ;
		polygon->GetPointIds()->InsertNextId(NUM_POINTS_ON_CIRCLE + i) ;
		
		polygons->InsertNextCell(polygon) ;
	}
	{
		//The last rectangles
		vtkSmartPointer<vtkPolygon> polygon = vtkSmartPointer<vtkPolygon>::New() ;
		
		polygon->GetPointIds()->InsertNextId(NUM_POINTS_ON_CIRCLE - 1) ;
		polygon->GetPointIds()->InsertNextId(0) ;
		polygon->GetPointIds()->InsertNextId(NUM_POINTS_ON_CIRCLE) ;
		polygon->GetPointIds()->InsertNextId(NUM_POINTS_ON_CIRCLE * 2 - 1) ;
		
		polygons->InsertNextCell(polygon) ;
	}

	vtkSmartPointer<vtkPolyData> data = vtkSmartPointer<vtkPolyData>::New() ;
	data->SetPoints(points) ;
	data->SetPolys(polygons) ;

	return data ;
}
