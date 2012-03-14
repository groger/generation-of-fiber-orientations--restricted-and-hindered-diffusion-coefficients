#ifndef CYLINDER_H
#define CYLINDER_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <Matrix.h>
#include <Vector.h>

class Cylinder{
	public:
		explicit Cylinder(const double height, const double radius) ;
		virtual vtkSmartPointer<vtkPolyData> GetVTKPolyData(const geometry::Vector &axis = geometry::Vector(0, 0, 1),
								    const geometry::Matrix &m = geometry::identity()) const ;
			//This method returns a cylinder whose axis is aligned with the input vector
		virtual ~Cylinder(){}
		
	protected:
		double height, radius ;
} ;

#endif
