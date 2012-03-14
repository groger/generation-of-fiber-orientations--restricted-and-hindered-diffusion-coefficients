#ifndef FIBERSEGMENTEXTRACTOR_H
#define FIBERSEGMENTEXTRACTOR_H

/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <FiberSegment.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

class FiberSegmentExtractor{
	/*
	 * 
	 * NOTICE that the FiberSegmentExtractor employs RAS coordinate system
	 * 
	 */
public:
	explicit FiberSegmentExtractor(const char *vtkFiberFile) ;
	
	int GetNumberOfFibers() const ;
	int GetNumberOfSegments() const ;
	
	operator bool() const ;		//To test if there are still segments left
	
	void RestartTraversal() ;
	
	FiberSegment Next() ;
	FiberSegment Next(bool &) ;
	
	virtual ~FiberSegmentExtractor(){}
	
protected:
	vtkSmartPointer<vtkPoints> points ;
	vtkSmartPointer<vtkCellArray> fibers ;
	
	std::vector<int> numSegmentsInEachFiber ;
	
	int traversalFiberIndex ;
	int traversalSegmentIndex ;
	vtkIdType *traversalPointIds ;
	
	friend std::ostream &operator<<(std::ostream &os, const FiberSegmentExtractor &fse) ;
} ;

#endif
