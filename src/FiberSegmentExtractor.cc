/*
 *	By Yinpeng Li, mousquetaires@unc.edu
 */

#include <FiberSegmentExtractor.h>
#include <FiberODF_Common.h>
#include <vtkPolyDataReader.h>
#include <vtkPolyData.h>

FiberSegmentExtractor::FiberSegmentExtractor(const char* vtkFiberFile){
	vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New() ;
	reader->SetFileName(vtkFiberFile) ;
	reader->Update() ;
	
	vtkSmartPointer<vtkPolyData> fiberData = reader->GetOutput() ;
	points = fiberData->GetPoints() ;
	fibers = fiberData->GetLines() ;
	
	numSegmentsInEachFiber.resize(GetNumberOfFibers()) ;
	
	//Get the number of segments in each fiber
	RestartTraversal() ;
	vtkIdType pointNum ;
	vtkIdType *points ;
	for (int i = 0; i < GetNumberOfFibers(); i++){
		fibers->GetNextCell(pointNum, points) ;
		numSegmentsInEachFiber[i] = (static_cast<int>(pointNum - 1) >= 0) ? static_cast<int>(pointNum - 1) : 0 ;
	}
	
	//Restart the traversal
	RestartTraversal() ;
}

int FiberSegmentExtractor::GetNumberOfFibers() const{
	return static_cast<int>(fibers->GetNumberOfCells()) ;
}

int FiberSegmentExtractor::GetNumberOfSegments() const{
	int accumu = 0 ;
	for (int i = 0; i < GetNumberOfFibers(); i++)
		accumu += numSegmentsInEachFiber[i] ;
	return accumu ;
}

FiberSegmentExtractor::operator bool() const{
	if (traversalFiberIndex >= GetNumberOfFibers()){
		return false ;
	}
	else{
		bool emptyFromHereOn = true ;
		for (int i = traversalFiberIndex; i < GetNumberOfFibers(); i++)
			if (numSegmentsInEachFiber[i] > 0){
				emptyFromHereOn = false ;
				break ;
			}
		return !emptyFromHereOn ;
	}
}

void FiberSegmentExtractor::RestartTraversal(){
	fibers->InitTraversal() ;
	traversalFiberIndex = 0 ;
	traversalSegmentIndex = 0 ;
	traversalPointIds = NULL ;
}

FiberSegment FiberSegmentExtractor::Next(){
	while (traversalSegmentIndex == 0){		//When the traversal segment index is 0, read in the new fiber
							//At this point, the traversal fiber index is guaranteed to be pointing to a new fiber
		if (traversalFiberIndex >= GetNumberOfFibers())
			throw EndOfTraversal("FiberSegmentExtractor::Next") ;
		
		vtkIdType temp ;
		fibers->GetNextCell(temp, traversalPointIds) ;
		
		if (numSegmentsInEachFiber[traversalFiberIndex] == 0)		//Skip all the empty fibers (fibers without segments)
			traversalFiberIndex++ ;
		else
			break ;
	}
	
	double ptarray[3] ;
	points->GetPoint(traversalPointIds[traversalSegmentIndex], ptarray) ;
	geometry::Point p1(ptarray[0], ptarray[1], ptarray[2]) ;
	points->GetPoint(traversalPointIds[traversalSegmentIndex + 1], ptarray) ;
	geometry::Point p2(ptarray[0], ptarray[1], ptarray[2]) ;
	traversalSegmentIndex++ ;
	
	if (traversalSegmentIndex >= numSegmentsInEachFiber[traversalFiberIndex]){
		traversalFiberIndex++ ;		//As long as the fiber is not fully processed, the traversal fiber index stays on the fiber
						//As soon as the fiber is finished, the traversal fiber index moves on to the next fiber
		traversalSegmentIndex = 0 ;
	}
	
	return FiberSegment(p1, p2) ;
}

FiberSegment FiberSegmentExtractor::Next(bool &fiberEnd){
	while (traversalSegmentIndex == 0){		//When the traversal segment index is 0, read in the new fiber
							//At this point, the traversal fiber index is guaranteed to be pointing to a new fiber
		if (traversalFiberIndex >= GetNumberOfFibers())
			throw EndOfTraversal("FiberSegmentExtractor::Next") ;
		
		vtkIdType temp ;
		fibers->GetNextCell(temp, traversalPointIds) ;
		
		if (numSegmentsInEachFiber[traversalFiberIndex] == 0)		//Skip all the empty fibers (fibers without segments)
			traversalFiberIndex++ ;
		else
			break ;
	}
	
	double ptarray[3] ;
	points->GetPoint(traversalPointIds[traversalSegmentIndex], ptarray) ;
	geometry::Point p1(ptarray[0], ptarray[1], ptarray[2]) ;
	points->GetPoint(traversalPointIds[traversalSegmentIndex + 1], ptarray) ;
	geometry::Point p2(ptarray[0], ptarray[1], ptarray[2]) ;
	traversalSegmentIndex++ ;
	
	if (traversalSegmentIndex >= numSegmentsInEachFiber[traversalFiberIndex]){
        //End of current fiber
        fiberEnd = true;
		traversalFiberIndex++ ;		//As long as the fiber is not fully processed, the traversal fiber index stays on the fiber
						//As soon as the fiber is finished, the traversal fiber index moves on to the next fiber
		traversalSegmentIndex = 0 ;
	}else{
        fiberEnd = false;
    }
	
	return FiberSegment(p1, p2) ;
}


std::ostream &operator<<(std::ostream &os, const FiberSegmentExtractor &fse){
	os << std::fixed ;
	os << "Fiber Segment Extractor:" << std::endl ;
	os << "Number of fibers:" << fse.GetNumberOfFibers() << std::endl ;
	os << "Number of segments:" << fse.GetNumberOfSegments() << std::endl ;
	os << "Number of segments in each fiber:" << std::endl ;
	for (int i = 0; i < fse.GetNumberOfFibers(); i++){
		os << std::setw(7) << i << ":" << std::setw(10) << fse.numSegmentsInEachFiber[i] ;
		if (i % 3 == 2)
			os << std::endl ;
		else
			os << "   " ;
	}
	if (fse.GetNumberOfFibers() % 3 != 0){
		os << std::endl ;
	}
	return os ;
}
