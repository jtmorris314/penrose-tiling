/* Package penrose implements a basic penrose tiling scheme based
   on subdividing two triangle types ("blue" and "red") that are 1/2
   of the fat and thin rhombi tiles, respectively.  After the triangles
   are subdivided, the similarly colored adjoing tiles are joined to
   create the rhombi.

   This go package is inspired from the python code and explanation
   found here:
      http://preshing.com/20110831/penrose-tiling-explained/

   The seed pattern is a wheel of 10 "red" triangles.

   The program currently iterates 10 times (10 generations of subdividing
   triangles.
*/
package main

import (
	"fmt"
	"github.com/ungerik/go-cairo"
	"math"
	"math/cmplx"
	"sort"
	"flag"
)

// ColorType represents the type of rhombi ("red" is thin and
// "blue" is thick).
type ColorType int

// Symbols for "red" and "blue" rhombi/triangles
const (
	RED ColorType = iota
	BLUE
)

const delta float64 = 0.000001

type Point complex128
type Polygon []Point

// GoldenRatio of yore is represented as complex128 number.
var goldenRatio Point = Point(complex((1+math.Sqrt(5))/2, 0))

// Triangle maintains the triangle type and the position
// of the three vertices: A, B, C.  The "x" coordinate is
// the real part and the "y" coordinate is the imaginary part.
type Triangle struct {
	color ColorType
	pts Polygon
}

// createSeedTriangles builds a wheel of red triangles around
// the origin and returns the wheel as a slice of triangles.
func createSeedTriangles() []Triangle {
	ts := make([]Triangle, 0, 10)
	for i := 0; i < 10; i++ {
		b := Point(cmplx.Rect(1, float64(2*i-1)*math.Pi/float64(10)))
		c := Point(cmplx.Rect(1, float64(2*i+1)*math.Pi/float64(10)))
		if i%2 == 0 {
			b, c = c, b
		}
		t := Triangle{color: RED, pts: Polygon{0i, b, c}}
		ts = append(ts, t)
	}
	return ts
}


// subdivideTriangles breaks each triangle into component triangles.
// The "red" triangles result in one "blue" and one "red" triangle.
// The "blue" triangles result in two "blue" and one "red" triangle.
// The dividing lines are calculated using the golden ratio.
func subdivideTriangles(ts []Triangle) []Triangle {
	result := make([]Triangle, 0, 10)
	for _, t := range ts {
		if t.color == RED {
			P := t.pts[0] + (t.pts[1]-t.pts[0])/goldenRatio
			result = append(result, 
					Triangle{RED, 
						Polygon{t.pts[2], 
							P, 
							t.pts[1]}},
					Triangle{BLUE, 
						Polygon{P, 
							t.pts[2], 
							t.pts[0]}})
		} else {
			Q := t.pts[1] + (t.pts[0]-t.pts[1])/goldenRatio
			R := t.pts[1] + (t.pts[2]-t.pts[1])/goldenRatio
			result = append(result, 
					Triangle{BLUE, 
						Polygon{R, 
							t.pts[2], 
							t.pts[0]}},
					Triangle{BLUE, Polygon{Q, R, t.pts[1]}},
					Triangle{RED, Polygon{R, Q, t.pts[0]}})
		}
	}
	return result
}

// ByBC represents a list of triangles, implements the sort
// interface and sorts by the B vertices first and C vertices
// second.
type ByBC []Triangle

// Len returns length of list.
func (ts ByBC) Len() int { return len(ts) }

// Swap exchanges positions of two triangles.
func (ts ByBC) Swap(i, j int) { ts[i], ts[j] = ts[j], ts[i] }

// Less compares two triangles and returns true if the first
// is less then the second if B1.x < B2.x, B1.y < B2.y,
// C1.x < C2.x and C1.y < C2.y.
//
// The coordinates are float64 and the comparison allows
// for <delta tolerance.
func (ts ByBC) Less(i, j int) bool {
	switch {
	case real(ts[i].pts[1])-real(ts[j].pts[1]) < -delta:
		return true
	case real(ts[i].pts[1])-real(ts[j].pts[1]) > delta:
		return false
	case imag(ts[i].pts[1])-imag(ts[j].pts[1]) < -delta:
		return true
	case imag(ts[i].pts[1])-imag(ts[j].pts[1]) > delta:
		return false
	case real(ts[i].pts[2])-real(ts[j].pts[2]) < -delta:
		return true
	case real(ts[i].pts[2])-real(ts[j].pts[2]) > delta:
		return false
	case imag(ts[i].pts[2])-imag(ts[j].pts[2]) < -delta:
		return true
	}
	return false
}

// Rhombus represents the type of rhombus ("red"/"blue" or
// "thick"/"thin" respectively) and it's four vertices.
type Rhombus struct {
	color ColorType
	pts Polygon
}

// mergeRhombi scans the sorted list of triangles and
// merges triangles with common BC segments.
//
// Unmatched triangles around the periphery are ignored.
//
// The BC matches are within <delta tolerance to account
// for float64 imprecision.
func mergeRhombi(ts []Triangle) []Rhombus {
	rs := make([]Rhombus, 0, 10)
	for i := 0; i < len(ts)-1; i++ {
		if ts[i].color == ts[i+1].color &&
			cmplx.Abs(complex128(ts[i].pts[1]-ts[i+1].pts[1])) < delta &&
			cmplx.Abs(complex128(ts[i].pts[2]-ts[i+1].pts[2])) < delta {
			rs = append(rs, Rhombus{color: ts[i].color,
				pts: Polygon{ts[i].pts[0], ts[i].pts[1],
						ts[i+1].pts[0], ts[i].pts[2]}})
			i = i + 1
		}
	}
	return rs
}


// AB <>  P <> CD and BC <> P <> DA
func polygonContains(bounds Polygon, point Point) bool {
	x := real(point)
	y := imag(point)

	inside := false
        j := len(bounds) - 1
	for i, v := range bounds {
		xi := real(v)
		xj := real(bounds[j])
		yi := imag(v)
		yj := imag(bounds[j])
		
		intersect := (yi > y) != (yj > y) && (x < (xj - xi) * (y - yi) / (yj - yi) + xi)
		if intersect {
			inside = !inside
		}

		j = i
	}
	return inside
/*
	return (real(point) >= real(bounds[0])+delta && real(point) <= real(bounds[2])-delta) && (imag(point) <= imag(bounds[0])-delta && imag(point) >= imag(bounds[2])+delta)
*/
}

// checkInbounds
//
// TBD: need to move call to this function to process rhomboids before drawing.
// TBD: need to calculate "truncated" rhomboids -- might want to key off
// "good" points in/around this function since we've just calculated them.  Could
// return good points and then have another func create new "truncated" rhomboid.
func checkInbounds(bounds Polygon, rhombus Polygon) (bool, bool) {
	contains := true
	overlaps := false
	for _, v := range rhombus {
		if polygonContains(bounds, v) {
			overlaps = true	
		} else {
			contains = false
		}
	}
	return contains, overlaps
}

func lineBoundaryLineIntersect(a1, a2, b1, b2 Point) (bool, Point) {
	ax := real(a1) - real(a2)
	ay := imag(a1) - imag(a2)
	bx := real(b1) - real(b2)
	by := imag(b1) - imag(b2)

	divisor := ax * by - ay * bx

	if (math.Abs(divisor) < delta) {
		return false, 0i
	}
	
	t1 := real(a1) * imag(a2) - imag(a1) * real(a2)
	t2 := real(b1) * imag(b2) - imag(b1) * real(b2)

	x := (t1 * bx - t2 * ax) / divisor
	y := (t1 * by - t2 * ay) / divisor 
	p := Point(complex(x, y))

	invertx := 1.0
	if real(b1) > real(b2) {
		invertx = -1.0
	}
	inverty := 1.0
	if imag(b1) < imag(b2) {
		inverty = -1.0
	}
	a := Point(complex(real(b1) - invertx * delta, imag(b1) + inverty * delta))
	b := Point(complex(real(b2) + invertx * delta, imag(b1) + inverty * delta))
	c := Point(complex(real(b2) + invertx * delta, imag(b2) - inverty * delta))
	d := Point(complex(real(b1) - invertx * delta, imag(b2) - inverty * delta))

	lineSegmentBounds := Polygon{a, b, c, d}
	if polygonContains(lineSegmentBounds, p) {
		return true, p
	}

	return false, 0i

}

func lineBoundsIntersect(bounds Polygon, p1, p2 Point) (bool, Point) {
        j := len(bounds) - 1
	for i, vi := range bounds {
		vj := bounds[j]
		if intersect, p := lineBoundaryLineIntersect(vi, vj, p1, p2); 
				intersect == true {
			return true, p
		}
		j = i;
	}
	return false, 0i
}

// truncateRhombus
//
// Find vertex out-of-bounds.  Find intersection with bounds between it's
// two connected sides.  Add intersection points to make new polygon.
// TBD: think about boundary cases (vertex on bounds, vertex and connecting
// side on bounds ... and corner case -- which we might just ignore)
func truncateRhombus(bounds Polygon, rhombus Polygon) Polygon {
	var polygon Polygon
        j := len(rhombus) - 1
	for i, vi := range rhombus {
		vj := rhombus[j]
		iIn := polygonContains(bounds, vi)
		jIn := polygonContains(bounds, vj)
		if iIn != jIn {
			intersects, pt := lineBoundsIntersect(bounds, vi, vj)
			if intersects {
				polygon = append(polygon, pt)
			} else {
				fmt.Println("ouch", bounds, polygon, pt, vi, vj)
			}
		}
		if iIn {
			polygon = append(polygon, vi)
		} 

		j = i;
	}
	return polygon
}

func clipToBounds(rs []Rhombus, width, height, scale float64) []Rhombus {
	var crs []Rhombus
	divisor := MM_PER_PT*2*scale
	bounds := Polygon{Point(complex(-width/divisor,+height/divisor)),
				Point(complex(+width/divisor,+height/divisor)),
				Point(complex(+width/divisor,-height/divisor)),
				Point(complex(-width/divisor,-height/divisor))}
	for _, r := range rs {
		contains, overlaps := checkInbounds(bounds, r.pts)
		if contains {
			crs = append(crs, r)
		} else if overlaps {
			newPoly := truncateRhombus(bounds, r.pts)
			crs = append(crs, Rhombus{color: r.color, pts: newPoly})
		}
	}

	return crs
}

func insetRhombi(rs []Rhombus, redInsetRatio, blueInsetRatio float64) []Rhombus {
	var irs []Rhombus
	for _, r := range rs {
		insetRatio := Point(complex(redInsetRatio, 0))
		if r.color == BLUE {
			insetRatio = Point(complex(blueInsetRatio, 0))
		}
		newA := r.pts[0] + (r.pts[2] - r.pts[0]) * insetRatio
		newB := r.pts[1] + (r.pts[3] - r.pts[1]) * insetRatio
		newC := r.pts[2] + (r.pts[0] - r.pts[2]) * insetRatio
		newD := r.pts[3] + (r.pts[1] - r.pts[3]) * insetRatio
		newRhombus := Rhombus{color: r.color, pts: Polygon{newA, newB, newC, newD}}
		irs = append(irs, newRhombus)
	}
	return irs
}

// drawRhombi draws the list of rhombi at the indicated scaling factor
//
func drawRhombi(fn string, rs []Rhombus, params Params) {
	s := cairo.NewPSSurface(fn, 
				(params.width+1)/MM_PER_PT, 
				(params.height+1)/MM_PER_PT, 
				cairo.PS_LEVEL_3)

	// Draw bounding box in red
	if params.boundingBox {
		s.SetSourceRGB(1, 0, 0)
		s.SetLineWidth(0.001)		// laser: has to be this value to cut 
		s.Rectangle(1/MM_PER_PT, 
				1/MM_PER_PT, 
				(params.width-1)/MM_PER_PT, 
				(params.height-1)/MM_PER_PT)
		s.Stroke()
	}

	s.Translate((params.width+1)/MM_PER_PT/2, 
			(params.height+1)/MM_PER_PT/2)
	s.Scale(params.scale, params.scale)
	s.SetLineWidth(0.001/params.scale)	// laser: has to be this value to cut 
	s.SetSourceRGB(0, 0, 0)

	for _, r := range rs {
		for i, v := range r.pts {
			if i == 0 {
				s.MoveTo(real(v), imag(v))
			} else {
				s.LineTo(real(v), imag(v))
			}
		}
		s.ClosePath()
		s.Stroke()
	}

	s.Finish()
}

// scaleToRhombusSide returns the scaling factor needed to produce
// rhombi with sides of size in points.
func scaleToRhombusSide(sideInGenericUnits, sideInPoints float64) float64 {
	return sideInPoints / sideInGenericUnits
}

// redBlueInsetRatios calculates the inset ratios from the vertices of the
// "red" and "blue" rhombi to achieve the specified inset from each side.
func redBlueInsetRatios(sideLength, inset float64) (redInsetRatio, blueInsetRatio float64) {
	redOpposite := sideLength * math.Sin(2*math.Pi/10)
	sizeRatio := (redOpposite - 2*inset) / redOpposite
	redInsetRatio = (1 - sizeRatio) / 2

	blueOpposite := sideLength * math.Sin(2*math.Pi/(10/3))
	sizeRatio = (blueOpposite - 2*inset) / blueOpposite
	blueInsetRatio = (1 - sizeRatio) / 2

	return redInsetRatio, blueInsetRatio
}

const MM_PER_PT float64 = 0.352777778                      // size of point in millimeters
const RHOMBUS_SIDE_IN_MM float64 = 3.0                     // length of sides for the rhombi
const SCREEN_WIDTH_IN_MM float64 = 0.75                    // width of screen's "wire"

type Params struct {
	width, height, margin float64
	boundingBox bool
	rhombusSide float64
        rhombusSeparation float64
	scale float64
}

// main creates the seed triangles, subdivides the triangles
// through 10 generations, sorts the triangles so correpsonding
// triangles are adjacent in the list of triangles,
// merges the triangles into rhombi, sets the
// target size of the rhombi and border size and finally
// draws the rhombi to "penrose.ps".
func main() {
        params := Params{}
	flag.Float64Var(&params.rhombusSide, "side", 
			RHOMBUS_SIDE_IN_MM, "mm's of space between rhombii")
	flag.Float64Var(&params.rhombusSeparation, "separation", 
			SCREEN_WIDTH_IN_MM, "mm's of space between rhombii")
	flag.BoolVar(&params.boundingBox, "outline?", false, "true to cut along bounds")
	flag.Float64Var(&params.height, "height", 610.0, "mm's of panel height")
	flag.Float64Var(&params.margin, "margin", 5, "mm's of panel margin")
	flag.Float64Var(&params.width, "width", 381.0, "mm's of panel width")
        fnamePtr := flag.String("file", "penrose.ps", "output file")
	gensPtr := flag.Int("gens", 10, "number of generations")

        flag.Parse()
        fmt.Println("fname", *fnamePtr)
	fmt.Println("gens", *gensPtr)
	fmt.Println("params", params)

	fmt.Printf("Starting penrose tiling.\n")
	ts := createSeedTriangles()
	for i := 0; i < *gensPtr; i++ {
		ts = subdivideTriangles(ts)
	}
	fmt.Println("ts len", len(ts))

	sort.Sort(ByBC(ts))
	rs := mergeRhombi(ts)
	fmt.Println("rs len", len(rs))

	length := math.Sqrt(math.Pow(real(ts[0].pts[0])-real(ts[0].pts[1]), 2) +
		math.Pow(imag(ts[0].pts[0])-imag(ts[0].pts[1]), 2))
	params.scale = scaleToRhombusSide(length, params.rhombusSide/MM_PER_PT)
	redInsetRatio, blueInsetRatio := redBlueInsetRatios(params.rhombusSide, 
								params.rhombusSeparation/2)
	rs = insetRhombi(rs, redInsetRatio, blueInsetRatio)
	fmt.Println("inset rs len", len(rs))

	rs = clipToBounds(rs, params.width - 2 * params.margin, 
				params.height - 2 * params.margin, params.scale)
	fmt.Println("clip rs len", len(rs))

	drawRhombi(*fnamePtr, rs, params)
	fmt.Printf("Ending penrose tiling.\n")
}
