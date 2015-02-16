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
)

// GoldenRatio of yore is represented as complex128 number.
var goldenRatio complex128 = complex((1+math.Sqrt(5))/2, 0)

// ColorType represents the type of rhombi ("red" is thin and
// "blue" is thick).
type ColorType int

// Symbols for "red" and "blue" rhombi/triangles
const (
	RED ColorType = iota
	BLUE
)

// Triangle maintains the triangle type and the position
// of the three vertices: A, B, C.  The "x" coordinate is
// the real part and the "y" coordinate is the imaginary part.
type Triangle struct {
	color ColorType
	A     complex128
	B     complex128
	C     complex128
}

// createSeedTriangles builds a wheel of red triangles around
// the origin and returns the wheel as a slice of triangles.
func createSeedTriangles() []Triangle {
	ts := make([]Triangle, 0, 10)
	for i := 0; i < 10; i++ {
		b := cmplx.Rect(1, float64(2*i-1)*math.Pi/float64(10))
		c := cmplx.Rect(1, float64(2*i+1)*math.Pi/float64(10))
		if i%2 == 0 {
			b, c = c, b
		}
		t := Triangle{color: RED, A: 0i, B: b, C: c}
		ts = append(ts, t)
	}
	return ts
}

// drawTriangles outputs a postscript file with the triangles
// filled and and the rhombi outlined.  The destination file
// name is the first parameter and the slice of triangles to
// draw is the second.
func drawTriangles(fn string, ts []Triangle) {
	s := cairo.NewPSSurface(fn, 1000, 1000, cairo.PS_LEVEL_3)
	s.Translate(500, 500)
	s.Scale(500, 500)
	s.SetLineWidth(0.001)
	s.SetLineJoin(cairo.LINE_JOIN_ROUND)

	for _, t := range ts {
		s.MoveTo(real(t.A), imag(t.A))
		s.LineTo(real(t.B), imag(t.B))
		s.LineTo(real(t.C), imag(t.C))
		s.ClosePath()
		if t.color == RED {
			s.SetSourceRGB(1.0, 0.35, 0.35)
		} else {
			s.SetSourceRGB(0.4, 0.4, 1.0)
		}
		s.FillPreserve()
		s.Stroke()

		s.MoveTo(real(t.C), imag(t.C))
		s.LineTo(real(t.A), imag(t.A))
		s.LineTo(real(t.B), imag(t.B))
		s.SetSourceRGB(0.2, 0.2, 0.2)
		s.SetLineWidth(0.001)
		s.Stroke()
	}
	s.Finish()
}

// printTriangle prints the triangle type (red/blue) and coordinates
// of each vertex.
func printTriangle(t Triangle) {
	fmt.Printf("Color=%d, A=%g, B=%g, C=%g\n",
		t.color, t.A, t.B, t.C)
}

// printTriangles dumps the values associated with every
// triangle in the slice of triangles.
func printTriangles(ts []Triangle) {
	for _, t := range ts {
		printTriangle(t)
	}
}

// subdivideTriangles breaks each triangle into component triangles.
// The "red" triangles result in one "blue" and one "red" triangle.
// The "blue" triangles result in two "blue" and one "red" triangle.
// The dividing lines are calculated using the golden ratio.
func subdivideTriangles(ts []Triangle) []Triangle {
	result := make([]Triangle, 0, 10)
	for _, t := range ts {
		if t.color == RED {
			P := t.A + (t.B-t.A)/goldenRatio
			result = append(result, Triangle{RED, t.C, P, t.B},
				Triangle{BLUE, P, t.C, t.A})
		} else {
			Q := t.B + (t.A-t.B)/goldenRatio
			R := t.B + (t.C-t.B)/goldenRatio
			result = append(result, Triangle{BLUE, R, t.C, t.A},
				Triangle{BLUE, Q, R, t.B},
				Triangle{RED, R, Q, t.A})
		}
	}
	return result
}

// ByBC represents a list of triangles, implements the sort
// interface and sorts by the B vertex first and C vertex
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
// for <0.0001 tolerance.
func (ts ByBC) Less(i, j int) bool {
	switch {
	case real(ts[i].B)-real(ts[j].B) < -0.0001:
		return true
	case real(ts[i].B)-real(ts[j].B) > 0.0001:
		return false
	case imag(ts[i].B)-imag(ts[j].B) < -0.0001:
		return true
	case imag(ts[i].B)-imag(ts[j].B) > 0.0001:
		return false
	case real(ts[i].C)-real(ts[j].C) < -0.0001:
		return true
	case real(ts[i].C)-real(ts[j].C) > 0.0001:
		return false
	case imag(ts[i].C)-imag(ts[j].C) < -0.0001:
		return true
	}
	return false
}

// Rhombus represents the type of rhombus ("red"/"blue" or
// "thick"/"thin" respectively) and it's four vertices.
type Rhombus struct {
	color ColorType
	A     complex128
	B     complex128
	C     complex128
	D     complex128
}

// mergeRhombi scans the sorted list of triangles and
// merges triangles with common BC segments.
//
// Unmatched triangles around the periphery are ignored.
//
// The BC matches are within <0.0001 tolerance to account
// for float64 imprecision.
func mergeRhombi(ts []Triangle) []Rhombus {
	rs := make([]Rhombus, 0, 10)
	for i := 0; i < len(ts)-1; i++ {
		if ts[i].color == ts[i+1].color &&
			cmplx.Abs(ts[i].B-ts[i+1].B) < 0.0001 &&
			cmplx.Abs(ts[i].C-ts[i+1].C) < 0.0001 {
			rs = append(rs, Rhombus{color: ts[i].color,
				A: ts[i].A,
				B: ts[i].B,
				C: ts[i+1].A,
				D: ts[i].C})
			i = i + 1
		}
	}
	return rs
}

// drawRhombi draws the rhombi in the list of rhombi.
//
// TBD: handle scaling like drawInsertRhombi.
func drawRhombi(fn string, rs []Rhombus, scaleFactor float64) {
	s := cairo.NewPSSurface(fn, 1000, 1000, cairo.PS_LEVEL_3)
	s.Translate(500, 500)
	s.Scale(scaleFactor, scaleFactor)
	s.SetLineWidth(0.001)
	s.SetLineJoin(cairo.LINE_JOIN_ROUND)

	for _, r := range rs {
		s.MoveTo(real(r.A), imag(r.A))
		s.LineTo(real(r.B), imag(r.B))
		s.LineTo(real(r.C), imag(r.C))
		s.LineTo(real(r.D), imag(r.D))
		s.ClosePath()
		if r.color == RED {
			s.SetSourceRGB(1.0, 0.35, 0.35)
		} else {
			s.SetSourceRGB(0.4, 0.4, 1.0)
		}
		s.FillPreserve()
		s.SetSourceRGB(0.2, 0.2, 0.2)
		s.Stroke()
	}
	s.Finish()
}

// drawInsetRhombi draws the list of rhombi at the indicated scaling factor
// inset according to the ratios for the "red" and "blue" rhombi.
//
// TBD: move sizing of surface to caller.
func drawInsetRhombi(fn string, rs []Rhombus, scaleFactor, redInsetRatio, blueInsetRatio float64) {
	s := cairo.NewPSSurface(fn, 250/MM_PER_PT, 500/MM_PER_PT, cairo.PS_LEVEL_3)
	s.Translate(250/MM_PER_PT/2, 500/MM_PER_PT/2)
	s.Scale(scaleFactor, scaleFactor)
	s.SetLineWidth(0.001)
	s.SetLineJoin(cairo.LINE_JOIN_ROUND)

	for _, r := range rs {
		insetRatio := redInsetRatio
		s.MoveTo(real(r.A), imag(r.A))
		s.LineTo(real(r.B), imag(r.B))
		s.LineTo(real(r.C), imag(r.C))
		s.LineTo(real(r.D), imag(r.D))
		s.ClosePath()
		s.SetSourceRGB(0.35, 0.33, 0.10)
		s.Fill()
		if r.color == BLUE {
			insetRatio = blueInsetRatio
		}
		newA := r.A + complex((real(r.C)-real(r.A)) * insetRatio,
			(imag(r.C)-imag(r.A)) * insetRatio)
		newB := r.B + complex((real(r.D)-real(r.B)) * insetRatio,
			(imag(r.D)-imag(r.B)) * insetRatio)
		newC := r.C + complex((real(r.A)-real(r.C)) * insetRatio,
			(imag(r.A)-imag(r.C)) * insetRatio)
		newD := r.D + complex((real(r.B)-real(r.D)) * insetRatio,
			(imag(r.B)-imag(r.D)) * insetRatio)
		s.MoveTo(real(newA), imag(newA))
		s.LineTo(real(newB), imag(newB))
		s.LineTo(real(newC), imag(newC))
		s.LineTo(real(newD), imag(newD))
		s.ClosePath()
		s.SetSourceRGB(1.0, 1.0, 1.0)
		s.FillPreserve()
		s.SetSourceRGB(0.2, 0.2, 0.2)
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
const RHOMBUS_INSET_IN_MM float64 = SCREEN_WIDTH_IN_MM / 2 // width of wire overlapping the rhombi's side

// main creates the seed triangles, subdivides the triangles
// through 10 generations, sorts the triangles so correpsonding
// triangles are adjacent in the list of triangles,
// merges the triangles into rhombi, sets the
// target size of the rhombi and border size and finally
// draws the rhombi to "penrose.ps".
func main() {
	fmt.Printf("Starting penrose tiling.\n")
	ts := createSeedTriangles()
	for i := 0; i < 10; i++ {
		ts = subdivideTriangles(ts)
	}
	sort.Sort(ByBC(ts))
	rs := mergeRhombi(ts)
	length := math.Sqrt(math.Pow(real(ts[0].A)-real(ts[0].B), 2) +
		math.Pow(imag(ts[0].A)-imag(ts[0].B), 2))
	scaleFactor := scaleToRhombusSide(length, RHOMBUS_SIDE_IN_MM/MM_PER_PT)
	redInsetRatio, blueInsetRatio := redBlueInsetRatios(RHOMBUS_SIDE_IN_MM, RHOMBUS_INSET_IN_MM)
	drawInsetRhombi("penrose.ps", rs, scaleFactor, redInsetRatio, blueInsetRatio)
	fmt.Printf("Ending penrose tiling.\n")
}
