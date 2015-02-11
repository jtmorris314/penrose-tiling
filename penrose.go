/* Package penrose implements a basic penrose tiling scheme based
   on subdividing two triangle types ("blue" and "red") that are 1/2
   of the fat and thin rhombii tiles, respectively.  After the triangles
   are subdivided, the similarly colored adjoing tiles are joined to
   create the rhombii.

   This go package is inspired from the python code and explanation
   found here:
      http://preshing.com/20110831/penrose-tiling-explained/

   The seed pattern is a wheel of 10 "red" triangles.

   The program currently iterates 8 times (8 generations of subdividing
   triangles.
*/
package main

import (
	"fmt"
	"github.com/ungerik/go-cairo"
	"math"
	"math/cmplx"
)

// GoldenRatio of yore is represented as complex128 number.
var goldenRatio complex128 = complex((1+math.Sqrt(5))/2, 0)

// ColorType represents the type of rhombii ("red" is thin and
// "blue" is thick).
type ColorType int

// Symbols for "red" and "blue" rhombii/triangles
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
// filled and and the rhombii outlined.  The destination file
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

// printTriangles dumps the values associated with every
// triangle in the slice of triangles.
func printTriangles(ts []Triangle) {
	for _, t := range ts {
		fmt.Printf("Color=%d, A=(%f,%f), B=(%f,%f), C=(%f,%f)\n",
			t.color, real(t.A), imag(t.A), real(t.B), imag(t.B), real(t.C), imag(t.C))
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

// main simply creates the seed triangles and then runs
// 8 generations of triangle subdivision.  The result
// are drawn to penrose.ps.
func main() {
	fmt.Printf("Starting penrose tiling.\n")
	ts := createSeedTriangles()
	for i := 0; i < 8; i++ {
		ts = subdivideTriangles(ts)
	}
	drawTriangles("penrose.ps", ts)
	fmt.Printf("Ending penrose tiling.\n")
}
