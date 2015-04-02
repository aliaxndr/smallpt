package main

import (
	"flag"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"runtime"
	"strconv"
	"time"
)

type Vec struct {
	x, y, z float64
}

var NULLVEC *Vec = &Vec{0.0, 0.0, 0.0}

func (this *Vec) Zero() {
	this.x = 0.0
	this.y = 0.0
	this.z = 0.0
}

func (this *Vec) String() string {
	return fmt.Sprintf("{%f, %f, %f}", this.x, this.y, this.z)
}

func (this *Vec) Add(v *Vec) *Vec {
	this.x += v.x
	this.y += v.y
	this.z += v.z
	return this
}

func (this *Vec) Sub(v *Vec) *Vec {
	this.x -= v.x
	this.y -= v.y
	this.z -= v.z
	return this
}

func (this *Vec) Mul(v *Vec) *Vec {
	this.x *= v.x
	this.y *= v.y
	this.z *= v.z
	return this
}

func (this *Vec) SMul(f float64) *Vec {
	this.x *= f
	this.y *= f
	this.z *= f
	return this
}

func (this *Vec) Norm() *Vec {
	l := math.Sqrt(this.x*this.x + this.y*this.y + this.z*this.z)
	this.x /= l
	this.y /= l
	this.z /= l
	return this
}

func (this *Vec) Dot(v *Vec) float64 {
	return this.x*v.x + this.y*v.y + this.z*v.z
}

func (this *Vec) Cross(v *Vec) *Vec {
	var x, y float64
	x = this.y*v.z - this.z*v.y
	y = this.z*v.x - this.x*v.z
	this.z = this.x*v.y - this.y*v.x
	this.x = x
	this.y = y
	return this
}

func (this *Vec) Clone() *Vec {
	return &Vec{this.x, this.y, this.z}
}

func (this *Vec) Clamp() *Vec {
	if this.x > 1.0 {
		this.x = 1.0
	} else if this.x < 0.0 {
		this.x = 0.0
	}
	if this.y > 1.0 {
		this.y = 1.0
	} else if this.y < 0.0 {
		this.y = 0.0
	}
	if this.z > 1.0 {
		this.z = 1.0
	} else if this.z < 0.0 {
		this.z = 0.0
	}
	return this
}

func Dot(v1, v2 Vec) float64 {
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z
}

func Add(v1, v2 Vec) (retv Vec) {
	retv.x = v1.x + v2.x
	retv.y = v1.y + v2.y
	retv.z = v1.z + v2.z
	return
}

func Sub(v1, v2 Vec) (retv Vec) {
	retv.x = v1.x - v2.x
	retv.y = v1.y - v2.y
	retv.z = v1.z - v2.z
	return
}

func Mul(v1, v2 Vec) (retv Vec) {
	retv.x = v1.x * v2.x
	retv.y = v1.y * v2.y
	retv.z = v1.z * v2.z
	return
}

func SMul(v Vec, f float64) (retv Vec) {
	retv.x = v.x * f
	retv.y = v.y * f
	retv.z = v.z * f
	return
}

func Norm(v Vec) Vec {
	return SMul(v, 1.0/math.Sqrt(v.x*v.x+v.y*v.y+v.z*v.z))
}

func Cross(v1, v2 Vec) (retv Vec) {
	retv.x = v1.y*v2.z - v1.z*v2.y
	retv.y = v1.z*v2.x - v1.x*v2.z
	retv.z = v1.x*v2.y - v1.y*v2.x
	return
}

type Ray struct {
	Origin, Direction Vec
}

func (this *Ray) String() string {
	return fmt.Sprintf("{{%s}, {%s}}", this.Origin, this.Direction)
}

const (
	DIFF = iota
	SPEC
	REFR
)

type Sphere struct {
	Radius                    float64
	Position, Emission, Color *Vec
	Reflection                int
}

func (this *Sphere) Intersect(ray *Ray) float64 {
	//op := this.Position.Clone().Sub(ray.Origin)
	op := Sub(*this.Position, ray.Origin)
	var eps float64 = 1e-4
	var b float64 = Dot(op, ray.Direction)
	var det float64 = b*b - Dot(op, op) + this.Radius*this.Radius
	if det < 0.0 {
		return 0.0
	} else {
		det = math.Sqrt(det)
	}
	if (b - det) > eps {
		return b - det
	} else if (b + det) > eps {
		return b + det
	}
	return 0.0
}

var spheres []*Sphere = []*Sphere{
	&Sphere{1e5, &Vec{1e5 + 1, 40.8, 81.6}, &Vec{0.0, 0.0, 0.0}, &Vec{.75, .25, .25}, DIFF},   //Left
	&Sphere{1e5, &Vec{-1e5 + 99, 40.8, 81.6}, &Vec{0.0, 0.0, 0.0}, &Vec{.25, .25, .75}, DIFF}, //Rght
	&Sphere{1e5, &Vec{50, 40.8, 1e5}, &Vec{0.0, 0.0, 0.0}, &Vec{.75, .75, .75}, DIFF},         //Back
	&Sphere{1e5, &Vec{50, 40.8, -1e5 + 170}, &Vec{0.0, 0.0, 0.0}, &Vec{0.0, 0.0, 0.0}, DIFF},  //Frnt
	&Sphere{1e5, &Vec{50, 1e5, 81.6}, &Vec{0.0, 0.0, 0.0}, &Vec{.75, .75, .75}, DIFF},         //Botm
	&Sphere{1e5, &Vec{50, -1e5 + 81.6, 81.6}, &Vec{0.0, 0.0, 0.0}, &Vec{.75, .75, .75}, DIFF}, //Top
	&Sphere{16.5, &Vec{27, 16.5, 47}, &Vec{0.0, 0.0, 0.0}, &Vec{.999, .999, .999}, SPEC},      //Mirr
	&Sphere{16.5, &Vec{73, 16.5, 78}, &Vec{0.0, 0.0, 0.0}, &Vec{.999, .999, .999}, REFR},      //Glas
	&Sphere{600, &Vec{50, 681.6 - .27, 81.6}, &Vec{12, 12, 12}, &Vec{0.0, 0.0, 0.0}, DIFF},    //Lite
}

func Clamp(x float64) float64 {
	if x > 1.0 {
		return 1.0
	} else if x < 0.0 {
		return 0.0
	}
	return x
}

func ToByte(x float64) byte {
	return byte(math.Pow(Clamp(x), 1/2.2)*255 + 0.5)
}

func Intersect(ray *Ray, t *float64, id *int) bool {
	var d, inf float64
	inf = 1e50
	*t = inf
	for i, s := range spheres {
		d = s.Intersect(ray)
		if d != 0.0 && d < *t {
			*t = d
			*id = i
		}
	}
	return *t < inf
}

func Radiance(ray *Ray, depth int, rnd *rand.Rand) Vec {

	var t, p float64
	var id int
	if !Intersect(ray, &t, &id) {
		return *NULLVEC
	}
	obj := spheres[id]
	x := Add(ray.Origin, SMul(ray.Direction, t))
	n := Norm(Sub(x, *obj.Position))
	var nl Vec
	if Dot(n, ray.Direction) < 0.0 {
		nl = n
	} else {
		nl = SMul(n, -1.0)
	}
	f := *obj.Color
	if f.x > f.y && f.x > f.z {
		p = f.x
	} else if f.y > f.z {
		p = f.y
	} else {
		p = f.z
	}
	if depth > 4 {
		if rnd.Float64() < p {
			f = SMul(f, 1.0/p)
		} else {
			return *obj.Emission
		}
	}
	if obj.Reflection == DIFF {
		var r1, r2, r2s float64
		r1 = 2 * math.Pi * rnd.Float64()
		r2 = rnd.Float64()
		r2s = math.Sqrt(r2)
		var w, u, v, d Vec
		w = nl
		if math.Abs(w.x) > 0.1 {
			u = Vec{0.0, 1.0, 0.0}
		} else {
			u = Vec{1.0, 0.0, 0.0}
		}
		v = Cross(w, u)
		d = Norm(Add(Add(SMul(u, math.Cos(r1)*r2s), SMul(v, math.Sin(r1)*r2s)), SMul(w, math.Sqrt(1-r2))))
		return Add(*obj.Emission, Mul(f, Radiance(&Ray{x, d}, depth+1, rnd)))
	} else if obj.Reflection == SPEC {
		return Add(*obj.Emission, Mul(f, Radiance(&Ray{x, Sub(ray.Direction, SMul(n, 2*Dot(n, ray.Direction)))}, depth+1, rnd)))
	}
	var reflRay *Ray = &Ray{x, Sub(ray.Direction, SMul(n, 2*Dot(n, ray.Direction)))}
	var into bool = Dot(n, nl) > 0
	var nc, nt float64 = 1, 1.5
	var nnt, ddn, cos2t float64
	if into {
		nnt = nc / nt
	} else {
		nnt = nt / nc
	}
	ddn = Dot(ray.Direction, nl)
	cos2t = 1 - nnt*nnt*(1-ddn*ddn)
	if cos2t < 0 {
		return Add(*obj.Emission, Mul(f, Radiance(reflRay, depth+1, rnd)))
	}
	var tdir Vec = SMul(ray.Direction, nnt)
	if into {
		tdir = Norm(Sub(tdir, SMul(n, ddn*nnt+math.Sqrt(cos2t))))
	} else {
		tdir = Norm(Sub(tdir, SMul(n, -(ddn*nnt+math.Sqrt(cos2t)))))
	}
	var a, b, c, R0, Re, Tr, P, RP, TP float64
	a, b = nt-nc, nt+nc
	R0 = a * a / (b * b)
	if into {
		c = 1 + ddn
	} else {
		c = 1 - Dot(tdir, n)
	}
	Re = R0 + (1-R0)*c*c*c*c*c
	Tr = 1 - Re
	P = .25 + .5*Re
	RP = Re / P
	TP = Tr / (1 - P)
	if depth > 1 {
		if rnd.Float64() < P {
			return Add(*obj.Emission, Mul(f, SMul(Radiance(reflRay, depth+1, rnd), RP)))
		} else {
			return Add(*obj.Emission, Mul(f, SMul(Radiance(&Ray{x, tdir}, depth+1, rnd), TP)))
		}
	}
	return Add(*obj.Emission, Mul(f, SMul(Add(Radiance(reflRay, depth+1, rnd),
		Radiance(&Ray{x, tdir}, depth+1, rnd)), Tr)))
}

var w, h, samps int = 1024, 768, 1
var cam *Ray
var colors []Vec

func main() {
	flag.Parse()
	var cx, cy Vec
	cam = &Ray{Vec{50, 52, 295.6}, Norm(Vec{0, -0.042612, -1})}
	colors = make([]Vec, h*w)
	if flag.NArg() > 0 {
		samps, _ = strconv.Atoi(flag.Arg(0))
		samps /= 4
	}

	start := time.Now()

	cx = Vec{float64(w) * 0.5135 / float64(h), 0.0, 0.0}
	cy = SMul(Norm(Cross(cx, cam.Direction)), 0.5135)

	fsamps := 1.0 / float64(samps)
	threads_num := runtime.NumCPU()
	runtime.GOMAXPROCS(threads_num)
	ch := make(chan int, threads_num)
	for i := 0; i < threads_num; i++ {
		go func(i int) {
			var (
				r1, r2    float64
				dx, dy    float64
				radiance  Vec
				direction Vec
			)
			rnd := rand.New(rand.NewSource(123 + int64(i)*13))
			for y := i; y < h; y += threads_num {
				for x := 0; x < w; x++ {
					for sy, i := 0, (h-y-1)*w+x; sy < 2; sy++ {
						for sx := 0; sx < 2; sx++ {
							radiance.x = 0
							radiance.y = 0
							radiance.z = 0
							for s := 0; s < samps; s++ {
								r1, r2 = 2*rnd.Float64(), 2*rnd.Float64()
								if r1 < 1 {
									dx = math.Sqrt(r1) - 1
								} else {
									dx = 1 - math.Sqrt(2-r1)
								}
								if r2 < 1 {
									dy = math.Sqrt(r2) - 1
								} else {
									dy = 1 - math.Sqrt(2-r2)
								}
								direction = Add(Add(SMul(cx, ((float64(sx)*.5+dx)/2+float64(x))/float64(w)-.5),
									SMul(cy, ((float64(sy)+.5+dy)/2+float64(y))/float64(h)-.5)), cam.Direction)
								radiance = Add(radiance, SMul(Radiance(&Ray{Add(cam.Origin, SMul(direction, 140.0)), Norm(direction)}, 0, rnd), fsamps))
							}
							colors[i] = Add(colors[i], SMul(radiance, 0.25))
						}
					}
				}
				ch <- 1
			}
		}(i)
	}
	for i := 0; i < h; i++ {
		<-ch
		fmt.Printf("\rRendering %d %%  ", 100*i/(h-1))
	}

	fmt.Println()
	log.Printf("\nRender took %s", time.Since(start))

	f, _ := os.OpenFile("image.ppm", os.O_CREATE|os.O_WRONLY|os.O_TRUNC, 0666)
	f.WriteString(fmt.Sprintf("P3\n%d %d\n%d\n", w, h, 255))
	for _, color := range colors {
		f.WriteString(fmt.Sprintf("%d %d %d ", ToByte(color.x), ToByte(color.y), ToByte(color.z)))
	}
	f.Close()
}
