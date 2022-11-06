package main

import (
	"encoding/json"
	"fmt"
	"io/ioutil"
	"log"
	"math"
	"math/rand"
	"os"
	"sort"
	"time"
)

const N = 50

var delta = 0.1

type Case struct {
	PHI   [N + 1]float64
	R     [N + 1]float64
	XAXIS []float64
	YAXIS []float64
	EPS   float64
	DIFF  float64
	DELTA float64
	H     float64
	M     int
	N     int
}

func main() {
	var M = 10
	var phi = normalDistribution(math.Pi, 2*math.Pi)
	var r = normalDistribution(5, 1)
	var eps = 0.1
	sort.Float64s(phi[0:N])
	phi[N] = phi[0] + 2*math.Pi
	r[N] = phi[0]
	h := hCalc(phi)
	log.Printf("h = %f", h)
	ZV(M, h, eps, phi, r)
}

func xAxis(h float64) []float64 {
	var res = make([]float64, 0)
	for i := 0.0; i < 2*math.Pi; i += h / 10 {
		res = append(res, i)
	}
	return res
}

func epsCalc(i int, k int, eps float64, fGov float64, r [N + 1]float64) float64 {
	switch k {
	case 1:
		return eps + r[i] - fGov
	case 2:
		return eps + delta
	default:
		return 0.0
	}
}

func ZV(M int, h float64, eps float64, phi [N + 1]float64, r [N + 1]float64) {
	var tEps = eps
	var tM = M
	for k := 1; k <= 2; k++ {
		for j := 1; j <= 2; j++ {
			M = tM
			eps = tEps
			start := time.Now()
			ok := true
			i := 0
			xs := xAxis(h)
			var ys []float64
			for ok {
				for i < N {
					fGov := f(M, h, phi[i], 0, phi, r) + eps
					if fGov < r[i] {
						eps = epsCalc(i, k, eps, fGov, r)
						if j == 1 {
							M++
						}
						i = 0
					} else {
						i++
					}
				}
				ok = false
				for i := 0; i < len(xs); i++ {
					ys = append(ys, f(M, h, xs[i], 0, phi, r)+eps)
				}
			}
			d := maxDiff(M, h, eps, phi, r)
			log.Printf(
				"M_%d_%d = %d eps_%d_%d = %.5f dif_%d_%d = %.5f \n",
				k, j, M,
				k, j, eps,
				k, j, d,
			)
			elapsed := time.Since(start)
			log.Printf("ZV_%d_%d took %s", k, j, elapsed)
			c := Case{
				EPS:   eps,
				H:     h,
				DIFF:  d,
				DELTA: delta,
				PHI:   phi,
				R:     r,
				M:     M,
				N:     N,
				XAXIS: xs,
				YAXIS: ys,
			}
			res, err := json.Marshal(c)
			if err != nil {
				log.Println(err)
			}
			ex := ioutil.WriteFile(fmt.Sprintf("json/case_%d_%d.json", k, j), res, os.ModePerm)
			if ex != nil {
				log.Println(ex)
			}
		}
	}
}

func maxDiff(M int, h float64, eps float64, phi [N + 1]float64, r [N + 1]float64) float64 {
	max := -1.0
	for i := 0; i < N; i++ {
		fGov := f(M, h, phi[i], 0, phi, r) - r[i] + eps
		if fGov > max {
			max = fGov
		}
	}
	return max
}

func hCalc(phi [N + 1]float64) float64 {
	var h = phi[1] - phi[0]
	for i := 1; i < N; i++ {
		t := phi[i+1] - phi[i]
		if t < h {
			h = t
		}
	}
	return h
}

func a0(h float64, eps float64, phi [N + 1]float64, r [N + 1]float64) float64 {
	var sum = 0.0
	for i := 0; i < N; i++ {
		sum += (r[i]+eps)*2*h/3 + C(r[i])*(phi[i+1]-phi[i]-2*h/3)
	}
	sum /= math.Pi
	return sum
}

func an(h float64, n int, eps float64, phi [N + 1]float64, r [N + 1]float64) float64 {
	var sum = 0.0
	for i := 0; i < N; i++ {
		sum += (r[i] + eps) * (math.Sin((phi[i]+h/3)*float64(n)) - math.Sin((phi[i]-h/3)*float64(n)))
		sum += C(r[i]) * (math.Sin((phi[i+1]-h/3)*float64(n)) - math.Sin((phi[i]+h/3)*float64(n)))
	}
	sum /= math.Pi * float64(n)
	return sum
}

func bn(h float64, n int, eps float64, phi [N + 1]float64, r [N + 1]float64) float64 {
	var sum = 0.0
	for i := 0; i < N; i++ {
		sum += (r[i] + eps) * (math.Cos((phi[i]-h/3)*float64(n)) - math.Cos((phi[i]+h/3)*float64(n)))
		sum += C(r[i]) * (math.Cos((phi[i]+h/3)*float64(n)) - math.Cos((phi[i]-h/3)*float64(n)))
	}
	sum /= math.Pi * float64(n)
	return sum
}

func f(M int, h float64, x float64, eps float64, phi [N + 1]float64, r [N + 1]float64) float64 {
	var sum = 0.0
	sum += a0(h, eps, phi, r) / 2
	for n := 1; n <= M; n++ {
		sum += an(h, n, eps, phi, r)*math.Cos(float64(n)*x) + bn(h, n, eps, phi, r)*math.Sin(float64(n)*x)
	}
	return sum
}

func C(r float64) float64 {
	return r / 10
}

func normalDistribution(mu float64, sigma float64) [N + 1]float64 {
	var accum = [N + 1]float64{}
	for i := 0; i < N; i++ {
		accum[i] = math.Mod(math.Abs(normalInverse(mu, sigma)), 2*math.Pi)
	}
	return accum
}

func normalInverse(mu float64, sigma float64) float64 {
	rand.Seed(time.Now().UnixNano())
	return rand.Float64()*sigma + mu
}
