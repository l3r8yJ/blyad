package main

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"sort"
	"time"
)

const N = 50

func main() {
	var M = 10
	var phi = normalDistribution(math.Pi, math.Pi)
	var r = normalDistribution(5, 1)
	var eps = 0.1
	log.Printf("M0 = %d, N = %d\n", M, N)
	sort.Float64s(phi[0:N])
	phi[N] = phi[0] + 2*math.Pi
	r[N] = phi[0]
	h := hCalc(phi)
	toConsole(phi)
	log.Printf("h = %f", h)
	ZV4(M, h, eps, phi, r)
	ZV5(M, h, eps, phi, r)
}

func ZV3(M int, h float64, eps float64, phi [N + 1]float64, r [N + 1]float64) {
	start := time.Now()
	ok := true
	i := 0
	for ok {
		for i < N+1 {
			if f(M, h, phi[i], eps, phi, r) < r[i] {
				eps += h
				M++
				i = 0
			} else {
				i++
			}
		}
		ok = false
	}
	d := maxDiff(M, h, eps, phi, r)
	log.Printf("M3 = %d eps3 = %.5f dif3 = %.5f\n", M, eps, d)
	elapsed := time.Since(start)
	log.Printf("ZV3 took %s", elapsed)
}

func ZV4(M int, h float64, eps float64, phi [N + 1]float64, r [N + 1]float64) {
	start := time.Now()
	ok := true
	i := 0
	for ok {
		for i < N+1 {
			fGov := f(M, h, phi[i], eps, phi, r)
			if fGov < r[i] {
				eps += r[i] - fGov
				M++
				i = 0
			} else {
				i++
			}
		}
		ok = false
	}
	d := maxDiff(M, h, eps, phi, r)
	log.Printf("M4 = %d eps4 = %.5f dif4 = %.5f\n", M, eps, d)
	elapsed := time.Since(start)
	log.Printf("ZV4 took %s", elapsed)
}

func ZV5(M int, h float64, eps float64, phi [N + 1]float64, r [N + 1]float64) {
	start := time.Now()
	ok := true
	i := 0
	for ok {
		for i < N+1 {
			fGov := f(M, h, phi[i], 0, phi, r) + eps
			if fGov < r[i] {
				eps += r[i] - fGov
				M++
				i = 0
			} else {
				i++
			}
		}
		ok = false
	}
	d := maxDiff1(M, h, eps, phi, r)
	log.Printf("M5 = %d eps5 = %.5f dif5 = %.5f\n", M, eps, d)
	elapsed := time.Since(start)
	log.Printf("ZV5 took %s", elapsed)
}

func maxDiff(M int, h float64, eps float64, phi [N + 1]float64, r [N + 1]float64) float64 {
	max := -1.0
	for i := 0; i < N; i++ {
		fGov := f(M, h, phi[i], eps, phi, r) - r[i]
		if fGov > max {
			max = fGov
		}
	}
	return max
}

func maxDiff1(M int, h float64, eps float64, phi [N + 1]float64, r [N + 1]float64) float64 {
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

func toConsole(list [N + 1]float64) {
	for _, value := range list {
		fmt.Printf("%f\n", value)
	}
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
