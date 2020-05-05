package main

import (
	"errors"
	"fmt"
	"log"
	"math"
	"myReader"
)

func main() {
	n, e, A, b, xFirst, err := myReader.ReadSystem("src.txt")
	fmt.Printf("Параметр диагонального преобладание матрицы A: %g\n", dp(A, n))
	if err != nil {
		log.Panic(err)
	}
	kJ, kS := predict(A, b, n, e, xFirst)
	if kJ > 00 {
		fmt.Printf("1. Метод Якоби\nОценочное число итераций: %d\n", kJ)
	} else {
		fmt.Printf("1. Метод Якоби\nОценочное число итераций отсутствует\n")
	}
	x, count, err := Jacobi(A, b, n, e, xFirst, math.MaxInt32)
	if err != nil{
		fmt.Println(err)
	} else {
		fmt.Println("Рещение:")
		fmt.Println(x)
		fmt.Printf("Число итераций: %d\n", count)
	}
	if kS > 0 {
		fmt.Printf("2. Метод Зейделя\nОценочное число итераций: %d\n", kS)
	} else {
		fmt.Printf("2. Метод Зейделя\nОценочное число итераций отсутствует\n")
	}
	x, count, err = Seidel(A, b, n, e, xFirst, math.MaxInt32)
	if err != nil{
		fmt.Println(err)
	} else {
		fmt.Println("Рещение:")
		fmt.Println(x)
		fmt.Printf("Число итераций: %d\n", count)
	}
}

//оценка числа итераций
func predict(A [][]float64, b []float64, n int, e float64, x0 []float64)(int, int){
	var normB float64 = 0
	for i := 0; i < n; i++ {
		var sum float64 = 0
		for j := 0; j < n; j++{
			if i != j {
				sum += math.Abs(A[i][j] / A[i][i])
			}
		}
		if sum > normB {
			normB = sum
		}
	}
	x1J := make([]float64, n)
	x1S := make([]float64, n)
	for i := 0; i < n; i++ {
		var sum float64 = 0
		for j := 0; j < n; j++{
			if i != j {
				sum += (A[i][j] / A[i][i]) * x0[j]
			}
		}
		x1J[i] = b[i] / A[i][i] - sum
	}
	for i := 0; i < n; i++ {
		var sum float64 = 0
		for j := 0; j < i; j++{
			sum += (A[i][j] / A[i][i]) * x1S[j]
		}
		for j := i + 1; j < n; j++{
			sum += (A[i][j] / A[i][i]) * x0[j]
		}
		x1S[i] = b[i] / A[i][i] - sum
	}
	_, normJ := crt1(x1J, x0, e)
	_, normS := crt1(x1S, x0, e)
	//j := math.Log((e * (1 - normB) / normJ)) / math.Log(normB)
	//fmt.Println(j)
	kJ := int(math.Log((e * (1 - normB) / normJ)) / math.Log(normB)) + 1
	kS := int(math.Log((e * (1 - normB) / normS)) / math.Log(normB)) + 1
	return kJ, kS
}

//метод Якоби
//вход: матрица системы, столбец свободных членов, кол-во неизвестных, точность, начальное приблежение, максимальное число итераций
//выход: вектор x, число итераций, ошибка
func Jacobi(A [][]float64, b []float64, n int, e float64, x0 []float64, maxCount int) ([]float64, int, error){
	// запуск итерационного процесса решения
	xk := make([]float64, n)
	for i := 0; i < n; i++{
		xk[i] = x0[i]
	}
	xk1 := make([]float64, n)
	var norm, oldNorm float64
	var count int //кол-во числа итераций
	isEnd := false
	for count = 0; !isEnd && count < maxCount; count++{
		for i := 0; i < n; i++{
			var sum float64 = 0
			for j := 0; j < n; j++{
				if i != j{
					sum += (A[i][j]/A[i][i]) * xk[j]
				}
			}
			xk1[i] = -sum + b[i] / A[i][i]
		}

		isEnd, norm = crt1(xk1, xk, e) //если выполнилось условие завершения
		if norm >= oldNorm && count != 0{ //если метод якоби не сходится
			return nil, -1, errors.New("Метод Якоби расходится")
		}
		// иначе свдинуться по приближениям
		a := xk
		xk = xk1
		xk1 = a
		oldNorm = norm
		//fmt.Printf("%d:", count)
		//fmt.Println(xk1)
	}
	return xk1, count, nil
}

//метод Зейделя
//вход: матрица системы, столбец свободных членов, кол-во неизвестных, точность, начальное приблежение, максимальное число итераций
//выход: вектор x, число итераций, ошибка
func Seidel(A [][]float64, b []float64, n int, e float64, x0 []float64, maxCount int) ([]float64, int, error){
	//запуск итерационного процесса решения
	xk := make([]float64, n)
	for i := 0; i < n; i++{
		xk[i] = x0[i]
	}
	xk1 := make([]float64, n)
	var norm, oldNorm float64
	var count int //кол-во числа итераций
	isEnd := false
	for count = 0; !isEnd && count < maxCount; count++{
		for i := 0; i < n; i++{
			var sum float64 = 0
			for j := 0; j < i; j++{
				sum += (A[i][j]/A[i][i]) * xk1[j]
			}
			for j := i + 1; j < n; j++{
				sum += (A[i][j]/A[i][i]) * xk[j]
			}
			xk1[i] = -sum + b[i] / A[i][i]
		}

		isEnd, norm = crt1(xk1, xk, e) //если выполнилось условие завершения
		if norm >= oldNorm && count != 0{ //если метод якоби не сходится
			return nil, -1, errors.New("Метод Зейделя расходится")
		}
		// иначе свдинуться по приближениям
		a := xk
		xk = xk1
		xk1 = a
		oldNorm = norm
	}
	return xk1, count, nil
}
//критериии завершения итераций
//1.
func crt1(xk1 []float64, xk []float64, e float64) (bool, float64){
	var norm float64 = -1
	for i, _ := range xk {
		if math.Abs(xk1[i] - xk[i]) > norm {
			norm = math.Abs(xk1[i] - xk[i])
		}
	}
	return norm < e, norm
}

//вычисления диагонального преобладания матрицы
func dp(A [][]float64, n int) float64 {
	var res float64 = 0
	for i := 0; i < n; i++{
		var sum float64 = 0
		for j := 0; j < n; j++{
			if i != j {
				sum += math.Abs(A[i][j] / A[i][i])
			}
		}
		if sum > res{
			res = sum
		}
	}
	return res
}

