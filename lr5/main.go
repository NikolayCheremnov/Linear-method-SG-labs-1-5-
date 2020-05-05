package main

import (
	"fmt"
	"math"
	"myReader"
)

func main() {
	n, A, e, err := myReader.ReadM5("big.txt")
	if err != nil {
		panic(err)
	}
	oldA := make([][]float64, n)
	for i := 0; i < n; i++{
		oldA[i] = make([]float64, n)
		copy(oldA[i], A[i])
	}
	fmt.Println("Исходные данные:\n", n, A, e)
	eigenValues, eigenVectors, rotates := rotate(n, A, e, 1000000)
	fmt.Println("Повортов сделано: ", rotates)
	fmt.Println("Собственное число; собственный вектор; невязка:")
	for i, elem := range eigenVectors {
		fmt.Println(eigenValues[i], ";", elem, ";", r(n, oldA, eigenValues[i], elem))
	}
}

//процедура вычисления невязки
func r(n int, A [][]float64, l float64, x []float64) []float64 {
	res := make([]float64, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			res[i] += A[i][j] * x[j]
		}
	}
	for i := 0; i < n; i++{
		res[i] -= l * x[i]
	}
	return res
}

//процедура метода вращения
//вход: размерность матрицы, матрица, точность
//выход:
func rotate(n int, A [][]float64, e float64, criticalIteration int) ([]float64, [][]float64, int){
	r := make([]float64, n) 	//суммы недиагональных элементов каждой строки
	for i := 0; i < n; i++ {	//заполнить суммы
		for j := 0; j < n; j++ {
			if i != j {
				r[i] += A[i][j] * A[i][j]
			}
		}
	}
	D := make([][] float64, n) //создать изначальную матрицу D
	for i := 0; i < n; i++ {
		D[i] = make([]float64, n)
		D[i][i] = 1
	}
	//основной итерационный цикл
	rotatesCount := 0
	for k, l := getkl(n, A, r, e); k != -1 && l != -1 && rotatesCount < criticalIteration; k, l = getkl(n, A, r, e) {
		U := getU(n, A, k, l) //получить очередную матрицу поворота
		D = du(n, D, U, k, l) //обновить матрицу D
		A = rotation(n, A, k, l) //повернуть матрицу
		r = rUpdate(n, A, r, k, l) //обновить суммы
		rotatesCount++
	}
	eigenValues, eigenVectors := recoverResults(n, A, D)
	return eigenValues, eigenVectors, rotatesCount
}

func recoverResults(n int, A [][]float64, D [][]float64)([]float64, [][]float64){
	eigenValues := make([]float64, n)
	for i := 0; i < n; i++{
		eigenValues[i] = A[i][i] //проинициализировать собственные числа
	}
	eigenVectors := make([][]float64, n)
	for i := 0; i < n; i++{
		eigenVectors[i] = make([]float64, n)	//очередной собственный вектор для iго собственного числа
		for j := 0; j < n; j++{
			eigenVectors[i][j] = D[j][i]
		}
	}
	return eigenValues, eigenVectors
}

//процедура корректирования сумм r
//вход: размерность матрицы, очередная матрица A, k, l, вектор сумм
//выход: обновленный вектор сумм
func rUpdate(n int, A [][]float64, r []float64, k int, l int) []float64{
	r[k] = 0
	r[l] = 0
	for j := 0; j < n; j++ {
		if j != k {
			r[k] += A[k][j] * A[k][j]
		}
	}
	for j := 0; j < n; j++ {
		if j != l {
			r[l] += A[l][j] * A[l][j]
		}
	}
	return r
}

//C=AU
func rotation(n int, A [][]float64, k int, l int)(C [][]float64){
	a, b := getAlphBeta(A, k, l)
	kCol := make([]float64, n)
	lCol := make([]float64, n)
	//сделали A = AU
	for i := 0; i < n; i++ {
		kCol[i] = A[i][k] * a + A[i][l] * b
		lCol[i] =  A[i][l] * a - A[i][k] * b
	}
	for i := 0; i < n; i++{
		A[i][k] = kCol[i]
		A[i][l] = lCol[i]
	}
	//сделаем A = UtA
	for i := 0; i < n; i++ {
		kCol[i] = A[k][i] * a + A[l][i] * b
		lCol[i] =  A[l][i] * a - A[k][i] * b
	}
	for i := 0; i < n; i++{
		A[k][i] = kCol[i]
		A[l][i] = lCol[i]
	}
	return A
}

//процедра D = DU
//вход: n, матрицы D и U, числа k и l
//выход:
func du(n int, D [][]float64, U [][]float64, k int, l int) [][]float64 {
	kCol := make([]float64, n)
	lCol := make([]float64, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++{
			kCol[i] += D[i][j] * U[j][k]
			lCol[i] += D[i][j] * U[j][l]
		}
	}
	for i := 0; i < n; i++ { //изменить матрицу D
		D[i][k] = kCol[i]
		D[i][l] = lCol[i]
	}
	return D
}

//формирования матрицы U
//вход:
//выход: матрица U
func getU(n int, A [][]float64, k int, l int) [][]float64 {
	U := make([][]float64, n)
	a, b := getAlphBeta(A, k, l)
	for i := 0; i < n; i++{ //заполним матрицу u
		U[i] = make([]float64, n)
		for j := 0; j < n; j++{
			if i == j {	//если диагональный элемент
				if i == k || i == l {
					U[i][j] = a
				} else {
					U[i][j] = 1
				}
			} else { //иначе недиагональный элемент
				if i == k && j == l {
					U[i][j] = -b
				} else if i == l && j == k {
					U[i][j] = b
				} else {
					U[i][j] = 0
				}
			}
		}
	}
	return U
}

//процедура получения a, b
func getAlphBeta(A [][]float64, k int, l int)(ap float64, bt float64){
	if A[k][k] == A[l][l]{
		return math.Sqrt(0.5), math.Sqrt(0.5)
	}

	ny := 2 * A[k][l] / (A[k][k] - A[l][l])
	ap = math.Sqrt(0.5 * (1 + 1/ (math.Sqrt(1 + ny*ny))))
	bt = math.Sqrt(0.5 * (1 - 1/(math.Sqrt(1 + ny*ny))))
	if ny < 0{
		bt *= -1
	}
	return ap, bt
}

//более топорный способ получения kl
func fooGetkl(n int, A [][]float64, e float64) (int, int) {
	k := -1
	l := -1
	for i := 0; i < n; i++{
		for j := 0; j < n; j++ {
			if i != j && (k == -1 || A[i][j] * A[i][j] > A[k][l] * A[k][l]){
				k = i
				l = j
			}
		}
	}
	if math.Abs(A[k][l]) < e {
		return -1, -1
	} else {
		return k, l
	}
}

//процедура выбора элемента akl
//вход: размерность матрицы, матрица, массив сумм r, тончость e
//выход: k, l
func getkl(n int, A [][]float64, r []float64, e float64) (int, int){
	k := -1
	for i := 0; i < n; i++ {
		if k == -1 || r[i] > r[k] {
			k = i
		}
	}
	l := -1
	for j := 0; j < n; j++ { //найти минимальное слагаемое
		if j != k && (l == -1 || A[k][j] * A[k][j] > A[k][l] * A[k][l]) { //если есть более большое слагаемое
			l = j
		}
	}
	if math.Abs(A[k][l]) < e { //если достигнута нужная точность
		return -1, -1
	} else {
		return k, l
	}
}