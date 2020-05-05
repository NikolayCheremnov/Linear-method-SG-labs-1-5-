package main

import (
	"fmt"
	"log"
	"math"
	"os"
)

func main(){
	cmdArgs := os.Args
	if len(cmdArgs) < 2 {
		log.Panic("Укажите файл с исходными данными")
	}
	var resFile *os.File = nil
	if len(cmdArgs) > 2 { //если указан файл для результатов
		resFile, _ = os.Create(cmdArgs[2])
	}
	Ab, n, err := readSourceData(cmdArgs[1]) //считать исходные данные
	if err != nil {
		if resFile != nil {
			fmt.Fprint(resFile, err)
			resFile.Close()
		}
		log.Panic(err)
	}

	x, r, A_, R := systemProcessing(Ab, n) //выполнить обработку системы
	if err != nil {
		if resFile != nil {
			fmt.Fprint(resFile, err)
			resFile.Close()
		}
		log.Panic(err)
	}

	if resFile != nil { //если указан так же файл куда записать результат
		defer resFile.Close()
		fmt.Fprint(resFile, "Решение: ")
		fmt.Fprintln(resFile, x)
		fmt.Fprint(resFile, "Невязка решения: ")
		fmt.Fprintln(resFile, r)
		fmt.Fprintln(resFile, "Обратная матрица:")
		for i := 0; i < n; i++ {
			fmt.Fprintln(resFile, A_[i])
		}
		fmt.Fprintln(resFile, "Невязка обратной матрицы:")
		for i := 0; i < n; i++ {
			fmt.Fprintln(resFile, R[i])
		}
		resFile.Close()
	} else { //иначе результаты в консоль
		fmt.Print("Решение: ")
		fmt.Println(x)
		fmt.Print("Невязка решения: ")
		fmt.Println(r)
		fmt.Println("Обратная матрица:")
		for i := 0; i < n; i++ {
			fmt.Println(A_[i])
		}
		fmt.Println("Невязка обратной матрицы:")
		for i := 0; i < n; i++ {
			fmt.Println(R[i])
		}
	}
}

//процедура обработки системы
//вход:
//выход: решение, невязка решения, обратная матрица, невязка обратной матрицы
func systemProcessing(Ab [][]float64, n int) ([]float64, []float64, [][]float64, [][]float64){
	//1. сформировать из Ab матрицу AbE (для нахождения обратной матрицы)
	AbE := Ab
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			if i == j {
				AbE[i] = append(AbE[i], 1) //1 на главной диагонали подматрицы E
			} else {
				AbE[i] = append(AbE[i], 0) //иначе 0
			}
		}
	}

	//2. получить матрицы S и D для матрицы A
	S, D := getSD(Ab, n)

	//3. решить исходную систему
	x := solve(AbE, n, n, S, D)

	//4. определить вектор невязки для полученного решения
	solutionDelta := make([]float64, n) //вектор невязки
	for i:= 0; i < n; i++{
		var sum float64 = 0
		for j := 0; j < n; j++ {
			sum += AbE[i][j] * x[j]
		}
		solutionDelta[i] = AbE[i][n] - sum
	}

	//5. получить решения являющиеся столбцами обратной матрицы и сформировать обратную матрицу
	A_ := make([][]float64, n) //обратная матрица
	for i := n + 1; i < n + 1 + n; i++ { // цикл по всем столбцам подматрицы E
		column := solve(AbE, n, i, S, D) //решить систему и получить очередной столбец обратной матрицы
		for j := 0; j < n; j++{
			A_[j] = append(A_[j], column[j]) //заполнить столбец обратной матрицы
		}
	}

	//6. определить невязку для обратной матрицы
	R := make([][]float64, n) //произведение прямой и обратной матрицы A, затем вектор невязки
	for i := 0; i < n; i++{ //цикл по строкам
		for j := 0; j < n; j++ { //цикл по столбцам
			R[i] = append(R[i], 0) //добавить элемент AA_[i][j]
			for k := 0; k < n; k++ { //цикл умножения iой строки A на jый столбец A_
				R[i][j] += AbE[i][k] * A_[k][j]
			}
		}
	}
	//преобразовать AA_ к вектору невязки (AA_ = E - AA_)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++{
			if i == j { //если элемент на главной диагонали
				R[i][j] = 1 - R[i][j]
			} else {
				R[i][j] = 0 - R[i][j]
			}
		}
	}
	return x, solutionDelta, A_, R
}

//процедура получения корней системы линейных уравнений по матрицам S, D
//вход: расширенная матрица системы, кол-во неизвестных, индекс столбца вектора правой части системы, матрицы S и D
//выход:
func solve(AbE [][]float64, n int, bIndex int, S [][]float64, D []float64) []float64{
	//1 решить уравнение StZ = B
	z := make([]float64, n) // вектор z
	for i := 0; i < n; i++ {
		var sum float64 = 0
		for j := 0; j < i; j++ {
			sum += S[j][i] * z[j]
		}
		z[i] = (AbE[i][bIndex] - sum) / S[i][i]
	}

	//2 решить уравнение Dy = z
	y := z //используем массив z для решений
	for i := 0; i < n; i++ {
		y[i] = z[i] / D[i]
	}

	//3. решить Sx = y
	x := y //используем всё тот же массив для сохранения решений
	for i := n - 1; i >= 0; i-- {
		var sum float64 = 0
		for j := i + 1; j < n; j++{
			sum += S[i][j] * x[j]
		}
		x[i] = (y[i] - sum) / S[i][i]
	}

	return x
}

//процедура получения матрицы S, D
//вход: исходная матрица А, кол-во неизвестынх (размерность A)
//выход: матрицы S, D
func getSD(A [][]float64, n int)  ([][]float64, []float64) {
	var d []float64 //диагональная матрица d (в виде массива)
	s := make([][]float64, n) //матрица s
	for i := 0; i < n; i++ { // цикл по размерности матрицы
		var sum float64 = 0
		for k := 0; k < i; k++ { //получить сумму
			sum += s[k][i]*s[k][i]*d[k]
		}
		d = append(d, sign(A[i][i] - sum)) //получить dii
		s[i] = make([]float64, n) //выделить память под строку s
		s[i][i] = math.Sqrt(math.Abs(A[i][i] - sum)) //получить sii
		for j := 0; j < n; j++ { //окончательно найти строку s[i]
			if i < j { //если элемент выше главной диагонали
				sum = 0
				for k := 0; k < i; k++ {
					sum += s[k][i] * s[k][j] * d[k]
				}
				s[i][j] = (A[i][j] - sum) / (s[i][i] * d[i])
			}
		}
	}

	return s, d //вернуть полученные матрицы
}

//вспомогательная процедура знака
func sign(val float64) float64 {
	if val < 0 {
		return -1
	} else{
		return 1
	}
}

//процедура чтения исходных данных:
//считывает расширенную матрицу уравнения из заданного файла
func readSourceData(fileName string) ([][]float64, int, error) {
	file, err := os.Open(fileName) // открыть файл с данными
	defer file.Close() // закрыть по завершению
	if err != nil { // если ошибка открытия, то завершения с ошибкой
		return nil, -1, err
	}
	var n int // кол-во неизвестных
	_, err = fmt.Fscanf(file, "%d\n", &n) //считать кол-во неихвестныъ
	if err != nil { // если ошибка, то завершение
		return nil, -1, err
	}
	var Ab [][]float64 = make([][]float64, n) //матрица n x (n + 1) - расширенная вида A|b
	//цикл чтения матрицы из файла
	for i := 0; i < n; i++{
		Ab[i] = make([]float64, n + 1) // выделить память под строку матрицы
		for j := 0; j <= n; j++ { // считтаь элементы строки
			_, err = fmt.Fscanf(file, "%g", &Ab[i][j])
			if err != nil {
				return nil, -1, err
			}
		}
		fmt.Fscanf(file, "\n") //пропустить newline
	}
	return Ab, n, nil // вернуть результат
}