package main

import (
	"errors"
	"fmt"
	"log"
	"math"
	"os"
)

//программа
//параметры комадной строки: имя файла с данными
func main() {
	cmdArgs := os.Args
	if len(cmdArgs) < 2 {
		log.Panic("Укажите файл с исходными данными")
	}
	var resFile *os.File = nil
	if len(cmdArgs) > 2 && cmdArgs[2] != "-wm" { //если указан файл для результатов
		resFile, _ = os.Create(cmdArgs[2])
	}
	wm := false //with major
	if (len(cmdArgs) > 2 && cmdArgs[2] == "-wm") ||(len(cmdArgs) > 3 && cmdArgs[3] == "-wm") {
		wm = true
	}
	Ab, n, err := readSourceData(cmdArgs[1]) //считать исходные данные
	if err != nil {
		if resFile != nil {
			fmt.Fprint(resFile, err)
			resFile.Close()
		}
		log.Panic(err)
	}

	x, r, detA, A_, R, err := systemProcessing(Ab, n, wm) //выполнить обработку системы
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
		fmt.Fprintf(resFile, "det = %g\n", detA)
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
		fmt.Printf("det = %g\n", detA)
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

//процедура обработки системы: получает решение системы, велечину невязки, определитель матрицы системы и обратную матрицу системы
//входные параметры: исходная матрица уравнения Ab, кол-во неизвестных, признак использования главного элемента в алгоритме приведения к треугольному виду
//результат: решение системы, вектор невязки для решения системы, определитель, обратная матрица, невязка для обратной матрицы, информация об ошибке
func systemProcessing(A [][]float64, n int, withMajor bool) ([]float64, []float64, float64, [][]float64, [][]float64, error){
	//скопировать матрицу A и вектор b для вычисления невязки обратной матрицы ниже
	AbE := A //сформируем матрицу вида A|b|E - расширенную марице уравнения
	A = make([][]float64, n)
	for i := 0; i < n; i++{
		A[i] = make([]float64, n)
		copy(A[i], AbE[i][:n])
	}
	b := make([]float64, n)
	for i := 0; i < n; i++{
		b[i] = AbE[i][n]
	}
	for i := 0; i < n; i++{ //цикл по строкам
		for j := n + 1; j <= n * 2; j++{ //цикл по столбцам
			if i == j - (n + 1) { //если элемент на главной диагонали подматрицы E
				AbE[i] = append(AbE[i], 1) //то записать 1
			} else {
				AbE[i] = append(AbE[i], 0) //иначе 0
			}
		}
	}
	AbE, swapStringNumb, err := reductionToTriangular(AbE, n, withMajor) // привести всё это дело к треугольному виду
	if err != nil { //проверка на возникновение ошибок
		return nil, nil, 0, nil, nil, err
	}

	solution, err := solve(AbE, n, n)//далее получить решение для системы уравнений
	if err != nil {
		return nil, nil, 0, nil, nil, err
	}

	//определить вектор невязки для полученного решения
	solutionDelta := make([]float64, n) //вектор невязки
	for i:= 0; i < n; i++{
		var sum float64 = 0
		for j := 0; j < n; j++ {
			sum += A[i][j] * solution[j]
		}
		solutionDelta[i] = b[i] - sum
	}

	//далее получить решения являющиеся столбцами обратной матрицы и сформировать обратную матрицу
	A_ := make([][]float64, n) //обратная матрица
	for i := n + 1; i < n + 1 + n; i++ { // цикл по всем столбцам подматрицы E
		column, err := solve(AbE, n, i) //решить систему и получить очередной столбец обратной матрицы
		if err != nil {
			return nil, nil, 0, nil, nil, err
		}
		for j := 0; j < n; j++{
			A_[j] = append(A_[j], column[j]) //заполнить столбец обратной матрицы
		}
	}

	//определить невязку для обратной матрицы
	R := make([][]float64, n) //произведение прямой и обратной матрицы A, затем вектор невязки
	for i := 0; i < n; i++{ //цикл по строкам
		for j := 0; j < n; j++ { //цикл по столбцам
			R[i] = append(R[i], 0) //добавить элемент AA_[i][j]
			for k := 0; k < n; k++ { //цикл умножения iой строки A на jый столбец A_
				R[i][j] += A[i][k] * A_[k][j]
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


	detA := det(AbE, n, swapStringNumb) // вычислить определитель матрицы A
	return solution, solutionDelta, detA, A_, R, nil
}


//процедура получения корней системы линейных уравнений по треугольной матрице системы
//входные параметры: расширенная матрица системы, кол-во неизвестных, индекс столбца вектора правой части системы
func solve(AbE [][]float64, n int, bIndex int) ([]float64, error) {
	res := make([]float64, n) // результат
	for i := n - 1; i >= 0; i-- { // цикл подъема по матрице вверх для нахождения корней
		var sum float64 = 0.0; // получить отнимаемую сумму
		for j := i + 1; j < n; j++ {
			sum += AbE[i][j] * res[j]
		}
		res[i] = (AbE[i][bIndex] - sum) / AbE[i][i]
		if math.IsNaN(res[i]) || math.IsInf(res[i], 1) || math.IsInf(res[i], -1) { // проверка на выход за пределы множества машинных чисел
			return nil, errors.New("Одно из чисел за пределами множества машинных чисел, получить решение невозможно. Попробуйте использовать алгоритм с выбором главного элемента.")
		}
	}
	return res, nil
}

//процедура нахождения определителя исходной матрицы A по приведенной треугольной матрице Ab
//входные параметры: приведенная к треугольному виду матрица Ab, число неизвестных, число перестановок строк в процессе приведения
//результат: определитель матрицы
func det(Ab [][]float64, n int, swapStringNumb int) float64 {
	var det float64 = 1
	for i := 0; i < n; i++ {
		det *= Ab[i][i]
	}
	if swapStringNumb % 2 == 1{ //если необходимо изменить знак
		det *= -1
	}
	return det
}
//процедура приведения к треугольному виду матрицы Ab
//треугольный вид имеет подматрица A
//входные параметры: матрица уравнения Ab, кол-во неизвестных n, признак использования главного элемента в алгоритме
//результат: матрица, приведённая к указанному виду, число перестановок строк в процессе алгоритма, сообщение об ошибках
func reductionToTriangular (Ab [][]float64, n int, withMajor bool) ([][]float64, int,  error){
	swapStringNumb := 0 // число перестановок строк в процессе алгоритма
	for i := 0; i < n; i++ { // цикл по столбцам матрицы A
		if withMajor { // если необходимо ставить на диагональ главный элемент
			major := i // индекс первичного значения
			for j := i + 1; j < n; j++ {
				if (math.Abs(Ab[major][i]) == 0.0) || ((math.Abs(Ab[j][i]) < math.Abs(Ab[major][i])) && math.Abs(Ab[j][i]) != 0.0) { // сравнить по модулю
					major = j
				}
			}
			if math.Abs(Ab[major][i]) == 0.0 { // если матрица вырожденая
				return nil, swapStringNumb, errors.New("Вырожденая матрица")
			}
			if i != major { //если выполняется перестановка
				Ab = swapString(Ab, i, major) // установить главный элемент если требуется
				swapStringNumb++ //увеличить число перестановок
			}
		}
		for j := i + 1; j < n; j++ { // цикл по строкам для столбца i
			Cji := Ab[j][i] / Ab[i][i] // строка i умножается на Cji и вычитается из строки j
			if math.IsNaN(Cji) || math.IsInf(Cji, 1) || math.IsInf(Cji, -1) { // проверка на выход за пределы множества машинных чисел
				return nil, swapStringNumb, errors.New("Одно из чисел за пределами множества машинных чисел, получить решение невозможно. Попробуйте использовать алгоритм с выбором главного элемента.")
			}
			for k := i; k <= n * 2; k++ {  // цикл вычитания строки i из строки j
				Ab[j][k] -= Ab[i][k] * Cji                                                       // вычесть очередной элемент строки
			}
		}
	}
	return Ab, swapStringNumb, nil
}

//вспомогательная процедура перестановки строк i и j в матрице m
func swapString(m [][]float64, i int, j int) [][]float64{
	buf := m[i]
	m[i] = m[j]
	m[j] = buf
	return m
}
//процедура чтения исодных данных:
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
