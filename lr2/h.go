package main

import (
	"fmt"
	"math"
	"os"
)

func main(){

	//цикл исследования по n
	f, _ := os.Create("h.dat")
	defer f.Close()
	for n := 20; n <= 20; n++ {
		fmt.Fprintf(f, "n = %d:\n", n)
		//1. сформировать матрицу гильберта
		fmt.Fprintln(f,"Матрица Гильберта:")
		H := make([][]float64, n)
		for i := 0; i < n; i++ {
			for j := 0; j < n; j++ {
				H[i] = append(H[i], 1/float64(i+j+1))
			}
			fmt.Fprintln(f, H[i]) //вывести строку матрицы
		}

		//2. получить H_ - обратную матрицу обработави фиктивную систему:
		for i := 0; i < n; i++ {
			H[i] = append(H[i], 0) //добавить пустой фиктивный столбец b
		}
		_, _, H_, _ := systemProcessing(H, n) //получить обратную матрицу
		fmt.Fprintln(f, "Обратная матрица Гильберта:")
		for i := 0; i < n; i++ {
			fmt.Fprintln(f, H_[i]) //вывод обратной матрицы
		}

		//3. рассчитать нормы матриц и получить число обусловленности
		var HNorm float64 = 0 //норма H (минимальное значение 0)
		var H_Norm float64 = 0
		for i := 0; i < n; i++ {
			var Hsum float64 = 0
			var H_sum float64 = 0
			for j := 0; j < n; j++{
				Hsum += math.Abs(H[i][j])
				H_sum += math.Abs(H_[i][j])
			}
			if Hsum > HNorm {
				HNorm = Hsum
			}
			if H_sum > H_Norm {
				H_Norm = H_sum
			}
		}
		fmt.Fprintf(f,"Норма H = %g\n", HNorm)
		fmt.Fprintf(f,"Норма H_ = %g\n", H_Norm)
		fmt.Fprintf(f,"Число обусловленности = %g\n", HNorm * H_Norm)
	}
}
