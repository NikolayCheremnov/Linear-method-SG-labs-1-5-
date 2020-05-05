package main

import (
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
)

//программка генерирования системы линейных уравнений
//параметры задаются аргументами командной строки:
//число неизвестных, минимальное значение элемента матрицы Ab, максимальное значение матрицы Ab, имя файла, в который записать исходные данные
func main() {
	cmdArgs := os.Args
	n, err := strconv.Atoi(cmdArgs[1])
	if err != nil {
		log.Fatal(err)
	}
	min, err := strconv.Atoi(cmdArgs[2])
	if err != nil {
		log.Fatal(err)
	}
	max, err := strconv.Atoi(cmdArgs[3])
	if err != nil {
		log.Fatal(err)
	}
	fName := cmdArgs[4]
	f, err := os.Create(fName)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	fmt.Fprintln(f, n)
	for i := 0; i < n; i++{
		for j := 0; j <= n; j++ {
			fmt.Fprintf(f, "%d ", rand.Intn(int(math.Abs((float64(max-min))))) + min)
		}
		fmt.Fprintln(f)
	}
}
