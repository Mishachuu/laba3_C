#include <iostream>
#include "Exteption.h"
#include <locale.h>
using namespace std;


void EInvalidSize::Print() {
	cout << endl << "Ошибка:\nОперация не возможна. Размеры матриц не совпадают" << endl;
}


void EInvalidIndex::Print() {
	cout << endl << "Ошибка:\nОперация не возможна. Нет элемента расположенного по этим индексам " << endl;
}

void EInvalidMull::Print() {
	cout << endl << "Ошибка:\nЧисло столбцов матрицы A не равно числу строк матрицы B" << endl;
}


void EInvalidQuadrate::Print() {
	cout << endl << "Ошибка:\nМатрица не квадратная" << endl;
}


void EdivisionZero::Print() {
	cout << endl << "Ошибка:\nДеление на ноль" << endl;
}