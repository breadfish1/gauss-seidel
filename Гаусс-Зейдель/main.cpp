#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <cmath>

#define E 0.0001
#define iterations 10
#define mon 5

double **create(int sz1, int sz2) { // создать двумерный массив
    double **t = (double**)malloc(sizeof(double*)*sz1);
    int i;
    
    for (i = 0; i < sz1; i++)
        t[i] = (double*)malloc(sizeof(double)*sz2);
    
    return t;
}

double *create1strDouble(int sz) { // выделить память под одномерный массив типа double и заполнить его нулями
    double *t = (double*)malloc(sizeof(double)*sz);
    
    for (int i = 0; i < sz; i++)
        t[i] = 0;
    
    return t;
}

int *create1strInt(int sz) { // выделить память под одномерный массив типа int и заполнить его натуральными числами с шагом 1, начиная с 0
    int *t = (int*)malloc(sizeof(int)*sz);
    
    for (int i = 0; i < sz; i++)
        t[i] = i;
    
    return t;
}

void del1str(double *x) { // очисить память от массива корней
    free(x);
}

void subPrint(double *ans, int sz) { // печать массива ответов
    int i;
    for (i = 0; i < sz - 1; i++)
        printf("X[%i] = %e, ", i, ans[i]);
    printf("X[%i] = %e\n", i, ans[i]);
}

double **reader(FILE *f, int *sz) { // считывает уравнения, возвращает указатель на считанную систему
    double **x = NULL;
    
    if ((f = fopen("eq.txt", "r")) != NULL) { // открываем файл на чтение, если он не пустой
        fscanf(f, "%i", sz);
        int sz1 = *sz + 1;
        x = create(*sz, sz1);
        int i = 0, j = 0;
        
        while (!feof(f)) {
            for (i = 0; i < *sz; i++)
                for (j = 0; j < sz1 && !feof(f); j++)
                    fscanf(f, "%lf", &(x[i][j]));
        }
        fclose(f); // закрываем
    }
    
    return x;
}

void del(double **x, int sz, int sz1) {
    for (int i = 0; i < sz; i++)
        free(x[i]);
    free(x);
}

void print(double **x, int *order, int sz1, int sz2) {
    for (int i = 0; i < sz1; i++) {
        for (int j = 0; j < sz2; j++)
            printf("%13e ", x[order[i]][j]);
        printf("\n");
    }
    printf("\n");
}


int zeros(double **x, int rows, int columns) { // проверка наличия нулей на диагонали, если нули есть, то возвращаем 0, иначе 1
    for (int i = 0; i < rows; i++)
        if (fabs(x[i][i]) < E)
            return 0;
    
    return 1;
}

int zerosOrder(double **x, int rows, int columns, int *order) { // проверка наличия нулей на диагонали в соответствии с порядком order, если нули есть, то возвращаем 0, иначе 1
    for (int i = 0; i < rows; i++)
        if (fabs(x[order[i]][i]) < E)
            return 0;
    
    return 1;
}

void copyDouble(double *curr, double *prev, int sz) { // копирует массив curr в prev (идентичные размерности)
    for (int i = 0; i < sz; i++)
        prev[i] = curr[i];
}

void copyInt(int *a, int *b, int sz) { // копирует массив a в b (идентичные размерности)
    for (int i = 0; i < sz; i++)
        b[i] = a[i];
}

void shift(int *str, int start, int end) { // сдвинуть массив str влево в пределах от start до end
    int temp = str[start];
    for (int i = start; i < end - 1; i++)
        str[i] = str[i + 1];
    
    str[end - 1] = temp;
}

double iter(double **x, double *curr, int *order, int rows, int columns) { // совершить одну итерацию
    double temp, temp1, max = 0.0;
    
    for (int i = 0; i < rows; i++) {
        temp = x[order[i]][columns - 1]; // свободный член
        
        for (int j = 0; j < i; j++) // отнимаем от свободного члена до искомого
            temp -= x[order[i]][j] * curr[j];
        for (int j = i + 1; j < columns - 1; j++) // после искомого
            temp -= x[order[i]][j] * curr[j];
        
        temp = temp / x[order[i]][i]; // считаем очередной Xi
        temp1 = fabs(temp - curr[i]); // модуль разности между предыдущей и текущей итерациями
        max = temp1 > max ? temp1 : max; // смотрим, не больше ли текущий модуль предыдущего
        
        curr[i] = temp;
    }
    
    return max;
}

double *gaussZeidel(double **x, int rows, int columns, int *order) { // решение методом Гаусса-Зейделя без проверки монотонности, возвращает массив корней размером rows
    double max;
    double *curr = create1strDouble(rows);
    
    do {
        max = iter(x, curr, order, rows, columns);
    } while (max > E);
    
    return curr;
}

double *gaussZeidelMon(double **x, int rows, int columns, int *order) { // решение методом Гаусса-Зейделя с проверкой монотоноости
    double max = 0.0, max_prev = 0.0;
    double *curr = create1strDouble(rows);
    int c = 0;
    
    for (int i = 0; i < iterations; i++) { // выполняем iterations итераций
        max = iter(x, curr, order, rows, columns); // находим максималный модуль разницы iой итерации
        if (max < max_prev) // сравниваем с предыдущим, если текущее меньше, то монотонность сохраняется, иначе обнуляем
            c++;
        else
            c = 0;
        max_prev = max;
    }
    
    del1str(curr);
    
    if (c >= mon) // если необходимое количество итераций монотонность сохранилась, то вызываем решение без проверки монотонности
        return gaussZeidel(x, rows, columns, order);
    else		  // иначе возвращаем NULL, что означает, что данным методом решить СЛАУ невозможно
        return NULL;
}

int DUS(double **x, int rows, int columns, double *amount, int *order) { // проверка ДУС, возвращает 1, если проверка пройдена, 0, если не пройдена
    double temp;
    int i, ind = 0;
    for (i = 0; i < rows; i++) {
        temp = 0;
        
        temp = amount[order[i]] - fabs(x[order[i]][i]); // отнимаем модуль диагонального элемента от суммы
        
        if (fabs(x[order[i]][i]) > fabs(temp) && fabs(fabs(temp) - fabs(x[order[i]][i])) > E) // сравниваем его с получившимся результатом, если диагональный элемент больше суммы, то ind = 1
            ind = 1;
        else
            if(fabs(fabs(temp) - fabs(x[order[i]][i])) > E) // модуль диагонального элемента не равен модулю суммы остальных элементов
                return 0;
    }
    
    return ind; // вернется 1, если хотя бы 1 диагональный элемент больше суммы всей строки, а остальные хотя бы равны. 0, если все равны
}

/* перестановка строк, избавление от 0 на диагонали, проверка ДУС, если проверка пройдена, решаем гауссом-зейделем, иначе с проверкой монотонности
 ind = 1, если найдена перестановка с ДУС, 2, если без ДУС, -1, если не найдена */
void zerosChanging(double **x, int rows, int *order, int *temp,  double *amount, int start, int *ind) {
    if (*ind != 1) { // если не нашли перестановку с ДУС
        for (int i = 0; i < rows - start; i++) {
            if (start < rows - 2)
                zerosChanging(x, rows, order, temp, amount, start + 1, ind);
            else {
                if (zerosOrder(x, rows, rows + 1, order)) { // зайдем, если нулей на диагонали с перестановкой order нет
                    if (DUS(x, rows, rows + 1, amount, order)) { // если нашли перестановку с ДУС, то сохраняем ее
                        copyInt(order, temp, rows);
                        *ind = 1;
                        return; // выходим из функции
                    }
                    
                    if (*ind != 2) { // если до этого такой перестановки не было найдено, то сохраняем ее в temp
                        copyInt(order, temp, rows);
                        *ind = 2;
                    }
                }
            }
            shift(order, start, rows);
        }
    }
}

double *AM(double **x, int rows, int columns) { // посчитать сумму строк
    double *a = create1strDouble(rows);
    
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < columns - 1; j++)
            a[i] += fabs(x[i][j]);
    
    return a;
}

double *solving(double **x, int rows, int columns) { // основная функция решения СЛАУ, возвращает решение системы в виде массива double
    double *amount = AM(x, rows, columns); // сумма строк
    int *order = create1strInt(rows), *temp = create1strInt(rows), ind = -1;
    if (zeros(x, rows, columns)) { // если нулей нет
        if (DUS(x, rows, columns, amount, order)) // вызываем функцию проверки ДУС, если прошли проверку, то решаем методом Гаусса-Зейделя
            return gaussZeidel(x, rows, columns, order);
        else
            return gaussZeidelMon(x, rows, columns, order); // иначе решаем с контролем монотонности
    }
    else {
        zerosChanging(x, rows, order, temp, amount, 0, &ind); // вызов функции перестановки для избавления от нулей с выполнением ДУС
        
        print(x, temp, rows, columns);
        
        if (ind == 1) // найдена перестановка с ДУС
            return gaussZeidel(x, rows, columns, temp);
        
        if (ind == 2) // найдена без ДУС
            return gaussZeidelMon(x, rows, columns, temp); 
        
        if (ind == -1) { // невозможно избавиться от нулей на диагонали
            return NULL;
        }
    }
    
    return NULL;
}

int main() {
    double **x, *ans;
    int sz, sz1;
    FILE *f = NULL;
    
    x = reader(f, &sz);
    sz1 = sz + 1;
    
    int *order = create1strInt(sz);
    
    print(x, order, sz, sz1);
    
    ans = solving(x, sz, sz1);
    
    if (ans) {
        subPrint(ans, sz);
        del1str(ans);
    }
    else
        printf("There is no solving\n");
    
    del(x, sz, sz1);
    
    system("pause");
}
