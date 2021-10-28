**Данная программа принимает 6 параметров:**<br>
x_init - x координата начальной точки алгоритма<br>
y_init - y координата начальной точки алгоритма<br>
alpha - коэфициэнт алгоритма<br>
beta - коэфициэнт алгоритма<br>
gamma - коэфициэнт алгоритма<br>
sigma - коэфициэнт алгоритма<br>

**Вывод программы** - три числа. <br>
Первые 2 - координаты найденнтой точки минимума, 3 - значение функции в данной точке.

Следующие запуски отвечают на вопросы задания:

**Запуск 1:**
<br>
./program.out 10 10 1.2 0.3 1.2 0.1

Вывод:
<br>
9.43613 10.9842 -106.765

**Запуск 2:**
<br>
./a.out 5 5 1 0.5 2 0.5 

Вывод:
<br>
3.15294 4.70104 -106.765

**Запуск 3:**
<br>
./a.out 0 0 1 0.5 2 0.5

Вывод:
<br>
0.905521 0.665279 1.48702

**Запуск 4:**
<br>
./a.out 0 0 1 0.5 2 0.5

Вывод:
<br>
0.905521 0.665279 1.48702
            
**Запуск 5:**
<br>
./a.out 0 0 10 0.5 2 0.5

Вывод:
<br>
3.97471 3.66754 -14.5609


Запуски 1, 2, 3 иллюстрируют, что при разных начальных точка метод может сходиться к разным точкам.<br>
Запуски 4, 5 иллюстрируют, что для одной и той же начальной точки метод может сойтись в разные конечные точки при различных наборах гиперпараметров.
