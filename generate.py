import sys
import random


def gen(v_row_1, v_row_2, v_column_1, v_column_2, row, column):
    table = [[0 for i in range(row)] for j in range(column)]
    lst_1 = list()
    lst_2 = list()
    i = 0
    while True:
        if i == (len(v_row_1)):
            i = 0
        choice_1 = random.randint(1, column-1)
        choice_2 = random.randint(1, column-1)
        while True:
            if choice_1 in lst_1:
                choice_1 += 1
            elif choice_1 == (column - 1):
                choice_1 = 1
            else:
                break
        while True:
            if choice_2 in lst_2:
                choice_2 += 1
            elif choice_2 == (column - 1):
                choice_2 = 1
            else:
                break
        table[0][choice_1] = v_row_1[i]
        table[column - 1][choice_2] = v_row_1[i]
        lst_1.append(choice_1)
        lst_2.append(choice_2)
        i += 1
        if len(lst_1) == (column - 2):
            break
    i = 0
    lst_1 = list()
    lst_2 = list()
    while True:
        if i == (len(v_column_1)):
            i = 0
        choice_1 = random.randint(1, row-1)
        choice_2 = random.randint(1, row-1)
        while True:
            if choice_1 in lst_1:
                choice_1 += 1
            elif choice_1 == (row - 1):
                choice_1 = 1
            else:
                break
        while True:
            if choice_2 in lst_2:
                choice_2 += 1
            elif choice_2 == (row - 1):
                choice_2 = 1
            else:
                break
        table[choice_1][0] = v_column_1[i]
        table[choice_2][row - 1] = v_column_2[i]
        lst_1.append(choice_1)
        lst_2.append(choice_2)
        i += 1
        if len(lst_1) == (row - 2):
            break
    with open("table.txt", 'w') as f:
        f.write(str(row))
        f.write(' ')
        f.write(str(column))
        f.write('\n')
        for i in range(len(table)):
            line = str()
            for j in range(len(table[i])):
                line += str(table[i][j])
                line += ' '
            line = line[:-1]
            f.write(line)
            f.write('\n')


if __name__ == "__main__":
    n = int(sys.argv[1])
    m = int(sys.argv[2])
    gen([100], [100], [100], [100], n, m)
