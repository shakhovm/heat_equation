import random


def gen(v_row_1, v_row_2, v_column_1, v_column_2):
    table = [[0 for i in range(10)] for j in range(10)]
    for i in range(len(v_row_1)):
        choice_1 = random.randint(1,len(table[0])-1)
        choice_2 = random.randint(1,len(table[0])-1)
        table[0][choice_1] = v_row_1[i]
        table[len(table) - 1][choice_2] = v_row_2[i]
    for i in range(len(v_column_1)):
        choice_1 = random.randint(1,len(table)-1)
        choice_2 = random.randint(1,len(table)-1)
        table[choice_1][0] = v_column_1[i]
        table[choice_2][len(table) - 1] = v_column_2[i]
    with open("table.txt", 'w') as f:
        for i in range(len(table)):
            line = str()
            for j in range(len(table[i])):
                line += str(table[i][j])
                line += ' '
            line = line[:-1]
            f.write(line)
            f.write('\n')


if __name__ == "__main__":
    gen([40,100,1,8], [0,0,0,0], [0,8,5,54,0,0], [0,0,0,0,0,0])
