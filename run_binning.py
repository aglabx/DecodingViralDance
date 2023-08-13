import os
from stack_frames import *


def calculate(file):
    print(file)
    f = read_and_bin(file)
    return f


def main():
    # Schedule three calls *concurrently*:
    list_of_files = os.listdir('/home/daria/Downloads/test/draft/')
    os.chdir('/home/daria/Downloads/test/draft/')
    L = [calculate(file) for file in list_of_files]
    df = ft.reduce(lambda left, right: pd.merge(left, right, on=['x', 'y'], how='outer'), L)
    df.to_csv('summary_table.csv')


main()
