__author__ = 'memery'

import xlwt
import pubmedgrabber
import os

def xlwt_explore():
    tester = xlwt.Workbook()
    test_sheet = tester.add_sheet('test')
    for n in range(10):
        test_sheet.write(n, 0, label=str(n))
    tester.save('tester.xls')


if __name__ == '__main__':
    xlwt_explore()

