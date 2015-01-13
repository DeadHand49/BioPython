from __future__ import unicode_literals
__author__ = 'memery'

from leadentry import Article, Author, lookup_up_title
from Bio import Entrez
from bs4 import BeautifulSoup
import unicodecsv as csv






if __name__ == '__main__':
    tester = ZoteroEntry(open('ZoteroTest.csv', 'rb'), 'Test')
    print tester.articles