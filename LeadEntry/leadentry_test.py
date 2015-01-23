# coding=utf-8
from __future__ import unicode_literals

__author__ = 'memery'

import leadentry
import unittest
from bs4 import BeautifulSoup, element


class LeadEntryTest(unittest.TestCase):

    def setUp(self):
        self.batch = leadentry.Batch('matthew.emery@stemcell.com')
        self.batch.pmids = ['24411336']
        self.batch.pubmed_xml = self.batch.create_pubmed_xml()
        self.batch.parse_pubmed_soup()

    def test_batch_title(self):
        self.assertEqual('A defined xeno-free and feeder-free culture system for the derivation, '
                         'expansion and direct differentiation of transgene-free patient-specific '
                         'induced pluripotent stem cells', self.batch.articles[0].info['Article Title'])

    def test_batch_abstract(self):
        self.assertIn('A defined xeno-free system for patient-specific iPSC derivation and differentiation',
                      self.batch.articles[0].info['Abstract'])

    def test_batch_pub_link(self):
        self.assertEqual('http://www.sciencedirect.com/science/article/pii/S0142961213015342',
                         self.batch.articles[0].info['Publication Link'])


class ArticleTest(unittest.TestCase):
    def setUp(self):
        tester = leadentry.Batch('matthew.emery@stemcell.com')
        with open('leadentry_test.xml') as xml:
            tester.pubmed_xml = BeautifulSoup(xml.read())
            self.test_article = leadentry.Article(info={'PMID': '25433608',
                                                        'Article Title': 'Mock'})
            tester.add_article(self.test_article)
            tester.parse_pubmed_soup()

    def test_find_pmid_in_xml(self):
        self.assertIsInstance(self.test_article.info['Tag'], element.Tag)

    def test_get_article(self):
        self.test_article.find_title()
        self.assertEqual('Mesenchymal stromal cells form vascular tubes when placed '
                         'in fibrin sealant and accelerate wound healing in vivo',
                         self.test_article.get_info('Article Title'))

    def test_get_date(self):
        self.assertEqual('11/26/2014', self.test_article.find_date())

if __name__ == "__main__":
    unittest.main()